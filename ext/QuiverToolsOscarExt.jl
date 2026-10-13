module QuiverToolsOscarExt

# Oscar-backed walls-and-chambers / VGIT functionality. Loaded automatically when a
# user runs `using Oscar` alongside QuiverTools. The docstrings for the exported
# entry points live on the stub definitions in src/WallsAndChambers.jl.
# To define new functionality, add it here, add a stub in src/WallsAndChambers.jl,
# import the stub here, add it to the list at the bottom of src/WallsAndChambers.jl,
# and export it from QuiverTools.jl.

using QuiverTools
using Oscar
using Singular: Singular
using Memoization: Memoization
using IterTools: IterTools

import Memoization: @memoize
import QuiverTools:
  is_special_subdimension_vector, all_special_subdimension_vectors, sst,
  vgit_walls, wall_system, vgit_chambers, vgit_fan, all_stability_parameters

# Disambiguate between Singular's and Oscar's overloads of `(::PolyRing)(::spoly)`,
# which collide once both packages are loaded.
(R::Singular.PolyRing)(f::Singular.spoly) = invoke(
  R, Tuple{Union{Singular.spoly,Oscar.AbstractAlgebra.MPolyRingElem}}, f
)

function is_special_subdimension_vector(
  Q::Quiver, e::AbstractVector{Int}, d::AbstractVector{Int}
)
  e = QuiverTools.coerce_vector(e)
  d = QuiverTools.coerce_vector(d)

  (all(ei -> ei == 0, e) || e == d) && return false
  is_general_subdimension_vector(Q, e, d) && return false

  eprimes = filter(
    eprime -> !all(epi -> epi == 0, eprime),
    all_general_subdimension_vectors(Q, e),
  )

  fprimes = filter(
    fprime -> !all(fpi -> fpi == 0, fprime),
    all_general_subdimension_vectors(Q, d - e),
  )

  cone = Oscar.polyhedron(
    vcat(eprimes, [e + f for f in fprimes]),
    zeros(Int, length(eprimes) + length(fprimes)),
  )
  cone = intersect(cone, sst(Q, d))

  return !any(
    issubset(cone, Oscar.polyhedron([eprime, -eprime], [0, 0]))
    for eprime in eprimes
  ) &&
         !any(
    issubset(cone, Oscar.polyhedron([e + f, - (e + f)], [0, 0]))
    for f in fprimes if f != d - e # this would check if cone ⋐ d^⟂, which is true of course
  )
end

@memoize Dict function all_special_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})
  d = QuiverTools.coerce_vector(d)
  candidates = filter(
    e -> !QuiverTools.is_general_subdimension_vector(Q, e, d),
    QuiverTools.all_subdimension_vectors(d, nonzero=true, strict=true),
  )
  return filter(
    e -> is_special_subdimension_vector(Q, e, d),
    candidates,
  )
end

"""
    __helper_accelerate(P::Oscar.Polyhedron)

Internal method.

Convert the polyhedron `P` to be spanned by its rays.
Polyhedra with a nontrivial lineality space or without rays are returned unchanged.
"""
function __helper_accelerate(P::Oscar.Polyhedron)
  # rays(P) omits the lineality space, so the positive hull would not recover P
  (Oscar.lineality_dim(P) > 0 || isempty(Oscar.rays(P))) && return P
  newP = Oscar.polyhedron(Oscar.positive_hull(Oscar.rays(P)))
  newP != P && return P
  return newP
end

@memoize Dict function sst(Q, e)
  data = QuiverTools.__semistable_cone_data(Q, e)
  equations = vcat(data.equations, .-data.equations)
  e_perp = Oscar.polyhedron(equations, zeros(Int, length(equations)))
  isempty(data.inequalities) && return e_perp
  return intersect(
    e_perp,
    Oscar.polyhedron(data.inequalities, zeros(Int, length(data.inequalities))),
  )
end

function vgit_walls(Q, d; inner=false, top_dimension=true)
  out = __vgit_walls(Q, d)
  if top_dimension
    dims = map(w -> Oscar.dim(w), out)
    top_dim = maximum(dims) # max dimension of WALLS
    out = out[findall(dims .== top_dim)]
  end
  if inner
    facettes = Oscar.facets(Oscar.Polyhedron, sst(Q, d))
    filter!(
      w -> !any(issubset(w, f) for f in facettes),
      out,
    )
  end
  return out
end

@memoize Dict function __vgit_walls(Q, d)
  walls_data = map(
    e -> QuiverTools.__vgit_wall_data(Q, d, e),
    QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  out = map(
    wall_datum ->
      Oscar.polyhedron(
        vcat(wall_datum.equations, .-wall_datum.equations, wall_datum.inequalities),
        zeros(Int, 2 * length(wall_datum.equations) + length(wall_datum.inequalities)),
      ),
    walls_data)
  return unique!(__helper_accelerate, out)
end

function wall_system(Q, d; inner=false, as_cones=true)
  sstd = sst(Q, d)

  all_walls = map(vgit_walls(Q, d; inner=inner, top_dimension=true)) do w
    Oscar.polyhedron(nothing, Oscar.affine_hull(w))   # only pass to the polyhedron constructor the affine hull's equations
  end
  as_cones && return map(w -> intersect(w, sstd), all_walls)
  return unique!(all_walls)
end

@memoize Dict function vgit_chambers(Q, d; verbose=false)
  sstd = __helper_accelerate(sst(Q, d))

  sstd_dim = Oscar.dim(sstd)
  # top-dimensional inner walls
  int_walls = vgit_walls(Q, d; inner=true, top_dimension=true)

  # @warn "we must treat the case of a wall W_e of dimension = Oscar.dim(sstd) here!"

  # we split the walls into two sets: the ones equal to the wall system hyperplane
  # they lay on, and the ones that are not.
  wallsyst = wall_system(Q, d; inner=true, as_cones=true)

  full_walls = filter(w -> w in wallsyst, int_walls)
  smaller_walls = filter(w -> !(w in wallsyst), int_walls)

  # n is a multiple of d iff all 2x2 minors of the matrix with rows n and d vanish
  is_multiple_of_d(n) =
    all(n[i] * d[j] == n[j] * d[i] for i in eachindex(d), j in eachindex(d))

  # a normal vector not proportional to d of a hyperplane in the affine hull of the wall;
  # together with d^⟂ it cuts out the hyperplane of d^⟂ spanned by the wall
  function wall_normal(w)
    normals = map(Oscar.normal_vector, Oscar.affine_hull(w))
    k = findfirst(!is_multiple_of_d, normals)
    k === nothing && error(
      "One wall spans d^⟂ for d = $d, so it does not cut out a hyperplane of d^⟂"
    )
    return normals[k]
  end
  full_normals = map(wall_normal, full_walls)
  smaller_normals = map(wall_normal, smaller_walls)

  top_chambers = [sstd]

  # helper function
  function helper_split(wall, chamber)
    return [
      intersect(chamber, Oscar.polyhedron([wall], [0])),
      intersect(chamber, Oscar.polyhedron([-wall], [0])),
    ]
  end
  verbose &&
    @info "There are $(length(full_walls)) full walls and $(length(smaller_walls)) smaller walls."

  # run `f`, timing it when verbose, so that both modes execute the same code
  timed(f) = verbose ? (@time f()) : f()

  verbose && @info "Treating the full walls..."
  for i in eachindex(full_normals)
    top_chambers = timed() do
      chambers = vcat(
        map(chamber -> helper_split(full_normals[i], chamber), top_chambers)...
      )
      map!(__helper_accelerate, chambers, chambers)
      return filter!(c -> Oscar.dim(c) == sstd_dim, chambers)
    end
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  length(smaller_walls) == 0 && return top_chambers

  verbose && @info "Treating the smaller walls..."
  for i in eachindex(smaller_walls)
    top_chambers = timed() do
      chambers = vcat(
        map(
          chamber -> if issubset(smaller_walls[i], chamber)
            helper_split(smaller_normals[i], chamber)
          else
            [chamber]
          end,
          top_chambers,
        )...,
      )
      map!(__helper_accelerate, chambers, chambers)
      return filter!(c -> Oscar.dim(c) == sstd_dim, chambers)
    end
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  verbose && @info "Removing fake walls..."

  # remove fake walls that appeared during the previous loop
  incomplete = true
  while incomplete
    incomplete = false
    for (ch1, ch2) in IterTools.subsets(top_chambers, 2)
      inters = intersect(ch1, ch2)
      # we are looking for common facets
      Oscar.dim(inters) != sstd_dim - 1 && continue
      # if their intersection (the common facet) lays on a W_e, good.
      any(issubset(inters, w) for w in int_walls) && continue
      verbose && @info "Found a fake wall, removing it..."
      incomplete = true
      # otherwise, we replace the two chambers in top_chambers by their union
      push!(top_chambers, Oscar.minkowski_sum(ch1, ch2))
      top_chambers = filter(ch -> (ch != ch1 && ch != ch2), top_chambers)
    end
  end
  verbose &&
    @info "Done. Found $(length(top_chambers)) chambers after removing the fake walls."
  return top_chambers
end

function vgit_fan(Q, d; verbose=false)
  chambers = vgit_chambers(Q, d; verbose=verbose)
  !all(Oscar.lineality_dim(ch) == 0 for ch in chambers) && throw(
    ArgumentError(
      """The VGIT fan is not strongly convex.
      Unsupported for now, see https://github.com/QuiverTools/QuiverTools.jl/issues/40""",
    ),
  )
  return Oscar.polyhedral_fan(
    map(ch -> Oscar.positive_hull(Oscar.rays(ch)),
      vgit_chambers(Q, d; verbose=verbose),
    ),
  )
end

"""
    __sst_cone(Q::Quiver, e::AbstractVector{Int})

Return the exact same polyhedral cone as `sst(Q, e)`
but seen as an Oscar `Cone` object.

For internal use only for now.
"""
function __sst_cone(Q::Quiver, e::AbstractVector{Int})
  data = QuiverTools.__semistable_cone_data(Q, e)
  isempty(data.inequalities) &&
    return Oscar.cone_from_inequalities(vcat(data.equations, .-data.equations))
  return Oscar.cone_from_inequalities(data.inequalities, data.equations)
end

"""
    __walls_cones(Q::Quiver, d::AbstractVector{Int})

Compute the polyhedral cones of the walls `W_e` of the quiver `Q` with dimension vector `d`.

For internal use only for now.
"""
function __walls_cones(Q::Quiver, d::AbstractVector{Int})
  out = map(
    e -> Oscar.intersect(__sst_cone(Q, e), __sst_cone(Q, d - e)),
    QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  return unique!(out)
end

"""
    __lower_fan(Q::Quiver, d::AbstractVector{Int})

Compute the polyhedral fan from the definition
using the fan construction via cones.

Does not contain the top dimensional chambers,
so it cannot be used to iterate over chambers.

It is way faster and more likely to be correct in G.P. opinion.

For internal use only for now.
"""
function __lower_fan(Q::Quiver, d::AbstractVector{Int})
  walls = __walls_cones(Q, d)
  isempty(walls) && return Oscar.polyhedral_fan(__sst_cone(Q, d)) # I guess
  return Oscar.polyhedral_fan(walls)
end

"""
    general_stability(cone)

Given a polyhedral cone of the VGIT fan,
return a stability parameter in its relative interior.

In practice this is a rescaled sum of the rays of the polyhedral cone input.

# Example

```jldoctests
julia> Q = Quiver("1-2---3,1--3"); d = [1, 1, 1];

julia> ch = vgit_chambers(Q, d); map(general_stability, ch)
2-element Vector{Vector{Int64}}:
 [2, -1, -1]
 [1, 1, -2]
```

This function allows to sample stability parameters from every equivalence class:

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> walls = vgit_walls(Q, d; top_dimension=false); map(general_stability, walls)
7-element Vector{Vector{Int64}}:
 [0, 1, 0, -1]
 [1, 0, 0, -1]
 [1, -1, 1, -1]
 [2, -1, 0, -1]
 [1, 0, -1, 0]
 [1, 1, -1, -1]
 [1, 0, -1, 0]
```
"""
function general_stability(cone)
  out = sum(Oscar.rays(cone))

  c = lcm(denominator.(out)...)
  out .*= c

  l = gcd(out...)
  return Int.(out ./ l)
end

function all_stability_parameters(Q::Quiver, d::AbstractVector{Int}; generic::Bool=false)
  F = vgit_fan(Q, d)
  if generic
    cones = Oscar.cones(F, Oscar.dim(F))
  else
    cones = vcat([collect(Oscar.cones(F, i)) for i in 1:Oscar.dim(F)]...)
  end
  return map(general_stability, cones) #do we add the semisimple condition?
end

# NB: no `@compile_workload` here. The VGIT functions return polymake C++ objects
# (via Oscar's polyhedral geometry) which cannot be serialized into a precompile
# image: cached objects come back "deleted". Precompiling this path would poison
# the Memoization caches with dead objects, so it is left to compile on first use.

end # module QuiverToolsOscarExt
