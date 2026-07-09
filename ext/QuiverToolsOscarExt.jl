module QuiverToolsOscarExt

# Oscar-backed walls-and-chambers / VGIT functionality. Loaded automatically when a
# user runs `using Oscar` alongside QuiverTools. The docstrings for the exported
# entry points live on the stub definitions in src/WallsAndChambers.jl.

using QuiverTools
using Oscar
using Singular: Singular
using Memoization: Memoization
using IterTools: IterTools

import Memoization: @memoize
import QuiverTools:
  is_special_subdimension_vector, all_special_subdimension_vectors, sst,
  vgit_walls, wall_system, vgit_chambers, vgit_fan, git_equivalent

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
Only works for strongly convex polyhedra.
"""
function __helper_accelerate(P::Oscar.Polyhedron)
  newP = Oscar.polyhedron(Oscar.positive_hull(Oscar.rays(P)))
  newP != P && return P
  return newP
end

@memoize Dict function sst(Q, e)
  e_perp = Oscar.polyhedron([e, -e], [0, 0]) #e^{\perp}
  all_gen = filter(
    eprime -> !all(ei == 0 for ei in eprime) && eprime != e,
    all_general_subdimension_vectors(Q, e),
  )
  isempty(all_gen) && return e_perp
  return intersect(e_perp, Oscar.polyhedron(all_gen, zeros(Int, length(all_gen))))
end

@memoize Dict function vgit_walls(Q, d; inner=false, top_dimension=true)
  all_subd = QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true)
  out = map(
    e -> reduce(intersect, [sst(Q, e), sst(Q, d - e), sst(Q, d)]), # definition of W_{e}
    all_subd)
  top_dim = maximum(Oscar.dim(w) for w in out) # max dimension of WALLS
  top_dimension && filter!(w -> Oscar.dim(w) == top_dim, out)
  inner && filter!(
    w -> !any(issubset(w, f) for f in Oscar.facets(Oscar.Polyhedron, sst(Q, d))), out
  )
  return unique!(__helper_accelerate, out)
end

function wall_system(Q, d; inner=false, as_cones=true)
  sstd = sst(Q, d)

  all_walls = vgit_walls(Q, d; top_dimension=true)
  if inner
    all_walls = filter(
      w -> !any(issubset(w, f) for f in Oscar.facets(Oscar.Polyhedron, sstd)), all_walls
    )
  end
  # each affine hull contains the defining hyperplanes that cut it out.
  all_walls = map(w -> Oscar.affine_hull(w), all_walls)
  all_walls = map(
    AH -> Oscar.polyhedron(
      vcat([hyp.a[1, :] for hyp in AH], [-hyp.a[1, :] for hyp in AH]),
      zeros(Int, 2 * length(AH)),
    ),
    all_walls,
  )
  as_cones && return map(w -> intersect(w, sstd), all_walls)
  return all_walls
end

@memoize Dict function vgit_chambers(Q, d; verbose=false)
  sstd = __helper_accelerate(sst(Q, d))

  sstd_dim = Oscar.dim(sstd)
  # top-dimensional inner walls
  int_walls = vgit_walls(Q, d; top_dimension=true)
  int_walls = filter(
    w -> !any(issubset(w, f) for f in Oscar.facets(Oscar.Polyhedron, sstd)), int_walls
  )

  # @warn "we must treat the case of a wall W_e of dimension = Oscar.dim(sstd) here!"

  # we split the walls into two sets: the ones equal to the wall system hyperplane
  # they lay on, and the ones that are not.
  wallsyst = wall_system(Q, d; inner=true, as_cones=true)

  full_walls = filter(w -> w in wallsyst, int_walls)
  smaller_walls = filter(w -> !(w in wallsyst), int_walls)

  full_walls = map(w -> Oscar.affine_hull(w), full_walls)
  smaller_walls_ah = map(w -> Oscar.affine_hull(w), smaller_walls)

  # find a vector not equal to d that cuts out the wall
  full_walls = map(AH -> AH[findfirst(hyp -> (hyp.a[1, :] != d), AH)].a[1, :], full_walls)
  smaller_walls_ah = map(
    AH -> AH[findfirst(hyp -> (hyp.a[1, :] != d), AH)].a[1, :], smaller_walls_ah
  )

  full_walls_iterate = eachindex(full_walls)
  smaller_walls_iterate = eachindex(smaller_walls_ah)

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

  verbose && @info "Treating the full walls..."
  for i in full_walls_iterate
    if verbose
      @time begin
        # in-place
        n = length(top_chambers)
        for j in 1:n
          push!(top_chambers,
            helper_split(full_walls[i], top_chambers[j])...,
          )
        end
        deleteat!(top_chambers, 1:n)

        map!(__helper_accelerate, top_chambers, top_chambers)
        filter!(c -> Oscar.dim(c) == sstd_dim, top_chambers)
      end
    else
      # in-place
      n = length(top_chambers)
      for j in 1:n
        push!(top_chambers,
          helper_split(full_walls[i], top_chambers[j])...,
        )
      end
      deleteat!(top_chambers, 1:n)

      map!(__helper_accelerate, top_chambers, top_chambers)
      filter!(c -> Oscar.dim(c) == sstd_dim, top_chambers)
    end
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  length(smaller_walls) == 0 && return top_chambers

  verbose && @info "Treating the smaller walls..."
  for i in smaller_walls_iterate
    if verbose
      @time begin
        top_chambers = vcat(
          map(
            chamber -> if issubset(smaller_walls[i], chamber)
              helper_split(smaller_walls_ah[i], chamber)
            else
              [chamber]
            end,
            top_chambers,
          )...,
        )
        map!(__helper_accelerate, top_chambers, top_chambers)
        filter!(c -> Oscar.dim(c) == length(d) - 1, top_chambers)
      end
    else
      top_chambers = vcat(
        map(
          chamber -> if issubset(smaller_walls[i], chamber)
            helper_split(smaller_walls_ah[i], chamber)
          else
            [chamber]
          end,
          top_chambers,
        )...,
      )
      map!(__helper_accelerate, top_chambers, top_chambers)
      filter!(c -> Oscar.dim(c) == sstd_dim, top_chambers)
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
  return Oscar.polyhedral_fan(
    map(ch -> Oscar.positive_hull(Oscar.rays(ch)),
      vgit_chambers(Q, d; verbose=verbose),
    ),
  )
end

function git_equivalent(Q, d, theta1, theta2)
  theta1 == theta2 && return true
  line = Oscar.convex_hull([theta1, theta2]) # 1-dimensional iif theta1 != theta2

  # either the line lies in a wall or it intersects none of them
  for w in vgit_walls(Q, d; top_dimension=false)
    if !Oscar.issubset(line, w) && Oscar.is_feasible(Oscar.intersect(line, w))
      return false
    end
  end
  return true
end

"""
    __sst_cone(Q::Quiver, e::AbstractVector{Int})

Return the exact same polyhedral cone as `sst(Q, e)`
but seen as an Oscar `Cone` object.

For internal use only for now.
"""
function __sst_cone(Q::Quiver, e::AbstractVector{Int})
  all_gen = filter(
    eprime -> !all(ei == 0 for ei in eprime) && eprime != e,
    all_general_subdimension_vectors(Q, e),
  )
  isempty(all_gen) && return Oscar.cone_from_inequalities([e, -e])
  return Oscar.cone_from_inequalities(all_gen, [e])
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
  walls = map(
    e -> Oscar.intersect(__sst_cone(Q, e), __sst_cone(Q, d - e)),
    QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  isempty(walls) && return Oscar.polyhedral_fan(__sst_cone(Q, d)) # I guess
  return Oscar.polyhedral_fan(walls)
end

# NB: no `@compile_workload` here. The VGIT functions return polymake C++ objects
# (via Oscar's polyhedral geometry) which cannot be serialized into a precompile
# image: cached objects come back "deleted". Precompiling this path would poison
# the Memoization caches with dead objects, so it is left to compile on first use.

end # module QuiverToolsOscarExt
