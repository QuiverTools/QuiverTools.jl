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
  vgit_walls, wall_system, vgit_chambers, vgit_fan, all_stability_parameters,
  project_to_d_perp

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
  isempty(out) && return out
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
  if gcd(d) != 1 || !is_schur_root(Q, d)
    throw(
      ArgumentError(
        "The dimension vector $d is not a primitive Schur root for the quiver $Q. This case is unsupported for now, see Issue #53. Also see Issue #40 and [Theorem 5.1, https://doi.org/10.1016/j.jpaa.2025.108153]."
      ),
    )
  end
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
  as_cones && (all_walls = map(w -> intersect(w, sstd), all_walls))
  return unique!(all_walls)
end

@memoize Dict function vgit_chambers(Q, d; verbose=false)
  if gcd(d) != 1 || !is_schur_root(Q, d)
    throw(
      ArgumentError(
        "The dimension vector $d is not a primitive Schur root for the quiver $Q. This case is unsupported for now, see Issue #53. Also see Issue #40 and [Theorem 5.1, https://doi.org/10.1016/j.jpaa.2025.108153]."
      ),
    )
  end
  sstd = __helper_accelerate(sst(Q, d))

  sstd_dim = Oscar.dim(sstd)
  # top-dimensional inner walls
  int_walls = vgit_walls(Q, d; inner=true, top_dimension=true)

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
      "There is at least one wall that spans d^⟂ for d = $d."
    )
    # under the sanity check above,
    # this can never happen, see  [Theorem 5.1, https://doi.org/10.1016/j.jpaa.2025.108153].
    return normals[k]
  end
  full_normals = map(wall_normal, full_walls)
  smaller_normals = map(wall_normal, smaller_walls)

  top_chambers = [sstd]

  # splits a top dimensional cone in two along an hyperplane
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
    unique!(top_chambers)
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  length(smaller_walls) == 0 && return top_chambers

  verbose && @info "Treating the smaller walls..."
  for i in eachindex(smaller_walls)
    top_chambers = timed() do
      chambers = vcat(
        map(
          # split when the wall meets the chamber in codimension 1
          chamber -> if Oscar.dim(intersect(smaller_walls[i], chamber)) == sstd_dim - 1
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
    unique!(top_chambers)
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  verbose && @info "Removing fake walls..."

  # Scan the list to find the chambers that must be joined
  # and classify them with their indices.
  # Indices are stored in `parent` and updated with root(),
  # apparently a "union find data structure".
  n = length(top_chambers)
  parent = collect(1:n)
  root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))

  timed() do
    for (i, j) in IterTools.subsets(1:n, 2)
      root(i) == root(j) && continue # already known to be in the same chamber

      inters = intersect(top_chambers[i], top_chambers[j])
      # we are looking for common facets
      Oscar.dim(inters) != sstd_dim - 1 && continue

      # If the common facet lies on a wall W_e we don't merge.
      # Note: the common facet is either in the interior
      # of a "true" chamber, or it is contained in the intersection C1 \cap C2 of two "true" chambers.
      # Such an intersection is contained in a unique wall W_e for some e,
      # so the error that commit 32179e8b374fa50bf13126a2413c6dbc6b6e5d39
      # "fixed" could not happen.
      # We keep this check because it is faster.
      p = Oscar.relative_interior_point(inters)
      any(w -> p in w, int_walls) && continue

      verbose && @info "Found a fake wall."
      parent[root(i)] = root(j)
    end
  end

  verbose && @info "Merging the fake walls now:"
  groups = [findall(i -> root(i) == r, 1:n) for r in unique(map(root, 1:n))]
  timed() do
    top_chambers = map(g -> reduce(Oscar.minkowski_sum, top_chambers[g]), groups)
  end
  verbose &&
    @info "Done. Found $(length(top_chambers)) chambers after removing the fake walls."
  return top_chambers
end

function vgit_fan(Q, d; verbose=false)
  chambers = vgit_chambers(Q, d; verbose=verbose)
  # Chambers can have a nontrivial lineality space, e.g. when d has a zero entry.
  # `Oscar.rays` would then be empty, so we must rebuild each chamber as cone(R) + span(L).
  return Oscar.polyhedral_fan(
    map(chambers) do ch
      R, L = Oscar.rays_modulo_lineality(ch)
      # if the chamber is a vector space, add the zero vector as its only ray
      isempty(R) && (R = zeros(Int, 1, Oscar.ambient_dim(ch)))
      Oscar.positive_hull(R, L)
    end,
  )
end

# internal method that samples a stability parameter
# from a cone without lineality and rescales it to be integer and coprime.
function __general_stability(cone)
  out = Oscar.relative_interior_point(cone)

  c = lcm(denominator.(out)...)
  out .*= c
  l = gcd(out...)
  out = Int.(out ./ l)
  return out
end

# The GIT equivalence classes are exactly the relative interiors
# of each cone in the fan.
# Note on the class of 0: a priori, one can show that
# the equivalence class of 0 is
# Z = {θ : θᵢ = 0 for every i with dᵢ > 0}.
# By knowing that the lineality space of the fan is its smallest
# cone, and relint(V) = V for a vector space V,
# we conclude that Z = the entire lineality space.
# A corollary is that a honest primitive Schur root
#  cannot have lineality space in the GIT fan.
function all_stability_parameters(
  Q::Quiver, d::AbstractVector{Int}; generic::Bool=false, nonzero::Bool=true
)
  F = vgit_fan(Q, d)
  # The smallest cone of the fan is its lineality space L. It always contains
  # Z = {θ : θ_i = 0 whenever d_i > 0}, which is the GIT class of θ = 0;
  l = Oscar.lineality_dim(F)
  # should never happen by the explanation above, but just in case
  l == count(iszero, d) || error("unexpected lineality space of dimension $l for d = $d")

  # dimensions on which we iterate
  dims = (l + 1):Oscar.dim(F)
  generic && (dims = filter(==(Oscar.dim(F)), dims)) # might be empty, can't just =[Oscar.dim(F)] this
  # the empty vector functions as init=[]
  out = vcat(Vector{Int}[], [map(__general_stability, Oscar.cones(F, i)) for i in dims]...)
  # the class of θ = 0 is L; it is a chamber only if L is the whole fan
  if !nonzero && (!generic || Oscar.dim(F) == l)
    pushfirst!(out, zeros(Int, length(d)))
  end
  return out
end

########################################################################
# Below lie functions to project all of the objects above to d^{\perp}.
# This is for user-friendliness and to allow Oscar.visualize()
# to display 3-dimensional fans in \mathbb{R}^4.
########################################################################

# The coordinates kept by `project_to_d_perp`: all but `drop`, which needs d[drop] != 0
# for the projection to be injective on d^⟂.
function __d_perp_coordinates(d::AbstractVector{Int}, drop::Int)
  drop in eachindex(d) && d[drop] != 0 || throw(
    ArgumentError("cannot drop coordinate $drop of d = $d: it must be a nonzero entry")
  )
  return [j for j in eachindex(d) if j != drop]
end

function __check_in_d_perp(generators, d::AbstractVector{Int}, X)
  for v in generators
    sum(d[j] * v[j] for j in eachindex(d)) == 0 ||
      throw(ArgumentError("$X does not lie in the hyperplane orthogonal to d = $d"))
  end
end

# the generators as the rows of a matrix, keeping only the coordinates `keep`
__rows(generators, keep) = [v[j] for v in generators, j in keep]

function project_to_d_perp(
  P::Oscar.Polyhedron, d::AbstractVector{Int}; drop::Int=findlast(!=(0), d)
)
  keep = __d_perp_coordinates(d, drop)
  K = Oscar.coefficient_field(P)
  # an infeasible inequality: the empty polyhedron
  Oscar.is_feasible(P) || return Oscar.polyhedron(K, zeros(Int, 1, length(keep)), [-1])
  # P = conv(V) + cone(R) + span(L); `vertices` and `rays` would be empty if P has lineality
  V, L = Oscar.minimal_faces(P)
  R = Oscar.rays_modulo_lineality(P).rays_modulo_lineality
  __check_in_d_perp(Iterators.flatten((V, R, L)), d, P)
  return Oscar.convex_hull(
    K,
    __rows(V, keep),
    isempty(R) ? nothing : __rows(R, keep),
    isempty(L) ? nothing : __rows(L, keep),
  )
end

function project_to_d_perp(
  C::Oscar.Cone, d::AbstractVector{Int}; drop::Int=findlast(!=(0), d)
)
  keep = __d_perp_coordinates(d, drop)
  # C = cone(R) + span(L); `rays` would be empty if C has lineality
  R, L = Oscar.rays_modulo_lineality(C)
  __check_in_d_perp(Iterators.flatten((R, L)), d, C)
  return Oscar.positive_hull(
    Oscar.coefficient_field(C),
    isempty(R) ? zeros(Int, 1, length(keep)) : __rows(R, keep), # the zero cone if no rays
    isempty(L) ? nothing : __rows(L, keep),
  )
end

function project_to_d_perp(
  F::Oscar.PolyhedralFan, d::AbstractVector{Int}; drop::Int=findlast(!=(0), d)
)
  return Oscar.polyhedral_fan(
    map(C -> project_to_d_perp(C, d; drop=drop), Oscar.maximal_cones(F))
  )
end

# a stability parameter
function project_to_d_perp(
  theta::AbstractVector{<:Union{Real,Oscar.FieldElem}}, d::AbstractVector{Int};
  drop::Int=findlast(!=(0), d),
)
  keep = __d_perp_coordinates(d, drop)
  __check_in_d_perp((theta,), d, theta)
  return theta[keep]
end

# a list of any of the above, e.g. the output of `vgit_chambers`
function project_to_d_perp(
  X::AbstractVector, d::AbstractVector{Int}; drop::Int=findlast(!=(0), d)
)
  return map(x -> project_to_d_perp(x, d; drop=drop), X)
end

# NB: no `@compile_workload` here. The VGIT functions return polymake C++ objects
# (via Oscar's polyhedral geometry) which cannot be serialized into a precompile
# image: cached objects come back "deleted". Precompiling this path would poison
# the Memoization caches with dead objects, so it is left to compile on first use.

end # module QuiverToolsOscarExt
