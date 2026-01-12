"""
    is_special_subdimension_vector(Q::Quiver, e::AbstractVector{Int}, d::AbstractVector{Int})

Compute whether `e` is a special subdimension vector of `d` for `Q`.

# Example

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-4"); d = [1, 1, 1, 1];

julia> is_special_subdimension_vector(Q, [1, 1, 1, 1], d)
false

julia> is_special_subdimension_vector(Q, [1, 0, 0, 1], d)
true
```
"""
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

  cone = polyhedron(
    vcat(eprimes, [e + f for f in fprimes]),
    zeros(Int, length(eprimes) + length(fprimes)),
  )
  cone = intersect(cone, sst(Q, d))

  return !any(
    issubset(cone, polyhedron([eprime, -eprime], [0, 0]))
    for eprime in eprimes
  ) &&
         !any(
    issubset(cone, polyhedron([e + f, - (e + f)], [0, 0]))
    for f in fprimes if f != d - e # this would check if cone ⋐ d^⟂, which is true of course
  )
end
"""
    all_special_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})

Compute all the special subdimension vectors of `d` for the quiver `Q`.

# Example

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> length(all_special_subdimension_vectors(Q, d))
7
```
"""
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
    __helper_accelerate(P::Polyhedron)

Internal method.

Convert the polyhedron `P` to be spanned by its rays.
Only works for strongly convex polyhedra.
"""
function __helper_accelerate(P::Polyhedron)
  newP = polyhedron(positive_hull(rays(P)))
  newP != P && return P
  return newP
end

"""
    sst(Q, e)

Compute the semistable cone `sst(e)` for a given quiver `Q` and a vector `e`.
This is the cone of all stability parameters for which semistable representations
of the quiver `Q` with dimension vector `e` exist.


# Example

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> collect(rays(sst(Q, d)))
3-element Vector{RayVector{QQFieldElem}}:
 [0, 0, 1, -1]
 [0, 1, -1, 0]
 [1, -1, 0, 0]
```
"""
@memoize Dict function sst(Q, e)
  e_perp = polyhedron([e, -e], [0, 0]) #e^{\perp}
  all_gen = filter(
    eprime -> !all(ei == 0 for ei in eprime) && eprime != e,
    all_general_subdimension_vectors(Q, e),
  )
  isempty(all_gen) && return e_perp
  return intersect(e_perp, polyhedron(all_gen, zeros(Int, length(all_gen))))
end

"""
    vgit_walls(Q, d; inner=false, top_dimension=true)

Compute all the walls `W_e` of the quiver `Q` with dimension vector `d`.
Defaults to only computing the top-dimensional walls, pass
the `top_dimension=false` keyword to compute all the `W_e` of the VGIT problem.
Defaults to computing all walls, pass the `inner=true` keyword
to not include the outer walls.

# Example

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> map(rays, vgit_walls(Q, d; inner=false, top_dimension=false))
7-element Vector{SubObjectIterator{RayVector{QQFieldElem}}}:
 [[0, 0, 1, -1], [0, 1, -1, 0]]
 [[0, 0, 1, -1], [1, 0, -1, 0]]
 [[0, 0, 1, -1], [1, -1, 0, 0]]
 [[1, 0, 0, -1], [1, -1, 0, 0]]
 [[1, 0, -1, 0]]
 [[1, 0, 0, -1], [0, 1, -1, 0]]
 [[1, -1, 0, 0], [0, 1, -1, 0]]
```
"""
@memoize Dict function vgit_walls(Q, d; inner=false, top_dimension=true)
  out = map(
    e -> reduce(intersect, [sst(Q, e), sst(Q, d - e), sst(Q, d)]), # definition of W_{e}
    QuiverTools.all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  top_dimension && (out = filter(w -> dim(w) == length(d) - 2, out))
  inner &&
    (out = filter(w -> !any(issubset(w, f) for f in facets(Polyhedron, sst(Q, d))), out))
  return map(__helper_accelerate, unique(out)) # remove duplicates, accelerate intersection
end

"""
    wall_system(Q, d; inner=false, as_cones=true)

Compute the wall system for the quiver `Q` with dimension vector `d`.
This is the set of hyperplanes `H_e` for which there exists at least one
vgit wall `W_e` that lays on `H_e`.
By default it returns the walls intersected with the semistable cone `sst(d)`,
pass the keyword `as_cones=false` to get the full hyperplanes.

"""
function wall_system(Q, d; inner=false, as_cones=true)
  sstd = sst(Q, d)

  all_walls = vgit_walls(Q, d; top_dimension=true)
  if inner
    all_walls = filter(
      w -> !any(issubset(w, f) for f in facets(Polyhedron, sstd)), all_walls
    )
  end
  # each affine hull contains the defining hyperplanes that cut it out.
  all_walls = map(w -> affine_hull(w), all_walls)
  all_walls = map(
    AH -> polyhedron(
      vcat([hyp.a[1, :] for hyp in AH], [-hyp.a[1, :] for hyp in AH]),
      zeros(Int, 2 * length(AH)),
    ),
    all_walls,
  )
  as_cones && return map(w -> intersect(w, sstd), all_walls)
  return all_walls
end

"""
    vgit_chambers(Q, d; verbose=false)

Compute all VGIT chambers for the quiver `Q` with dimension vector `d`.

VGIT chambers are the top-dimensional equivalence classes in the VGIT problem;
they have dimension `length(d) - 1` and are defined by the walls `W_e` of dimension
`length(d) - 2`.

# Example

The following example has three inner walls `W_e`, and all of them are strictly smaller
than the corresponding `H_e \\cap sst(d)` in the wall system.

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-4"); d = [1, 1, 1, 1];

julia> map(rays, vgit_chambers(Q, d; verbose=false))
3-element Vector{SubObjectIterator{RayVector{QQFieldElem}}}:
 [[1, -1, 0, 0], [1, 0, 0, -1], [0, 0, 1, -1]]
 [[0, 1, -1, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
 [[1, -1, 0, 0], [1, 0, 0, -1], [0, 1, -1, 0]]
```

The following example has three inner walls as well, but one is equal to its wall system
hyperplane and the two others are not.

```jldoctests
julia> Q = Quiver("1-2,2-3,3-4,1-3,1-4"); d = [1, 1, 1, 1];

julia> map(rays, vgit_chambers(Q, d; verbose=false))
4-element Vector{SubObjectIterator{RayVector{QQFieldElem}}}:
 [[1, -1, 0, 0], [1, 0, 0, -1], [0, 0, 1, -1]]
 [[1, -1, 0, 0], [1, 0, -1, 0], [1, 0, 0, -1]]
 [[0, 1, -1, 0], [1, 0, 0, -1], [1, 0, -1, 0]]
 [[0, 1, -1, 0], [0, 0, 1, -1], [1, 0, 0, -1]]
```
"""
@memoize Dict function vgit_chambers(Q, d; verbose=false)
  sstd = __helper_accelerate(sst(Q, d))

  # top-dimensional inner walls
  int_walls = vgit_walls(Q, d; top_dimension=true)
  int_walls = filter(w -> !any(issubset(w, f) for f in facets(Polyhedron, sstd)), int_walls)

  # we split the walls into two sets: the ones equal to the wall system hyperplane
  # they lay on, and the ones that are not.
  wallsyst = wall_system(Q, d; inner=true, as_cones=true)

  full_walls = filter(w -> w in wallsyst, int_walls)
  smaller_walls = filter(w -> !(w in wallsyst), int_walls)

  full_walls = map(w -> affine_hull(w), full_walls)
  smaller_walls_ah = map(w -> affine_hull(w), smaller_walls)

  # find a vector not equal to d that cuts out the wall
  full_walls = map(AH -> AH[findfirst(hyp -> (hyp.a[1, :] != d), AH)].a[1, :], full_walls)
  smaller_walls_ah = map(
    AH -> AH[findfirst(hyp -> (hyp.a[1, :] != d), AH)].a[1, :], smaller_walls_ah
  )

  full_walls_iterate = eachindex(full_walls)
  smaller_walls_iterate = eachindex(smaller_walls_ah)

  top_chambers = [sstd]

  # helper function
  helper_split(wall, chamber) = [
    intersect(chamber, polyhedron([wall], [0])),
    intersect(chamber, polyhedron([-wall], [0])),
  ]

  verbose &&
    @info "There are $(length(full_walls)) full walls and $(length(smaller_walls)) smaller walls."

  verbose && @info "Treating the full walls..."
  for i in full_walls_iterate
    top_chambers = vcat(
      map(
        chamber -> helper_split(full_walls[i], chamber),
        top_chambers)...,
    )
    if verbose
      @time top_chambers = filter(
        c -> dim(c) == length(d) - 1, __helper_accelerate.(unique(top_chambers))
      )
    else
      top_chambers = filter(
        c -> dim(c) == length(d) - 1, __helper_accelerate.(unique(top_chambers))
      )
    end
    verbose && @info "Found $(length(top_chambers)) unique chambers after $(i) steps.\n"
  end

  length(smaller_walls) == 0 && return top_chambers

  verbose && @info "Treating the smaller walls..."
  for i in smaller_walls_iterate
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
    if verbose
      @time top_chambers = filter(
        c -> dim(c) == length(d) - 1, __helper_accelerate.(unique(top_chambers))
      )
    else
      top_chambers = filter(
        c -> dim(c) == length(d) - 1, __helper_accelerate.(unique(top_chambers))
      )
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
      dim(inters) != length(d) - 2 && continue
      # if their intersection (the common facet) lays on a W_e, good.
      any(issubset(inters, w) for w in int_walls) && continue
      verbose && @info "Found a fake wall, removing it..."
      incomplete = true
      # otherwise, we replace the two chambers in top_chambers by their union
      push!(top_chambers, minkowski_sum(ch1, ch2))
      top_chambers = filter(ch -> (ch != ch1 && ch != ch2), top_chambers)
    end
  end
  verbose &&
    @info "Done. Found $(length(top_chambers)) chambers after removing the fake walls."
  return top_chambers
end

"""
    vgit_fan(Q, d; verbose=false)

Compute the VGIT fan for the quiver `Q` with dimension vector `d`.
"""
function vgit_fan(Q, d; verbose=false)
  return polyhedral_fan(
    map(ch -> positive_hull(rays(ch)),
      vgit_chambers(Q, d; verbose=verbose),
    ),
  )
end

"""
    git_equivalent(Q, d, theta1, theta2)
Check if the two stability parameters `theta1` and `theta2` are equivalent.

By [Corollary 4.4, arXiv:2506.20568], this is equivalent to their convex hull
either lying in a wall or not intersecting any of them.
"""
function git_equivalent(Q, d, theta1, theta2)
  theta1 == theta2 && return true
  line = convex_hull(theta1, theta2) # 1-dimensional iif theta1 != theta2
  # either the line lies in a wall or it intersects none of them
  return all(
    dim(intersect(w, line)) in [1, -1]
    for w in vgit_walls(Q, d; top_dimension=false)
  )
end
