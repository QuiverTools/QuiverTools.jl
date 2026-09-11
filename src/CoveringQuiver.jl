########################################################################################
# Covering quivers of a finite quiver Q
#
# Given a finite quiver Q with arrow set Q_1 = {a_1, ..., a_m}, the *covering quiver*
# Q(w) is the (infinite) quiver
#
#   vertices:  Q_0 x Z^m,
#   arrows:    (s(a_k), xi) -> (t(a_k), xi + e_k)   for each a_k in Q_1, each xi in Z^m,
#
# where (e_k)_{k=1}^m is the standard basis of Z^m. The projection (i, xi) -> i is the
# universal abelian cover of Q with respect to the free abelian group on Q_1.
#
# These covers describe the fixed loci of the natural action of the full-rank torus
# T = G_m^{Q_1} on M^{theta}(Q, d), see
#
#   * Boos--Franzen, *Weight spaces and attracting sets for torus actions on
#     quiver moduli*, Bull. Lond. Math. Soc. 54 (2022), 1658--1682,
#     [doi:10.1112/blms.12649](https://doi.org/10.1112/blms.12649),
#
# which builds on Weist's localisation
#
#   * Weist, *Localization in quiver moduli spaces*,
#     Represent. Theory 17 (2013), 382--425,
#     [doi:10.1090/S1088-4165-2013-00436-3](https://doi.org/10.1090/S1088-4165-2013-00436-3).
#
# A *compatible dimension vector* for d in N^{Q_0} is a function
# beta: Q_0 x Z^m -> N with finite support such that
# sum_{xi} beta(i, xi) = d_i for each i in Q_0. The group Z^m acts on compatible
# dimension vectors by s_chi(beta)(i, xi) = beta(i, xi + chi). A connected-support
# shift class is a candidate for a component of the T-fixed locus. It contributes
# precisely when the lifted stable moduli space
# F_beta = M^{theta_hat}(Q(w), beta), where theta_hat_{i, xi} = theta_i, is nonempty.
########################################################################################

"""
# Summary

`const CoveringDimVector = Dict{Tuple{Int, Vector{Int}}, Int}`

A finitely-supported dimension vector on a covering quiver `Q(w)`.

Keys are pairs `(i, xi)` where `i` is a vertex of `Q` and `xi` is a lattice point
in `Z^{n_arrows(Q)}`; values are positive multiplicities. Only nonzero entries
are stored.

The lattice coordinates are mutable `Vector`s for compatibility with the rest of
the package. Do not mutate a coordinate while it is used as a dictionary key;
doing so invalidates the dictionary's hash table. Public covering-quiver
operations validate vertex indices, coordinate lengths, and multiplicities.

For the role of `CoveringDimVector` in the description of the natural torus fixed
locus of `M^{theta}(Q, d)`, see
[[Theorem 3.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)].
The shift action of `Z^{n_arrows(Q)}` is realised by [`shift_beta`](@ref);
the enumeration of equivalence classes by [`compatible_dimension_vectors`](@ref).

# Examples

```jldoctest
julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> beta[(1, [0, 0])]
1

julia> length(beta)
2
```
"""
const CoveringDimVector = Dict{Tuple{Int,Vector{Int}},Int}

_covering_arrows(Q::Quiver) = n_arrows(Q) == 0 ? Tuple{Int,Int}[] : arrows(Q)

function _validate_covering_dimension_vector(
  beta::CoveringDimVector,
  lattice_rank::Int;
  vertex_count::Union{Nothing,Int}=nothing,
  name::String="beta",
)
  for ((v, xi), multiplicity) in beta
    v > 0 || throw(ArgumentError("$name contains the nonpositive vertex index $v"))
    if vertex_count !== nothing && v > vertex_count
      throw(
        ArgumentError(
          "$name contains vertex $v, but the quiver has $vertex_count vertices"
        ),
      )
    end
    length(xi) == lattice_rank ||
      throw(
        DimensionMismatch(
          "$name contains a lattice point of length $(length(xi)); " *
          "expected $lattice_rank",
        ),
      )
    multiplicity > 0 ||
      throw(ArgumentError("$name contains the nonpositive multiplicity $multiplicity"))
  end
  return nothing
end

function _validate_covering_dimension_vector(
  Q::Quiver,
  beta::CoveringDimVector;
  name::String="beta",
)
  return _validate_covering_dimension_vector(
    beta,
    n_arrows(Q);
    vertex_count=n_vertices(Q),
    name=name,
  )
end

"""
    shift_beta(beta::CoveringDimVector, chi::AbstractVector{Int})

Compute the shift `s_chi(beta)` of `beta` on the covering quiver `Q(w)`.

For `chi` in `Z^{n_arrows(Q)}`, the shift action is
```math
s_\\chi(\\beta)_{i, \\xi} = \\beta_{i, \\xi + \\chi}.
```
At the level of the underlying `Dict`, this maps each key `(i, xi)` to
`(i, xi - chi)`, so that reading the shifted vector at `(i, eta)` returns the
original multiplicity at `(i, eta + chi)`.

See [[Section 3, Boos--Franzen](https://doi.org/10.1112/blms.12649)].

# Input

- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.
- `chi::AbstractVector{Int}` a lattice shift in `Z^{n_arrows(Q)}`.

# Output

- the shifted dimension vector `s_chi(beta)`.

# Examples

```jldoctest
julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> shift_beta(beta, [1, 0]) == CoveringDimVector((1, [-1, 0]) => 1, (2, [0, 0]) => 1)
true

julia> shift_beta(beta, [0, 0]) == beta
true

julia> shift_beta(shift_beta(beta, [1, 0]), [-1, 0]) == beta
true
```
"""
function shift_beta(beta::CoveringDimVector, chi::AbstractVector{Int})
  _validate_covering_dimension_vector(beta, length(chi))
  shifted = CoveringDimVector()
  sizehint!(shifted, length(beta))
  for ((v, xi), count) in beta
    shifted[(v, xi .- chi)] = count
  end
  return shifted
end

"""
    covering_euler_form(Q::Quiver, beta::CoveringDimVector, gamma::CoveringDimVector)

Compute the Euler form `\\langle \\beta, \\gamma\\rangle_{Q(w)}` on the covering quiver.

For finitely-supported `beta`, `gamma`, this is
```math
\\langle \\beta, \\gamma\\rangle_{Q(w)}
  = \\sum_{i \\in Q_0} \\sum_{\\xi} \\beta_{i, \\xi}\\, \\gamma_{i, \\xi}
  - \\sum_{a\\colon i \\to j} \\sum_{\\xi}
        \\beta_{i, \\xi}\\, \\gamma_{j, \\xi + e_a},
```
where `e_a` is the standard basis vector of `Z^{n_arrows(Q)}` indexing the arrow `a`.

This is the natural lift of [`euler_form`](@ref) along the covering map
`Q(w) -> Q`. If `beta` and `gamma` push forward to dimension vectors `d`
and `e` on `Q`, respectively, then
```math
\\langle d, e\\rangle_Q
  = \\sum_{\\chi \\in \\mathbb Z^{Q_1}}
      \\langle \\beta, s_{-\\chi}\\gamma\\rangle_{Q(w)}.
```
Only finitely many summands on the right are nonzero.

See [[Section 6, Boos--Franzen](https://doi.org/10.1112/blms.12649)].

# Input

- `Q::Quiver` a quiver.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.
- `gamma::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.

# Output

- the Euler form `\\langle \\beta, \\gamma\\rangle_{Q(w)}` as an `Int`.

# Examples

For a real-root fixed point of the moduli space `M^{theta}(K_2, (1, 1))`
the covering Euler form is `1`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> covering_euler_form(Q, beta, beta)
1
```

For finitely-supported `beta`, the covering Euler form equals the underlying
Euler form on the finite subquiver induced on `supp(beta)`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> (sub_Q, sub_d, _) = extract_finite_subquiver(Q, beta);

julia> covering_euler_form(Q, beta, beta) == euler_form(sub_Q, sub_d, sub_d)
true
```
"""
function covering_euler_form(
  Q::Quiver, beta::CoveringDimVector, gamma::CoveringDimVector
)
  _validate_covering_dimension_vector(Q, beta)
  _validate_covering_dimension_vector(Q, gamma; name="gamma")
  arrow_list = _covering_arrows(Q)
  m = n_arrows(Q)

  vertex_sum = 0
  for ((v, xi), b_val) in beta
    g_val = get(gamma, (v, xi), 0)
    vertex_sum += b_val * g_val
  end

  arrow_sum = 0
  for (a_idx, (s, t)) in enumerate(arrow_list)
    e_a = zeros(Int, m)
    e_a[a_idx] = 1
    for ((v, xi), b_val) in beta
      v == s || continue
      g_val = get(gamma, (t, xi .+ e_a), 0)
      arrow_sum += b_val * g_val
    end
  end

  return vertex_sum - arrow_sum
end

"""
    extract_finite_subquiver(Q::Quiver, beta::CoveringDimVector)

Build the finite subquiver of `Q(w)` induced on the support of `beta`.

The returned `sub_Q` is a finite `Quiver` whose vertices are the elements of
`supp(beta)`, ordered lexicographically, and whose arrows are the arrows of `Q(w)`
whose endpoints both lie in `supp(beta)`. The companion `sub_d` records the
multiplicities `beta(i, xi)`.

A stability parameter `theta` on `Q` lifts to `Q(w)` via the projection
`(i, xi) -> i`, i.e. `theta_hat_{i, xi} = theta_i`; this lift is
```julia
sub_theta = [theta[k[1]] for k in sort(collect(keys(beta)))]
```
and is the construction used to identify `F_beta` with `M^{theta_hat}(Q(w), beta)`
in [[Theorem 3.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)].

# Input

- `Q::Quiver` a quiver.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.

# Output

A 3-tuple `(sub_Q, sub_d, vertex_map)`:

- `sub_Q::Quiver` the induced finite subquiver.
- `sub_d::Vector{Int}` the dimension vector on `sub_Q` matching `beta`.
- `vertex_map::Dict{Tuple{Int, Vector{Int}}, Int}` mapping each support point
  `(i, xi)` to its vertex index in `sub_Q`.

# Examples

A real-root fixed point of `M^{theta}(K_2, (1, 1))` has induced subquiver `A_2`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> (sub_Q, sub_d, _) = extract_finite_subquiver(Q, beta);

julia> (n_vertices(sub_Q), n_arrows(sub_Q), sub_d)
(2, 1, [1, 1])
```

For a fixed point of `M^{theta}(K_3, (1, 1))` corresponding to one of the three
arrows, the subquiver is the linear `A_2`:

```jldoctest
julia> Q = kronecker_quiver(3);

julia> beta = CoveringDimVector(
           (1, [0, 0, 0]) => 1,
           (2, [1, 0, 0]) => 1,
       );

julia> (sub_Q, sub_d, _) = extract_finite_subquiver(Q, beta);

julia> (n_vertices(sub_Q), n_arrows(sub_Q), sub_d)
(2, 1, [1, 1])
```
"""
function extract_finite_subquiver(Q::Quiver, beta::CoveringDimVector)
  _validate_covering_dimension_vector(Q, beta)
  arrow_list = _covering_arrows(Q)
  m = n_arrows(Q)

  support_verts = sort(collect(keys(beta)))
  vertex_map = Dict{Tuple{Int,Vector{Int}},Int}()
  for (idx, v) in enumerate(support_verts)
    vertex_map[v] = idx
  end
  n_sub = length(support_verts)

  adj = zeros(Int, n_sub, n_sub)
  for (a_idx, (s, t)) in enumerate(arrow_list)
    e_a = zeros(Int, m)
    e_a[a_idx] = 1
    for (v, xi) in support_verts
      v == s || continue
      target = (t, xi .+ e_a)
      if haskey(vertex_map, target)
        adj[vertex_map[(v, xi)], vertex_map[target]] += 1
      end
    end
  end

  sub_Q = Quiver(adj)
  sub_d = [beta[v] for v in support_verts]
  return (sub_Q, sub_d, vertex_map)
end

"""
    extract_finite_subquiver(
      Q::Quiver,
      beta::CoveringDimVector,
      theta::AbstractVector{Int},
    )

Build the finite subquiver induced on the support of `beta`, together with the
lift of the stability parameter `theta`.

The lifted parameter is defined by `theta_hat[(i, xi)] = theta[i]`. Its entries
are ordered compatibly with the dimension vector returned by
[`extract_finite_subquiver(Q, beta)`](@ref).

# Input

- `Q::Quiver` a quiver.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.
- `theta::AbstractVector{Int}` a stability parameter on `Q`.

# Output

A 4-tuple `(sub_Q, sub_d, sub_theta, vertex_map)`, where the first, second, and
fourth entries are as in [`extract_finite_subquiver(Q, beta)`](@ref), and
`sub_theta` is the lifted stability parameter.

# Examples

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> (_, sub_d, sub_theta, _) = extract_finite_subquiver(Q, beta, [1, -1]);

julia> (sub_d, sub_theta)
([1, 1], [1, -1])
```
"""
function extract_finite_subquiver(
  Q::Quiver,
  beta::CoveringDimVector,
  theta::AbstractVector{Int},
)
  length(theta) == n_vertices(Q) ||
    throw(ArgumentError("length(theta) must equal n_vertices(Q)"))

  sub_Q, sub_d, vertex_map = extract_finite_subquiver(Q, beta)
  sub_theta = Vector{Int}(undef, length(vertex_map))
  for ((v, _), i) in vertex_map
    sub_theta[i] = theta[v]
  end
  return (sub_Q, sub_d, sub_theta, vertex_map)
end

"""
    compatible_dimension_vectors(Q::Quiver, d::AbstractVector{Int})

Enumerate the connected-support dimension vectors on `Q(w)` compatible with `d`,
up to the `Z^{n_arrows(Q)}`-shift action.

A `beta::CoveringDimVector` is *compatible* with `d` if
`sum_{xi} beta(i, xi) = d_i` for each vertex `i`. By
[[Theorem 3.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)] and
[[Theorem 3.8, Weist](https://doi.org/10.1090/S1088-4165-2013-00436-3)], each
nonempty fixed component of the natural torus action on `M^{theta}(Q, d)` arises
from a shift-equivalence class of such `beta`.

This function returns one representative per shift-equivalence class
*with connected support*. A dimension vector with disconnected support in
(the underlying graph of) `Q(w)` decomposes any representation as a direct sum,
so its `theta_hat`-stable moduli is empty and it never contributes a fixed point.
The connected-support classes returned here are candidates: use
[`torus_fixed_components`](@ref) to retain exactly those with nonempty lifted
stable moduli. Compatible `beta` with disconnected support are intentionally
omitted from the output.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector on `Q`.

# Output

- `Vector{CoveringDimVector}` of canonical representatives, one per
  shift-equivalence class of connected-support compatible dimension vectors.

# Examples

For the `r`-Kronecker quiver with `d = (1, 1)` the moduli space is
`P^{r-1}`, so the enumeration recovers the `r` torus-fixed points:

```jldoctest
julia> length(compatible_dimension_vectors(kronecker_quiver(2), [1, 1]))
2

julia> length(compatible_dimension_vectors(kronecker_quiver(3), [1, 1]))
3

julia> length(compatible_dimension_vectors(kronecker_quiver(5), [1, 1]))
5
```

For the 3-Kronecker quiver with `d = (2, 3)` there are 55 shift classes:

```jldoctest
julia> length(compatible_dimension_vectors(kronecker_quiver(3), [2, 3]))
55
```

The zero dimension vector has a single (empty) class:

```jldoctest
julia> compatible_dimension_vectors(kronecker_quiver(2), [0, 0]) == [CoveringDimVector()]
true
```
"""
function compatible_dimension_vectors(Q::Quiver, d::AbstractVector{Int})
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("length(d) must equal n_vertices(Q)"))
  any(<(0), d) && throw(ArgumentError("d must be nonnegative"))

  total_d = sum(d; init=0)
  total_d == 0 && return [CoveringDimVector()]

  m = n_arrows(Q)
  arrow_list = _covering_arrows(Q)
  e_vecs = [_standard_basis_vector(m, k) for k in 1:m]

  # Anchor: place one unit of i_0 at the origin and recurse. After
  # `_normalize`, every shift class has its lex-min i_0 position at the
  # origin, so each class is reached.
  i_0 = findfirst(>(0), d)
  origin = zeros(Int, m)

  initial_beta = CoveringDimVector((i_0, origin) => 1)
  remaining = collect(d)
  remaining[i_0] -= 1

  memo = Set{Tuple{Vector{Tuple{Tuple{Int,Vector{Int}},Int}},Vector{Int}}}()
  results = Set{Vector{Tuple{Tuple{Int,Vector{Int}},Int}}}()

  _build_betas!(results, remaining, initial_beta, memo, i_0, arrow_list, e_vecs)

  canonical_results = sort!(collect(results))
  return [CoveringDimVector(pairs) for pairs in canonical_results]
end

function _standard_basis_vector(m::Int, k::Int)
  v = zeros(Int, m)
  v[k] = 1
  return v
end

"""
    _canonicalize(beta::CoveringDimVector)

Sorted list of `((vertex, lattice_point), multiplicity)` triples used as a
hashable canonical form of `beta`.
"""
_canonicalize(beta::CoveringDimVector) = sort!([(k, v) for (k, v) in beta])

"""
    _memo_key(beta::CoveringDimVector, remaining::Vector{Int})

Hashable key for memoising partial states in `_build_betas!`. The `remaining`
vector is copied defensively because the caller mutates a private copy
between recursive invocations.
"""
function _memo_key(beta::CoveringDimVector, remaining::Vector{Int})
  return (_canonicalize(beta), copy(remaining))
end

"""
    _normalize(beta::CoveringDimVector, i_0::Int)

Shift `beta` so that the lexicographically smallest `xi` with `beta(i_0, xi) > 0`
lands at the origin. Produces the unique canonical representative of `beta`'s
`Z^{n_arrows(Q)}`-shift equivalence class.

The function is a no-op when `beta` has no `i_0` support (e.g. when `beta` is
empty). All callers in this module guarantee `beta(i_0, .) != 0`.
"""
function _normalize(beta::CoveringDimVector, i_0::Int)
  xi_min = nothing
  for ((v, xi), _) in beta
    v == i_0 || continue
    if xi_min === nothing || xi < xi_min
      xi_min = xi
    end
  end
  xi_min === nothing && return beta
  all(==(0), xi_min) && return beta
  return shift_beta(beta, xi_min)
end

"""
    _neighbor_points(beta::CoveringDimVector, arrow_list, e_vecs)

Return the set of `(v, xi)` adjacent in `Q(w)` to `supp(beta)` but not in
`supp(beta)` itself. Used to expand the frontier during enumeration.
"""
function _neighbor_points(
  beta::CoveringDimVector,
  arrow_list::Vector{Tuple{Int,Int}},
  e_vecs::Vector{Vector{Int}},
)
  neighbors = Set{Tuple{Int,Vector{Int}}}()
  for (a_idx, (s, t)) in enumerate(arrow_list)
    e_a = e_vecs[a_idx]
    for (v, xi) in keys(beta)
      if v == s
        candidate = (t, xi .+ e_a)
        haskey(beta, candidate) || push!(neighbors, candidate)
      end
      if v == t
        candidate = (s, xi .- e_a)
        haskey(beta, candidate) || push!(neighbors, candidate)
      end
    end
  end
  return neighbors
end

"""
    _build_betas!(results, remaining, current_beta, memo, i_0, arrow_list, e_vecs)

Recursive worker of `compatible_dimension_vectors`. Places one unit per step,
branching on every `(v, xi)` in `supp(beta) ∪ frontier(beta)` with `remaining[v] > 0`,
and uses memoisation to prune branches that revisit a partial state.

Connected-support completeness: every compatible `beta` with connected support
can be assembled one unit at a time along the connectivity graph of `Q(w)`
starting from `(i_0, 0)`, so every connected-support shift class is reached by
some branch.
"""
function _build_betas!(
  results::Set{Vector{Tuple{Tuple{Int,Vector{Int}},Int}}},
  remaining::Vector{Int},
  current_beta::CoveringDimVector,
  memo::Set{Tuple{Vector{Tuple{Tuple{Int,Vector{Int}},Int}},Vector{Int}}},
  i_0::Int,
  arrow_list::Vector{Tuple{Int,Int}},
  e_vecs::Vector{Vector{Int}},
)
  if all(==(0), remaining)
    push!(results, _canonicalize(_normalize(current_beta, i_0)))
    return nothing
  end

  key = _memo_key(current_beta, remaining)
  key in memo && return nothing
  push!(memo, key)

  candidates = Tuple{Int,Vector{Int}}[]
  for nbr in _neighbor_points(current_beta, arrow_list, e_vecs)
    remaining[nbr[1]] > 0 && push!(candidates, nbr)
  end
  for (v, xi) in keys(current_beta)
    remaining[v] > 0 && push!(candidates, (v, xi))
  end
  isempty(candidates) && return nothing

  for (v, xi) in candidates
    new_beta = copy(current_beta)
    new_beta[(v, xi)] = get(new_beta, (v, xi), 0) + 1
    new_remaining = copy(remaining)
    new_remaining[v] -= 1
    _build_betas!(results, new_remaining, new_beta, memo, i_0, arrow_list, e_vecs)
  end
end

"""
    weight_space_dimension(
      M::QuiverModuliSpace,
      beta::CoveringDimVector,
      chi::AbstractVector{Int},
    )

Return the dimension of the `chi`-weight space of the tangent space along the
stable fixed component indexed by `beta`.

Here `beta` must be compatible with the dimension vector of `M`, and its lifted
stability condition must admit stable representations. For a stable
representation `N` of `Q(w)` with dimension vector `beta`, the weight space is
`Ext^1_{Q(w)}(N, s_{-chi} N)`. By
[[Theorem 6.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)],
```math
\\dim (T_{[M]} \\mathcal M)_\\chi
  = \\delta_{\\chi, 0} - \\langle \\beta, s_{-\\chi}\\beta\\rangle_{Q(w)}.
```

# Input

- `M::QuiverModuliSpace` the ambient stable moduli space, or a semistable moduli
  space whose stable and semistable loci agree.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.
- `chi::AbstractVector{Int}` a character in `Z^{n_arrows(Q)}`.

# Output

- the dimension as an `Int`.

# Examples

At a real-root fixed point of `M^{theta}(K_2, (1, 1)) \\cong \\mathbb P^1`, the
`chi = 0` weight space is `0`-dimensional (the fixed point is isolated):

```jldoctest
julia> Q = kronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1], [1, -1]);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> weight_space_dimension(M, beta, [0, 0])
0
```

The tangent space at the same fixed point has its single nonzero weight
at `chi = [-1, 1]`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1], [1, -1]);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> weight_space_dimension(M, beta, [-1, 1])
1

julia> weight_space_dimension(M, beta, [1, -1])
0
```
"""
function weight_space_dimension(
  M::QuiverModuliSpace,
  beta::CoveringDimVector,
  chi::AbstractVector{Int},
)
  _require_stable_fixed_locus(M)
  _fixed_component_moduli(M, beta) === nothing &&
    throw(ArgumentError("beta does not admit a stable lift for M"))
  return _stable_weight_space_dimension(M.Q, beta, chi)
end

function _stable_weight_space_dimension(
  Q::Quiver,
  beta::CoveringDimVector,
  chi::AbstractVector{Int},
)
  m = n_arrows(Q)
  length(chi) == m || throw(ArgumentError("length(chi) must equal n_arrows(Q)"))

  vertex_sum = 0
  for ((v, xi), beta_value) in beta
    vertex_sum += beta_value * get(beta, (v, xi .- chi), 0)
  end

  arrow_sum = 0
  for (a_idx, (s, t)) in enumerate(_covering_arrows(Q))
    e_a = _standard_basis_vector(m, a_idx)
    for ((v, xi), beta_value) in beta
      v == s || continue
      arrow_sum += beta_value * get(beta, (t, xi .+ e_a .- chi), 0)
    end
  end

  delta = all(iszero, chi) ? 1 : 0
  return delta - vertex_sum + arrow_sum
end

"""
    tangent_weight_multiplicities(
      M::QuiverModuliSpace,
      beta::CoveringDimVector,
    )

Return the characters `chi in Z^{n_arrows(Q)}` with
`weight_space_dimension(M, beta, chi) > 0`, together with their multiplicities.
The zero character is included when the fixed component has positive dimension.

By [[Theorem 6.1, Boos--Franzen](https://doi.org/10.1112/blms.12649)],
the support of the character of `T_{[M]} \\mathcal M` for the full-rank torus
is finite, with the only possible nonzero weights of the form
`chi = xi + e_a - xi'` for an arrow `a: s -> t` and points
`(s, xi), (t, xi')` in `supp(beta)`, or `chi = xi - xi'` for points
`(v, xi), (v, xi')` of the same vertex in `supp(beta)` (including `chi = 0`).

# Input

- `M::QuiverModuliSpace` the ambient stable moduli space, or a semistable moduli
  space whose stable and semistable loci agree.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.

# Output

- `Vector{Tuple{Vector{Int}, Int}}` of `(chi, dim)` pairs with `dim > 0`.

# Examples

At each of the two fixed points of `\\mathbb P^1 \\cong M^{theta}(K_2, (1, 1))`,
there is exactly one nonzero weight, of dimension `1`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1], [1, -1]);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> tangent_weight_multiplicities(M, beta)
1-element Vector{Tuple{Vector{Int64}, Int64}}:
 ([-1, 1], 1)
```

For every stable fixed component, all weight multiplicities are nonnegative and
sum to `dim M = 1 - \\langle d, d\\rangle_Q`:

```jldoctest
julia> Q = kronecker_quiver(3); d = [1, 1];

julia> M = QuiverModuliSpace(Q, d, [1, -1]);

julia> all(torus_fixed_components(M)) do component
           sum(last, tangent_weight_multiplicities(M, component.beta); init=0) ==
             1 - euler_form(Q, d, d)
       end
true
```
"""
function tangent_weight_multiplicities(
  M::QuiverModuliSpace,
  beta::CoveringDimVector,
)
  _require_stable_fixed_locus(M)
  _fixed_component_moduli(M, beta) === nothing &&
    throw(ArgumentError("beta does not admit a stable lift for M"))

  out = Tuple{Vector{Int},Int}[]
  for chi in _weight_candidates(M.Q, beta)
    dim = _stable_weight_space_dimension(M.Q, beta, chi)
    dim > 0 && push!(out, (chi, dim))
  end
  return sort!(out; by=first)
end

"""
    _weight_candidates(Q::Quiver, beta::CoveringDimVector)

Return the (finite) set of characters `chi in Z^{n_arrows(Q)}` outside of which
the tangent-weight expression vanishes. Used internally by
[`tangent_weight_multiplicities`](@ref).
"""
function _weight_candidates(Q::Quiver, beta::CoveringDimVector)
  arrow_list = _covering_arrows(Q)
  m = n_arrows(Q)

  support_by_vertex = Dict{Int,Vector{Vector{Int}}}()
  for ((v, xi), count) in beta
    count > 0 || continue
    push!(get!(support_by_vertex, v, Vector{Int}[]), xi)
  end

  candidates = Set{Vector{Int}}()
  for (a_idx, (s, t)) in enumerate(arrow_list)
    e_a = zeros(Int, m)
    e_a[a_idx] = 1
    s_points = get(support_by_vertex, s, Vector{Int}[])
    t_points = get(support_by_vertex, t, Vector{Int}[])
    for xi in s_points, xi_prime in t_points
      push!(candidates, xi .+ e_a .- xi_prime)
    end
  end
  for (_, points) in support_by_vertex
    for xi in points, xi_prime in points
      push!(candidates, xi .- xi_prime)
    end
  end
  return candidates
end

function _lift_denominator(M::QuiverModuliSpace, vertex_map::Dict)
  function lifted_denominator(sub_d::AbstractVector{Int})
    length(sub_d) == length(vertex_map) ||
      throw(ArgumentError("length(sub_d) must equal length(vertex_map)"))
    projected_d = zeros(Int, n_vertices(M.Q))
    for ((v, _), i) in vertex_map
      projected_d[v] += sub_d[i]
    end
    return M.denom(projected_d)
  end
  return lifted_denominator
end

function _require_stable_fixed_locus(M::QuiverModuliSpace)
  if M.condition == "semistable" && !semistable_equals_stable(M)
    throw(
      ArgumentError(
        "semistable and stable loci must agree; use condition=\"stable\" " *
        "to compute the fixed locus of the stable moduli space",
      ),
    )
  end
  return nothing
end

function _fixed_component_moduli(
  M::QuiverModuliSpace,
  beta::CoveringDimVector,
)
  _validate_covering_dimension_vector(M.Q, beta)

  projected_d = zeros(Int, n_vertices(M.Q))
  for ((v, _), multiplicity) in beta
    projected_d[v] += multiplicity
  end
  projected_d == M.d ||
    throw(ArgumentError("beta is not compatible with the dimension vector of M"))

  sub_Q, sub_d, sub_theta, vertex_map = extract_finite_subquiver(M.Q, beta, M.theta)
  sub_denom = _lift_denominator(M, vertex_map)
  has_stables(sub_Q, sub_d, sub_theta, sub_denom) || return nothing
  return QuiverModuliSpace(sub_Q, sub_d, sub_theta, "stable", sub_denom)
end

"""
    torus_fixed_components(M::QuiverModuliSpace)

Return the connected components of the fixed locus of the natural arrow-scaling
torus on `M`.

The natural torus is `T = G_m^{Q_1}`. Each returned entry is a named tuple
`(beta=beta, moduli=F_beta)`, where `beta` is the canonical representative of a
shift class of compatible covering dimension vectors and `F_beta` is the stable
moduli space on the finite subquiver induced by `supp(beta)`.
Both the stability parameter and a custom denominator are lifted along the
covering projection.

If `M` parametrizes semistable representations, its semistable and stable loci
must agree. For a moduli space constructed with `condition="stable"`, the
function computes the fixed components of the stable locus directly.

# Input

- `M::QuiverModuliSpace` a quiver moduli space.

# Output

- a vector of named tuples `(beta, moduli)` describing the nonempty fixed
  components.

# Examples

The natural torus action on `M^theta(K_2, (1, 1))` has two isolated fixed
points:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(2), [1, 1], [1, -1]);

julia> components = torus_fixed_components(M);

julia> length(components)
2

julia> all(component -> dimension(component.moduli) == 0, components)
true
```

Compatible covering dimension vectors whose lifted stable moduli are empty do
not occur in the output:

```jldoctest
julia> Q = Quiver([0 1; 0 0]);

julia> M = QuiverModuliSpace(Q, [1, 2], [2, -1]);

julia> isempty(torus_fixed_components(M))
true
```
"""
function torus_fixed_components(M::QuiverModuliSpace)
  _require_stable_fixed_locus(M)

  component_type = NamedTuple{(:beta, :moduli),Tuple{CoveringDimVector,QuiverModuliSpace}}
  components = component_type[]
  for beta in compatible_dimension_vectors(M.Q, M.d)
    moduli = _fixed_component_moduli(M, beta)
    moduli === nothing && continue
    push!(components, (beta=beta, moduli=moduli))
  end
  return components
end
