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
#   * Boos--Franzen, *GKM-theory for quiver moduli*,
#     [arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049),
#
# which builds on Weist's localisation
#
#   * Weist, *Localization in quiver moduli spaces*,
#     Represent. Theory 17 (2013), 382--425.
#
# A *compatible dimension vector* for d in N^{Q_0} is a function
# beta: Q_0 x Z^m -> N with finite support such that
# sum_{xi} beta(i, xi) = d_i for each i in Q_0. The group Z^m acts on compatible
# dimension vectors by s_chi(beta)(i, xi) = beta(i, xi + chi). The connected
# components of the T-fixed locus are indexed by equivalence classes of compatible
# beta under this action, with the identification
# F_beta = M^{theta_hat}(Q(w), beta), where theta_hat_{i, xi} = theta_i.
########################################################################################

"""
# Summary

`const CoveringDimVector = Dict{Tuple{Int, Vector{Int}}, Int}`

A finitely-supported dimension vector on a covering quiver `Q(w)`.

Keys are pairs `(i, xi)` where `i` is a vertex of `Q` and `xi` is a lattice point
in `Z^{n_arrows(Q)}`; values are positive multiplicities. Only nonzero entries
are stored.

For the role of `CoveringDimVector` in the description of the natural torus fixed
locus of `M^{theta}(Q, d)`, see Section 3 of
[[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)].
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

See Section 3 of [[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)].

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
  shifted = CoveringDimVector()
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
`Q(w) -> Q`: if `beta` and `gamma` have all support at `xi = 0`, the result equals
the Euler form of the underlying dimension vectors on `Q`.

See Section 3 of [[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)].

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
  arrow_list = arrows(Q)
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
in Section 3 of [[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)].

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

For the `(2, 3)` fixed point of `M^{theta}(K_3, (1, 1))` corresponding to the
diagonal arrow, the subquiver is the linear `A_3`:

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
  arrow_list = arrows(Q)
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
    compatible_dimension_vectors(Q::Quiver, d::AbstractVector{Int})

Enumerate the connected-support dimension vectors on `Q(w)` compatible with `d`,
up to the `Z^{n_arrows(Q)}`-shift action.

A `beta::CoveringDimVector` is *compatible* with `d` if
`sum_{xi} beta(i, xi) = d_i` for each vertex `i`. By Section 3 of
[[arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)], the connected
components of the natural-torus fixed locus of `M^{theta}(Q, d)` are indexed by
the shift-equivalence classes of such `beta`.

This function returns one representative per shift-equivalence class
*with connected support*. A dimension vector with disconnected support in
(the underlying graph of) `Q(w)` decomposes any representation as a direct sum,
so its `theta_hat`-stable moduli is empty and it never contributes a fixed point.
Restricting the enumeration to connected supports is mathematically sound for
indexing the fixed locus, but means that compatible `beta` with disconnected
support are intentionally omitted from the output.

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
  arrow_list = arrows(Q)
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

  return [CoveringDimVector(pairs) for pairs in results]
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
    weight_space_dimension(Q::Quiver, beta::CoveringDimVector, chi::AbstractVector{Int})

Dimension of the `chi`-weight space of `Ext^1_{Q(w)}(N, s_{-chi} N)`,
where `N` is a representation of `Q(w)` with dimension vector `beta`.

This is the dimension of the `chi`-weight space of the tangent space at the
fixed point `[M]` in `M^{theta}(Q, d)` corresponding to `beta`, under the
full-rank torus action `T = G_m^{Q_1}`. By
[[Theorem 6.1, arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)],
```math
\\dim (T_{[M]} \\mathcal M)_\\chi
  = \\delta_{\\chi, 0} - \\langle \\beta, s_{-\\chi}\\beta\\rangle_{Q(w)}.
```

Summing the right-hand side over all `chi in Z^{n_arrows(Q)}` recovers
`dim M^{theta}(Q, d) = 1 - \\langle d, d\\rangle_Q` (the proof reduces to
recognising that the double sums of `beta(i, xi) beta(i, xi - chi)` and
`beta(s(a), xi) beta(t(a), xi + e_a - chi)` over `chi` collapse to
`d_i^2` and `d_{s(a)} d_{t(a)}` respectively).

# Input

- `Q::Quiver` a quiver.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.
- `chi::AbstractVector{Int}` a character in `Z^{n_arrows(Q)}`.

# Output

- the dimension as an `Int`.

# Examples

At a real-root fixed point of `M^{theta}(K_2, (1, 1)) \\cong \\mathbb P^1`, the
`chi = 0` weight space is `0`-dimensional (the fixed point is isolated):

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> weight_space_dimension(Q, beta, [0, 0])
0
```

The tangent space at the same fixed point has its single nonzero weight
at `chi = [-1, 1]`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> weight_space_dimension(Q, beta, [-1, 1])
1

julia> weight_space_dimension(Q, beta, [1, -1])
0
```
"""
function weight_space_dimension(
  Q::Quiver, beta::CoveringDimVector, chi::AbstractVector{Int}
)
  delta = all(==(0), chi) ? 1 : 0
  s_neg_chi = shift_beta(beta, .-chi)
  return delta - covering_euler_form(Q, beta, s_neg_chi)
end

"""
    nonzero_weights(Q::Quiver, beta::CoveringDimVector)

Return the characters `chi in Z^{n_arrows(Q)}` with
`weight_space_dimension(Q, beta, chi) > 0`, together with the corresponding dimensions.

By [[Theorem 6.1, arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)],
the support of the character of `T_{[M]} \\mathcal M` for the full-rank torus
is finite, with the only possible nonzero weights of the form
`chi = xi + e_a - xi'` for an arrow `a: s -> t` and points
`(s, xi), (t, xi')` in `supp(beta)`, or `chi = xi - xi'` for points
`(v, xi), (v, xi')` of the same vertex in `supp(beta)` (including `chi = 0`).

# Input

- `Q::Quiver` a quiver.
- `beta::CoveringDimVector` a finitely-supported dimension vector on `Q(w)`.

# Output

- `Vector{Tuple{Vector{Int}, Int}}` of `(chi, dim)` pairs with `dim > 0`.

# Examples

At each of the two fixed points of `\\mathbb P^1 \\cong M^{theta}(K_2, (1, 1))`,
there is exactly one nonzero weight, of dimension `1`:

```jldoctest
julia> Q = kronecker_quiver(2);

julia> beta = CoveringDimVector((1, [0, 0]) => 1, (2, [1, 0]) => 1);

julia> nonzero_weights(Q, beta)
1-element Vector{Tuple{Vector{Int64}, Int64}}:
 ([-1, 1], 1)
```

At a `beta` whose fixed-point component is smooth (e.g. any real-root `beta`
representing an isolated fixed point, as in `M^{theta}(K_m, (1, 1)) \\cong
\\mathbb P^{m-1}`), all weight space dimensions are nonnegative and they sum to
`dim M = 1 - \\langle d, d\\rangle_Q`:

```jldoctest
julia> Q = kronecker_quiver(3); d = [1, 1];

julia> all(compatible_dimension_vectors(Q, d)) do beta
           sum(dim for (_, dim) in nonzero_weights(Q, beta)) == 1 - euler_form(Q, d, d)
       end
true
```

In general, the algebraic identity
```math
\\sum_{\\chi} \\bigl(\\delta_{\\chi, 0} - \\langle \\beta, s_{-\\chi}\\beta\\rangle_{Q(w)}\\bigr)
  = 1 - \\langle d, d\\rangle_Q
```
holds for every compatible `beta`, but individual terms can be negative when
`beta` is not a fixed point of a smooth moduli component.
"""
function nonzero_weights(Q::Quiver, beta::CoveringDimVector)
  out = Tuple{Vector{Int},Int}[]
  for chi in _weight_candidates(Q, beta)
    dim = weight_space_dimension(Q, beta, chi)
    dim > 0 && push!(out, (chi, dim))
  end
  return out
end

"""
    _weight_candidates(Q::Quiver, beta::CoveringDimVector)

Return the (finite) set of characters `chi in Z^{n_arrows(Q)}` outside of which
`weight_space_dimension(Q, beta, chi) = 0`. Used by `nonzero_weights`; exposed for
testing the algebraic sum identity of
[[Theorem 6.1, arXiv:2002.12049](https://doi.org/10.48550/arXiv.2002.12049)].
"""
function _weight_candidates(Q::Quiver, beta::CoveringDimVector)
  arrow_list = arrows(Q)
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
