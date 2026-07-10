########################################################################################
# Domokos reductions of quiver-dimension vector pairs
#
# Implements the two reduction operations of
# [[Domokos, Quiver moduli spaces of a given dimension, J. Comb. Algebra (2024)](https://doi.org/10.4171/JCA/97)]:
#
#   - τ_u : merge and delete a *large* vertex u (Definition 2.1, Lemma 3.1);
#   - σ_u : reflect a *small* source or sink u (Definition 2.2, Lemma 3.3).
#
# Both preserve the moduli space up to isomorphism and preserve stability
# (Theorem 2.5): M(Q, α, θ) ≅ M(τ_u Q, τ_u α, τ_u θ), and likewise for σ_u. Iterating
# them toward a `τσ`-minimal representative (Definition 2.3) is Domokos' route to the
# finiteness of moduli spaces of a fixed dimension.
#
# SIGN CONVENTION. Domokos follows King: a representation is θ-semistable when
# θ·dim(S) ≥ 0 for every subrepresentation S, so destabilizing means θ·dim(S) < 0.
# QuiverTools' `slope`/`has_stables` declare subdimension vectors of strictly larger
# slope destabilizing, which for θ·d = 0 means θ·dim(S) > 0. Hence
# θ_QuiverTools = −θ_paper. Every weight
# formula below is *linear* in θ, so negation commutes with it and the formulas
# are written directly in QuiverTools' sign. The one visible consequence is in
# `tau_reduction`: the paper selects its two cases by sign(θ_paper(u)), so under
# negation cases (a)⇄(b) swap. The code keys on sign(θ(u)) in QuiverTools' sign
# and is annotated with the matching paper case.
########################################################################################

# arrows into u weighted by d:  Σ_{tb=u} d(sb) = Σ_v (#v→u)·d(v)   (column u · d)
_in_sum(Q::Quiver, d::AbstractVector{Int}, u::Int) = sum(Q.adjacency[:, u] .* d)
# arrows out of u weighted by d: Σ_{sc=u} d(tc) = Σ_v (#u→v)·d(v)   (row u · d)
_out_sum(Q::Quiver, d::AbstractVector{Int}, u::Int) = sum(Q.adjacency[u, :] .* d)

"""
    is_large(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether vertex `u` is *large* for the pair `(Q, d)` in the sense of
[Definition 2.1, [Domokos](https://doi.org/10.4171/JCA/97)]: `Q` has no loop at `u`,
`u` has positive degree, and

```math
d(u) \\ge \\max\\Bigl\\{\\sum_{b\\colon tb=u} d(sb),\\ \\sum_{c\\colon sc=u} d(tc)\\Bigr\\}.
```

A large vertex can be removed by [`tau_reduction`](@ref).

# Examples

```jldoctest
julia> Q = Quiver("1-2,2---3");  # 1 → 2, and 2 ⇉⇉⇉ 3

julia> is_large(Q, [1, 3, 1], 2)
true

julia> is_large(Q, [1, 2, 1], 2)  # d(2) = 2 < 3 = Σ_{sc=2} d(tc)
false
```
"""
function is_large(Q::Quiver, d::AbstractVector{Int}, u::Int)
  Q.adjacency[u, u] == 0 || return false                  # no loop at u
  indegree(Q, u) + outdegree(Q, u) > 0 || return false    # deg_Q(u) > 0
  return d[u] >= max(_in_sum(Q, d, u), _out_sum(Q, d, u))
end

"""
    is_small_source(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether `u` is a *small source* for `(Q, d)`: a source of `Q` with
`∑_{sa=u} d(ta) > d(u)` [Definition 2.2, [Domokos](https://doi.org/10.4171/JCA/97)].

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v1 is a source

julia> is_small_source(Q, [2, 1, 3], 1)
true
```
"""
is_small_source(Q::Quiver, d::AbstractVector{Int}, u::Int) =
  is_source(Q, u) && _out_sum(Q, d, u) > d[u]

"""
    is_small_sink(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether `u` is a *small sink* for `(Q, d)`: a sink of `Q` with
`∑_{ta=u} d(sa) > d(u)` [Definition 2.2, [Domokos](https://doi.org/10.4171/JCA/97)].

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v2 is a sink

julia> is_small_sink(Q, [2, 1, 3], 2)
true
```
"""
is_small_sink(Q::Quiver, d::AbstractVector{Int}, u::Int) =
  is_sink(Q, u) && _in_sum(Q, d, u) > d[u]

"""
    tau_reduction(Q::Quiver, d, theta, u::Int)
    tau_reduction(M::QuiverModuliSpace, u::Int)

Apply the `τ_u` reduction at a large vertex `u`
[Definition 2.1 and Lemma 3.1, [Domokos](https://doi.org/10.4171/JCA/97)].

The vertex `u` and its adjacent arrows are deleted; for every pair of arrows
`b : v → u` and `c : u → w` a new arrow `v → w` is added (so the number of new arrows
`v → w` is `(#v→u)·(#u→w)`). The dimension vector is restricted to the remaining
vertices. The weight transforms, in QuiverTools' sign convention, by

- `θ(u) > 0` (paper case (b), requires `d(u) = ∑_{sc=u} d(tc)`):
  `(τθ)(v) = θ(v) + (#u→v)·θ(u)`;
- `θ(u) < 0` (paper case (a), requires `d(u) = ∑_{tb=u} d(sb)`):
  `(τθ)(v) = θ(v) + (#v→u)·θ(u)`;
- `θ(u) = 0` (case (c)): `(τθ)(v) = θ(v)`.

Returns the triple `(Q', d', θ')`, or a new `QuiverModuliSpace` when called on `M`.
By [Theorem 2.5, [Domokos](https://doi.org/10.4171/JCA/97)] the moduli space is
unchanged up to isomorphism.

# Examples

Reduce a three-vertex quiver whose moduli space is `ℙ²` to the `3`-Kronecker quiver with
dimension vector `[1, 1]`:

```jldoctest
julia> Q = Quiver("1-2,2---3");

julia> Qr, dr, thetar = tau_reduction(Q, [1, 3, 1], [3, 1, -6], 2);

julia> Qr
Quiver with adjacency matrix [0 3; 0 0]

julia> dr, thetar
([1, 1], [3, -3])
```
"""
function tau_reduction(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, u::Int
)
  is_large(Q, d, u) || throw(ArgumentError("vertex $u is not large for (Q, d)"))
  n = n_vertices(Q)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)

  # new arrows v → w with multiplicity (#v→u)·(#u→w); row/col u are dropped below.
  B = A .+ col * row'
  keep = deleteat!(collect(1:n), u)
  Qnew = Quiver(B[keep, keep])
  dnew = collect(d)[keep]

  tu = theta[u]
  if tu > 0                                       # paper case (b): outgoing side
    d[u] == _out_sum(Q, d, u) ||
      throw(ArgumentError("(Q, d, θ) is not θ-semistable at the large vertex $u"))
    theta_full = collect(theta) .+ row .* tu      # (τθ)(v) = θ(v) + (#u→v)·θ(u)
  elseif tu < 0                                   # paper case (a): incoming side
    d[u] == _in_sum(Q, d, u) ||
      throw(ArgumentError("(Q, d, θ) is not θ-semistable at the large vertex $u"))
    theta_full = collect(theta) .+ col .* tu      # (τθ)(v) = θ(v) + (#v→u)·θ(u)
  else                                            # case (c)
    theta_full = collect(theta)
  end
  return Qnew, dnew, theta_full[keep]
end

"""
    sigma_reduction(Q::Quiver, d, theta, u::Int)
    sigma_reduction(M::QuiverModuliSpace, u::Int)

Apply the `σ_u` reflection at a small source or small sink `u`
[Definition 2.2 and Lemma 3.3, [Domokos](https://doi.org/10.4171/JCA/97)].

All arrows adjacent to `u` are reversed (so a source becomes a sink and vice versa).
The dimension changes only at `u`:

```math
(\\sigma_u d)(u) = -d(u) + \\begin{cases}
  \\sum_{sa=u} d(ta) & u \\text{ a source},\\\\
  \\sum_{ta=u} d(sa) & u \\text{ a sink}.
\\end{cases}
```

The weight transforms by `(σθ)(u) = -θ(u)` and, for `v ≠ u`,
`(σθ)(v) = θ(v) + (#u→v)·θ(u)` if `u` is a source, `θ(v) + (#v→u)·θ(u)` if `u` is a sink.
(These formulas are sign-convention independent.)

`σ_u` is an involution. Returns `(Q', d', θ')`, or a new `QuiverModuliSpace` on `M`.

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v1 a small source for [2, 1, 3]

julia> Qr, dr, thetar = sigma_reduction(Q, [2, 1, 3], [2, -1, -1], 1);

julia> dr, thetar
([5, 1, 3], [-2, 1, 3])
```
"""
function sigma_reduction(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, u::Int
)
  n = n_vertices(Q)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)
  tu = theta[u]

  if is_small_source(Q, d, u)
    dnew_u = -d[u] + _out_sum(Q, d, u)
    theta_new = collect(theta) .+ row .* tu    # source: θ(v) + (#u→v)·θ(u)
  elseif is_small_sink(Q, d, u)
    dnew_u = -d[u] + _in_sum(Q, d, u)
    theta_new = collect(theta) .+ col .* tu    # sink:   θ(v) + (#v→u)·θ(u)
  else
    throw(ArgumentError("vertex $u is not a small source or small sink for (Q, d)"))
  end

  B = copy(A)
  B[:, u] = row    # incoming arrows ← old outgoing
  B[u, :] = col    # outgoing arrows ← old incoming
  dnew = collect(d)
  dnew[u] = dnew_u
  theta_new[u] = -tu                            # (σθ)(u) = -θ(u)
  return Quiver(B), dnew, theta_new
end

for reduction in (:tau_reduction, :sigma_reduction)
  @eval function $reduction(M::QuiverModuliSpace, u::Int)
    Q, d, theta = $reduction(M.Q, M.d, M.theta, u)
    return QuiverModuliSpace(Q, d, theta, M.condition, M.denom)
  end
end

"""
    tau_sigma_reduce(Q::Quiver, d, theta)

Greedily reduce `(Q, d, θ)` to a smaller, moduli-isomorphic representative by repeatedly

1. applying [`tau_reduction`](@ref) at any large vertex (this drops a vertex), then
2. applying [`sigma_reduction`](@ref) at any small source/sink whose reflection *shrinks*
   the dimension vector (i.e. `∑ neighbours < 2·d(u)`).

Both steps strictly decrease `(|Q₀|, |d|)` lexicographically, so this terminates. The
result has no large vertex and no dimension-shrinking small source/sink. Note this greedy
descent is weaker than `τσ`-minimality (Definition 2.3), which also permits `σ`-steps that
temporarily *raise* `|d|`; see [`is_taus_minimal`](@ref).

Returns the reduced triple `(Q', d', θ')`.

# Examples

The `2`-cycle with dimension `[1, 1]` reduces to a single vertex with one loop:

```jldoctest
julia> Q = Quiver("1-2,2-1");

julia> Qr, dr, thetar = tau_sigma_reduce(Q, [1, 1], [1, -1]);

julia> Qr, dr
(Quiver with adjacency matrix [1;;], [1])
```
"""
function tau_sigma_reduce(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}
)
  d = collect(d)
  theta = collect(theta)
  while true
    n = n_vertices(Q)
    u = findfirst(v -> is_large(Q, d, v), 1:n)
    if !isnothing(u)
      Q, d, theta = tau_reduction(Q, d, theta, u)
      continue
    end
    u = findfirst(1:n) do v
      (is_small_source(Q, d, v) && _out_sum(Q, d, v) < 2 * d[v]) ||
        (is_small_sink(Q, d, v) && _in_sum(Q, d, v) < 2 * d[v])
    end
    if !isnothing(u)
      Q, d, theta = sigma_reduction(Q, d, theta, u)
      continue
    end
    return Q, d, theta
  end
end

"""
    is_taus_minimal(Q::Quiver, d, theta; max_states::Int = 10_000)

Decide whether `(Q, d)` is `τσ`-minimal among all sincere quiver-dimension vector pairs
[Definition 2.3, [Domokos](https://doi.org/10.4171/JCA/97)]: whether *no* sequence of
`τ`/`σ` reductions reaches a pair `(Q', d')` with `|Q'₀| < |Q₀|`, or `|Q'₀| = |Q₀|` and
`|d'| < |d|`.

This explores the reduction graph breadth-first. Since `τ` drops a vertex and a
dimension-lowering `σ` is immediately witnessed, a *negative* answer (`false`) is always a
genuine witness. A `σ`-orbit can be infinite (Section 9 of the reference is `τσ`-minimal
yet has `σ`-steps that raise `|d|` without bound), so the search is capped at `max_states`
pairs; if the cap is hit the function returns `true` with a warning, meaning "not disproved
within the search bound".

# Examples

A pair with a large vertex is never minimal (its `τ` reduction drops a vertex):

```jldoctest
julia> is_taus_minimal(Quiver("1-2,2-1"), [1, 1], [1, -1])
false
```

A single vertex carrying a loop is trivially minimal:

```jldoctest
julia> is_taus_minimal(Quiver([1;;]), [1], [0])
true
```
"""
function is_taus_minimal(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}; max_states::Int=10_000
)
  start = (n_vertices(Q), sum(d))
  seen = Set{Tuple{Matrix{Int},Vector{Int}}}()
  queue = [(Q, collect(d), collect(theta))]
  push!(seen, (Matrix{Int}(Q.adjacency), collect(d)))

  while !isempty(queue)
    Qc, dc, tc = popfirst!(queue)
    n = n_vertices(Qc)
    neighbours = Tuple{Quiver,Vector{Int},Vector{Int}}[]
    for v in 1:n
      is_large(Qc, dc, v) && push!(neighbours, tau_reduction(Qc, dc, tc, v))
      (is_small_source(Qc, dc, v) || is_small_sink(Qc, dc, v)) &&
        push!(neighbours, sigma_reduction(Qc, dc, tc, v))
    end
    for (Qn, dn, tn) in neighbours
      (n_vertices(Qn), sum(dn)) < start && return false   # genuine witness
      key = (Matrix{Int}(Qn.adjacency), dn)
      key in seen && continue
      if length(seen) >= max_states
        @warn "is_taus_minimal: search truncated at $max_states pairs; " *
          "returning `true` unproven"
        return true
      end
      push!(seen, key)
      push!(queue, (Qn, dn, tn))
    end
  end
  return true
end
