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
# The paper assumes throughout that the dimension vector is sincere, i.e. nonzero at
# every vertex. The public predicates and reductions below assert this hypothesis.
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

function __check_sincere_dimension_vector(Q::Quiver, d::AbstractVector{Int})
  __check_dimension_vector(Q, d)
  @assert all(>(0), d) "dimension vector must be sincere"
  return nothing
end

function __checked_dot(a::AbstractVector{Int}, b::AbstractVector{Int})
  length(a) == length(b) || throw(DimensionMismatch("vectors must have equal lengths"))
  return foldl(eachindex(a); init=0) do total, i
    Base.Checked.checked_add(total, Base.Checked.checked_mul(a[i], b[i]))
  end
end

__checked_sum(a::AbstractVector{Int}) =
  foldl(Base.Checked.checked_add, a; init=0)

function __check_domokos_weight(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}
)
  __check_sincere_dimension_vector(Q, d)
  length(theta) == n_vertices(Q) ||
    throw(ArgumentError("stability parameter must have length $(n_vertices(Q))"))
  __checked_dot(theta, d) == 0 ||
    throw(ArgumentError("the stability parameter theta must satisfy theta . d = 0"))
  has_semistables(Q, d, theta) ||
    throw(ArgumentError("the dimension vector must be theta-semistable"))
  return nothing
end

# arrows into u weighted by d:  Σ_{tb=u} d(sb) = Σ_v (#v→u)·d(v)   (column u · d)
__in_sum(Q::Quiver, d::AbstractVector{Int}, u::Int) =
  __checked_dot(Q.adjacency[:, u], d)
# arrows out of u weighted by d: Σ_{sc=u} d(tc) = Σ_v (#u→v)·d(v)   (row u · d)
__out_sum(Q::Quiver, d::AbstractVector{Int}, u::Int) =
  __checked_dot(Q.adjacency[u, :], d)

"""
    is_large(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether vertex `u` is *large* for the pair `(Q, d)` in the sense of
[Definition 2.1, [Domokos](https://doi.org/10.4171/JCA/97)]: `Q` has no loop at `u`,
`u` has positive degree, and

```math
d(u) \\ge \\max\\Bigl\\{\\sum_{b\\colon tb=u} d(sb),\\ \\sum_{c\\colon sc=u} d(tc)\\Bigr\\}.
```

A large vertex can be removed by [`tau_reduction`](@ref).
As in the cited definition, `d` must be sincere; an `AssertionError` is thrown
otherwise.

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
  __check_sincere_dimension_vector(Q, d)
  Q.adjacency[u, u] == 0 || return false                  # no loop at u
  indegree(Q, u) + outdegree(Q, u) > 0 || return false    # deg_Q(u) > 0
  return d[u] >= max(__in_sum(Q, d, u), __out_sum(Q, d, u))
end

"""
    is_small_source(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether `u` is a *small source* for `(Q, d)`: a source of `Q` with
``\\sum_{a\\colon sa=u} d(ta) > d(u)``
[Definition 2.2, [Domokos](https://doi.org/10.4171/JCA/97)].
As in the cited definition, `d` must be sincere; an `AssertionError` is thrown
otherwise.

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v1 is a source

julia> is_small_source(Q, [2, 1, 3], 1)
true
```
"""
function is_small_source(Q::Quiver, d::AbstractVector{Int}, u::Int)
  __check_sincere_dimension_vector(Q, d)
  return is_source(Q, u) && __out_sum(Q, d, u) > d[u]
end

"""
    is_small_sink(Q::Quiver, d::AbstractVector{Int}, u::Int)

Check whether `u` is a *small sink* for `(Q, d)`: a sink of `Q` with
``\\sum_{a\\colon ta=u} d(sa) > d(u)``
[Definition 2.2, [Domokos](https://doi.org/10.4171/JCA/97)].
As in the cited definition, `d` must be sincere; an `AssertionError` is thrown
otherwise.

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v2 is a sink

julia> is_small_sink(Q, [2, 1, 3], 2)
true
```
"""
function is_small_sink(Q::Quiver, d::AbstractVector{Int}, u::Int)
  __check_sincere_dimension_vector(Q, d)
  return is_sink(Q, u) && __in_sum(Q, d, u) > d[u]
end

"""
    tau_reduction(Q::Quiver, d, u::Int)
    tau_reduction(Q::Quiver, d, theta, u::Int)
    tau_reduction(M::QuiverModuliSpace, u::Int)

Apply the ``\\tau_u`` reduction at a large vertex `u`
[Definition 2.1 and Lemma 3.1, [Domokos](https://doi.org/10.4171/JCA/97)].

The pair-only method implements Domokos's operation on `(Q, d)`.
The vertex `u` and its adjacent arrows are deleted; for every pair of arrows
``b\\colon v \\to u`` and ``c\\colon u \\to w`` a new arrow ``v \\to w`` is added, so
that ``r_{vu} r_{uw}`` new arrows ``v \\to w`` appear, where ``r_{vu}`` denotes the
number of arrows ``v \\to u``. The dimension vector is restricted to the remaining
vertices. The weight transforms, in QuiverTools' sign convention, by

- ``\\theta(u) > 0`` (paper case (b), requires ``d(u) = \\sum_{c\\colon sc=u} d(tc)``):
  ``(\\tau\\theta)(v) = \\theta(v) + r_{uv}\\theta(u)``;
- ``\\theta(u) < 0`` (paper case (a), requires ``d(u) = \\sum_{b\\colon tb=u} d(sb)``):
  ``(\\tau\\theta)(v) = \\theta(v) + r_{vu}\\theta(u)``;
- ``\\theta(u) = 0`` (case (c)): ``(\\tau\\theta)(v) = \\theta(v)``.

Returns the pair `(Q', d')`, the triple `(Q', d', θ')`, or a new
`QuiverModuliSpace`, according to the method called.
By [Theorem 2.5, [Domokos](https://doi.org/10.4171/JCA/97)] the moduli space is
unchanged up to isomorphism. For the weight-aware methods, `d` must be
`theta`-semistable and `theta` must satisfy ``\\theta \\cdot d = 0``, i.e. King's
normalization. The transformed weight is normalized as well. A reduced moduli object
uses the standard denominator `sum`; a denominator closure on the original vertex set
is not reused after that set changes.

All methods require a sincere dimension vector, as in the cited results; an
`AssertionError` is thrown otherwise.

# Examples

Reduce a three-vertex quiver whose moduli space is ``\\mathbb{P}^2`` to the
`3`-Kronecker quiver with dimension vector `[1, 1]`:

```jldoctest
julia> Q = Quiver("1-2,2---3");

julia> Qr, dr, thetar = tau_reduction(Q, [1, 3, 1], [3, 1, -6], 2);

julia> Qr
Quiver with adjacency matrix [0 3; 0 0]

julia> dr, thetar
([1, 1], [3, -3])
```
"""
function tau_reduction(Q::Quiver, d::AbstractVector{Int}, u::Int)
  __check_sincere_dimension_vector(Q, d)
  is_large(Q, d, u) || throw(ArgumentError("vertex $u is not large for (Q, d)"))
  n = n_vertices(Q)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)

  # new arrows v → w with multiplicity (#v→u)·(#u→w); row/col u are dropped below.
  B = similar(A)
  for v in axes(A, 1), w in axes(A, 2)
    B[v, w] = Base.Checked.checked_add(
      A[v, w], Base.Checked.checked_mul(col[v], row[w])
    )
  end
  keep = deleteat!(collect(1:n), u)
  Qnew = Quiver(B[keep, keep])
  dnew = collect(d)[keep]
  return Qnew, dnew
end

function tau_reduction(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, u::Int
)
  __check_domokos_weight(Q, d, theta)
  Qnew, dnew = tau_reduction(Q, d, u)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)
  keep = deleteat!(collect(1:n_vertices(Q)), u)

  tu = theta[u]
  if tu > 0                                       # paper case (b): outgoing side
    d[u] == __out_sum(Q, d, u) ||
      throw(ArgumentError("(Q, d, θ) is not θ-semistable at the large vertex $u"))
    theta_full = [
      Base.Checked.checked_add(theta[v], Base.Checked.checked_mul(row[v], tu)) for
      v in eachindex(theta)
    ]                                             # (τθ)(v) = θ(v) + (#u→v)·θ(u)
  elseif tu < 0                                   # paper case (a): incoming side
    d[u] == __in_sum(Q, d, u) ||
      throw(ArgumentError("(Q, d, θ) is not θ-semistable at the large vertex $u"))
    theta_full = [
      Base.Checked.checked_add(theta[v], Base.Checked.checked_mul(col[v], tu)) for
      v in eachindex(theta)
    ]                                             # (τθ)(v) = θ(v) + (#v→u)·θ(u)
  else                                            # case (c)
    theta_full = collect(theta)
  end
  theta_new = theta_full[keep]
  @assert __checked_dot(theta_new, dnew) == 0 "tau reduction must preserve King normalization"
  return Qnew, dnew, theta_new
end

"""
    sigma_reduction(Q::Quiver, d, u::Int)
    sigma_reduction(Q::Quiver, d, theta, u::Int)
    sigma_reduction(M::QuiverModuliSpace, u::Int)

Apply the ``\\sigma_u`` reflection at a small source or small sink `u`
[Definition 2.2 and Lemma 3.3, [Domokos](https://doi.org/10.4171/JCA/97)].

All arrows adjacent to `u` are reversed (so a source becomes a sink and vice versa).
The dimension changes only at `u`:

```math
(\\sigma_u d)(u) = -d(u) + \\begin{cases}
  \\sum_{sa=u} d(ta) & u \\text{ a source},\\\\
  \\sum_{ta=u} d(sa) & u \\text{ a sink}.
\\end{cases}
```

The pair-only method implements Domokos's operation on `(Q, d)`.
The weight transforms by ``(\\sigma\\theta)(u) = -\\theta(u)`` and, for ``v \\neq u``,
``(\\sigma\\theta)(v) = \\theta(v) + r_{uv}\\theta(u)`` if `u` is a source and
``\\theta(v) + r_{vu}\\theta(u)`` if `u` is a sink, where ``r_{uv}`` denotes the number
of arrows ``u \\to v``. (These formulas are sign-convention independent.)

``\\sigma_u`` is an involution. Returns `(Q', d')`, `(Q', d', θ')`, or a new
`QuiverModuliSpace`, according to the method called.
For the weight-aware methods, `d` must be `theta`-semistable and the stability
parameter must satisfy ``\\theta \\cdot d = 0``, i.e. King's normalization. The
transformed weight is normalized as well, and a reduced moduli object uses the
standard denominator `sum`.

All methods require a sincere dimension vector, as in the cited results; an
`AssertionError` is thrown otherwise.

# Examples

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");  # v1 a small source for [2, 1, 3]

julia> Qr, dr, thetar = sigma_reduction(Q, [2, 1, 3], [2, -1, -1], 1);

julia> dr, thetar
([5, 1, 3], [-2, 1, 3])
```
"""
function sigma_reduction(Q::Quiver, d::AbstractVector{Int}, u::Int)
  __check_sincere_dimension_vector(Q, d)
  n = n_vertices(Q)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)

  if is_small_source(Q, d, u)
    dnew_u = Base.Checked.checked_sub(__out_sum(Q, d, u), d[u])
  elseif is_small_sink(Q, d, u)
    dnew_u = Base.Checked.checked_sub(__in_sum(Q, d, u), d[u])
  else
    throw(ArgumentError("vertex $u is not a small source or small sink for (Q, d)"))
  end

  B = copy(A)
  B[:, u] = row    # incoming arrows ← old outgoing
  B[u, :] = col    # outgoing arrows ← old incoming
  dnew = collect(d)
  dnew[u] = dnew_u
  return Quiver(B), dnew
end

function sigma_reduction(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, u::Int
)
  __check_domokos_weight(Q, d, theta)
  Qnew, dnew = sigma_reduction(Q, d, u)
  A = Matrix{Int}(Q.adjacency)
  col = A[:, u]    # col[v] = #(v → u)
  row = A[u, :]    # row[v] = #(u → v)
  tu = theta[u]
  theta_new = if is_small_source(Q, d, u)
    [
      Base.Checked.checked_add(theta[v], Base.Checked.checked_mul(row[v], tu)) for
      v in eachindex(theta)
    ]                                # source: θ(v) + (#u→v)·θ(u)
  else
    [
      Base.Checked.checked_add(theta[v], Base.Checked.checked_mul(col[v], tu)) for
      v in eachindex(theta)
    ]                                # sink:   θ(v) + (#v→u)·θ(u)
  end
  theta_new[u] = Base.Checked.checked_neg(tu)    # (σθ)(u) = -θ(u)
  @assert __checked_dot(theta_new, dnew) == 0 "sigma reduction must preserve King normalization"
  return Qnew, dnew, theta_new
end

for reduction in (:tau_reduction, :sigma_reduction)
  @eval function $reduction(M::QuiverModuliSpace, u::Int)
    Q, d, theta = $reduction(M.Q, M.d, M.theta, u)
    return QuiverModuliSpace(Q, d, theta, M.condition)
  end
end

"""
    tau_sigma_reduce(Q::Quiver, d)
    tau_sigma_reduce(Q::Quiver, d, theta)

Greedily reduce `(Q, d, θ)` to a smaller, moduli-isomorphic representative by repeatedly

1. applying [`tau_reduction`](@ref) at any large vertex (this drops a vertex), then
2. applying [`sigma_reduction`](@ref) at any small source or sink whose reflection
   *shrinks* the dimension vector, i.e. whose weighted neighbour sum is less than
   ``2d(u)``.

Both steps strictly decrease ``(\\#Q_0, |d|)`` lexicographically, so this terminates.
The result has no large vertex and no dimension-shrinking small source or sink. Note
this greedy descent is weaker than ``\\tau\\sigma``-minimality (Definition 2.3), which
also permits ``\\sigma``-steps that temporarily *raise* ``|d|``; see
[`is_taus_minimal`](@ref).

Returns the reduced pair `(Q', d')` or triple `(Q', d', θ')`, according to the
method called. The dimension vector must be sincere. The weight-aware method also
requires `d` to be `theta`-semistable and `theta` to be King-normalized; every
intermediate weight remains normalized.

# Examples

The `2`-cycle with dimension `[1, 1]` reduces to a single vertex with one loop:

```jldoctest
julia> Q = Quiver("1-2,2-1");

julia> Qr, dr, thetar = tau_sigma_reduce(Q, [1, 1], [1, -1]);

julia> Qr, dr
(Quiver with adjacency matrix [1;;], [1])
```
"""
function tau_sigma_reduce(Q::Quiver, d::AbstractVector{Int})
  __check_sincere_dimension_vector(Q, d)
  d = collect(d)
  while true
    n = n_vertices(Q)
    u = findfirst(v -> is_large(Q, d, v), 1:n)
    if !isnothing(u)
      Q, d = tau_reduction(Q, d, u)
      continue
    end
    u = findfirst(1:n) do v
      (is_small_source(Q, d, v) && __out_sum(Q, d, v) < 2 * d[v]) ||
        (is_small_sink(Q, d, v) && __in_sum(Q, d, v) < 2 * d[v])
    end
    if !isnothing(u)
      Q, d = sigma_reduction(Q, d, u)
      continue
    end
    return Q, d
  end
end

function tau_sigma_reduce(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}
)
  __check_domokos_weight(Q, d, theta)
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
      (is_small_source(Q, d, v) && __out_sum(Q, d, v) < 2 * d[v]) ||
        (is_small_sink(Q, d, v) && __in_sum(Q, d, v) < 2 * d[v])
    end
    if !isnothing(u)
      Q, d, theta = sigma_reduction(Q, d, theta, u)
      continue
    end
    return Q, d, theta
  end
end

"""
    is_taus_minimal(Q::Quiver, d; max_states::Int = 10_000)

Decide whether `(Q, d)` is ``\\tau\\sigma``-minimal in the class of all sincere
quiver-dimension vector pairs [Definition 2.3,
[Domokos](https://doi.org/10.4171/JCA/97)]. That is, decide whether no finite sequence
of ``\\tau`` and ``\\sigma`` reductions reaches a pair `(Q', d')` with fewer vertices,
or with the same number of vertices and smaller total dimension.

This is a property of `(Q, d)` alone; a stability parameter is not part of the cited
definition. The function explores the reduction graph breadth-first. A `false` result
comes with a reduction witness, and `true` is returned only after the reachable graph
has been exhausted. A ``\\sigma``-orbit can be infinite, so the search may instead
reach `max_states`. In that case an `ArgumentError` is thrown because minimality is
still undetermined. The paper does not supply a general search bound: for its Section 9
example, it proves minimality separately by inequalities that describe every possible
sequence of reflections.

The dimension vector must be sincere, as in the cited definition; an `AssertionError`
is thrown otherwise. `max_states` must be positive.

# Examples

A pair with a large vertex is never minimal (its ``\\tau`` reduction drops a vertex):

```jldoctest
julia> is_taus_minimal(Quiver("1-2,2-1"), [1, 1])
false
```

A single vertex carrying a loop is trivially minimal:

```jldoctest
julia> is_taus_minimal(Quiver([1;;]), [1])
true
```

A state limit is not a proof of minimality. In this four-vertex example, the first
reflection does not make the pair smaller, but it exposes a later ``\\tau`` reduction:

```jldoctest
julia> Q = Quiver([0 0 1 0; 0 0 2 2; 0 1 0 1; 0 2 1 0]);

julia> d = [3, 2, 6, 1];

julia> try
           is_taus_minimal(Q, d; max_states=1)
       catch error
           error isa ArgumentError
       end
true

julia> is_taus_minimal(Q, d; max_states=100)
false
```

The pair in Section 9 of the paper is minimal, but its infinite reflection orbit also
shows why bounded breadth-first search need not prove that fact:

```jldoctest
julia> Q = Quiver("1--3,1-2,3-2");

julia> try
           is_taus_minimal(Q, [2, 1, 3]; max_states=2)
       catch error
           error isa ArgumentError
       end
true
```
"""
function is_taus_minimal(
  Q::Quiver, d::AbstractVector{Int}; max_states::Int=10_000
)
  __check_sincere_dimension_vector(Q, d)
  max_states > 0 || throw(ArgumentError("max_states must be positive"))
  start = (n_vertices(Q), __checked_sum(d))
  seen = Set{Tuple{Matrix{Int},Vector{Int}}}()
  queue = [(Q, collect(d))]
  head = 1
  push!(seen, (Matrix{Int}(Q.adjacency), collect(d)))

  while head <= length(queue)
    Qc, dc = queue[head]
    head += 1
    n = n_vertices(Qc)
    neighbours = Tuple{Quiver,Vector{Int}}[]
    for v in 1:n
      is_large(Qc, dc, v) && push!(neighbours, tau_reduction(Qc, dc, v))
      (is_small_source(Qc, dc, v) || is_small_sink(Qc, dc, v)) &&
        push!(neighbours, sigma_reduction(Qc, dc, v))
    end
    for (Qn, dn) in neighbours
      (n_vertices(Qn), __checked_sum(dn)) < start && return false   # genuine witness
      key = (Matrix{Int}(Qn.adjacency), dn)
      key in seen && continue
      if length(seen) >= max_states
        throw(
          ArgumentError(
            "max_states=$max_states reached before the reduction graph was exhausted; " *
            "minimality is undetermined",
          ),
        )
      end
      push!(seen, key)
      push!(queue, (Qn, dn))
    end
  end
  return true
end
