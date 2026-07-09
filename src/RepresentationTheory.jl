##########################################################################
# Methods dealing with various representation-theoretic aspects of quivers
##########################################################################

"""
    euler_matrix(Q::Quiver)

Compute the Euler matrix of `Q`.

The Euler matrix of a quiver ``Q`` is defined as
```math
E = I - A,
```
where ``A`` is the adjacency matrix of ``Q`` and ``I``
is the identity matrix of the same size as ``A``.

# Input

- `Q::Quiver` a quiver.

# Output

- the Euler matrix of the quiver.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> euler_matrix(Q) == [1 -4; 0 1]
true
```
"""
@memoize Dict euler_matrix(Q::Quiver) = identity_matrix(n_vertices(Q)) - Q.adjacency

"""
    euler_form(Q::Quiver, x::AbstractVector{Int}, y::AbstractVector{Int})

Compute the Euler form of `Q` for `x` and `y`.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where ``E`` is the Euler matrix of the quiver.

# Input

- `Q::Quiver` a quiver.
- `x::AbstractVector{Int}`: a vector.
- `y::AbstractVector{Int}`: a vector.

# Output

- the Euler form ``\\langle x, y\\rangle_{Q}``.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> euler_form(Q, [1, 1], [1, 1]) == -2
true
```
"""
euler_form(Q::Quiver, x::AbstractVector{Int}, y::AbstractVector{Int}) =
  x' * y - x' * (Q.adjacency * y)
# this inlines x' * (I - adjacency) * y: retrieving the memoized Euler matrix
# costs more than recomputing the two products

########################################################################################
# Canonical decomposition
########################################################################################

"""
    general_ext(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})

Compute the dimension of the ``\\mathrm{Ext}^1`` group between general representations
of dimension vectors `a` and `b`.

According to [[Theorem 5.4, MR1162487]
(https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487)],
we have

```math
ext(a,b)=max\\{-\\langle c,b\\rangle~~|~~c~\\text{is a general subdimension vector of }a\\}.
```

# Input

- `Q::Quiver` a quiver.
- `a::AbstractVector{Int}`: a vector.
- `b::AbstractVector{Int}`: a vector.

# Output

- the dimension of the general extensions ``\\mathrm{ext}^1(a, b)``.

# Examples

```jldoctest
julia> Q1 = kronecker_quiver(3);

julia> general_ext(Q1, [2, 3], [6, 7])
9

julia> general_ext(Q1, [1, 1], [1, 0])
0

julia> Q2 = three_vertex_quiver(1, 6, 7);

julia> general_ext(Q2, [5, 6, 7], [6, 7, 8])
483
```
"""
function general_ext(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
  if sum(a) <= sum(b)
    return maximum(-euler_form(Q, c, b) for c in all_general_subdimension_vectors(Q, a))
  else
    return maximum(-euler_form(Q, a, b - c) for c in all_general_subdimension_vectors(Q, b))
  end
end

"""
    general_hom(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})

Compute the dimension of the ``\\mathrm{Hom}`` group between general representations
of dimension vectors `a` and `b`.

# Input

- `Q::Quiver` a quiver.
- `a::AbstractVector{Int}`: a vector.
- `b::AbstractVector{Int}`: a vector.

# Output

- the dimension of the general homomorphisms ``\\mathrm{hom}(a, b)``.

# Examples

```jldoctest
julia> Q1 = kronecker_quiver(3);

julia> general_hom(Q1, [2, 3], [6, 7])
0

julia> general_hom(Q1, [1, 1], [1, 0])
1

julia> Q2 = three_vertex_quiver(1, 6, 7);

julia> general_hom(Q2, [5, 6, 7], [6, 7, 8])
0
```
"""
function general_hom(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
  return euler_form(Q, a, b) + general_ext(Q, a, b)
end

"""
    canonical_decomposition(Q::Quiver, d)

Compute the canonical decomposition of `d` for the quiver `Q`.

If ``\\beta_1, \\dots, \\beta_{\\ell}`` is a sequence
of Schur roots such that, for all ``i \\neq j``, one has

```math
\\mathrm{ext}(\\beta_i, \\beta_j) = \\mathrm{ext}(\\beta_j, \\beta_i) = 0,
```

then the general representation of dimension ``\\sum_i \\beta_i`` is
isomorphic to the direct sum of irreducible representations
of dimension vectors ``\\beta_i``.

Such a decomposition is called the canonical decomposition.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.

# Output

- a list of dimension vectors representing the canonical decomposition of `d`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> canonical_decomposition(Q, [6, 7]) == [[6, 7]]
true

julia> canonical_decomposition(Q, [1, 1]) == [[1, 1]]
true

julia> canonical_decomposition(Q, [6, 2]) == [[3, 1], [3, 1]]
true

julia> Q = kronecker_quiver(2);

julia> canonical_decomposition(Q, [8, 8]) == [[1, 1] for i in 1:8]
true
```
"""
function canonical_decomposition(Q::Quiver, d::AbstractVector{Int})
  general_subdimension_vectors = filter(e -> e != d, all_general_subdimension_vectors(Q, d))
  for e in general_subdimension_vectors
    if d - e in general_subdimension_vectors
      return vcat(canonical_decomposition(Q, e), canonical_decomposition(Q, d - e))
    end
  end
  return [d] # if nothing above worked then d is a Schur root.
end

"""
    in_fundamental_domain(Q::Quiver, d::AbstractVector{Int}; interior::Bool=false)

Check if the dimension vector `d` is in the fundamental domain of the quiver `Q`.

The fundamental domain is the cone of dimension vectors in ``\\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i,

```math
(s_i, d) := \\langle d, s_i\\rangle + \\langle s_i, d\\rangle  \\leq 0,
```

where ``s_i`` is the dimension vector with all entries set to ``0`` and the i-th
set to ``1``.

# Input

- `Q::Quiver` a quiver.
- `d::AbstractVector{Int}` a dimension vector.

Keyword arguments:

- `interior`: if `true` checks whether `d` belongs to the interior of the fundamental
domain. Default is `false`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> in_fundamental_domain(Q, [2, 3])
true

julia> in_fundamental_domain(Q, [1, 1])
true

julia> in_fundamental_domain(Q, [2, 2])
true

julia> in_fundamental_domain(Q, [1, 2])
false
```
"""
function in_fundamental_domain(Q::Quiver, d::AbstractVector{Int}; interior::Bool=false)
  # https://arxiv.org/abs/2209.14791 uses a strict inequality,
  # while https://arxiv.org/abs/2310.15927 uses a non-strict.
  # here we set it to non-strict by default.

  n = n_vertices(Q)
  simples = [unit_vector(n, i) for i in 1:n]
  bound = interior ? -1 : 0
  return all(
    simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) <= bound, simples
  )
end

########################################################################################
# Bocklandt's reduction algorithm
########################################################################################

# One pass of the reduction steps R_I, R_II, R_III of [MR1929191] on a strongly
# connected quiver setting. Returns the new setting, or `nothing` if no step applies,
# i.e., if the setting is reduced in the sense of [Definition 3.1, MR1929191].
#
# The setting is given by a plain adjacency matrix and dimension vector rather than a
# Quiver: the adjacency of a Quiver is an immutable static matrix whose size is a type
# parameter, so the repeated resizing done here would allocate a new type at every
# step, and calling the memoized euler_form on such throwaway quivers would pollute
# its cache.
function __bocklandt_step(A::Matrix{Int}, d::Vector{Int})
  n = length(d)
  # \chi(d, e_v) and \chi(e_v, d), for e_v the unit vector at the vertex v
  chi_in(v) = d[v] - sum(d[w] * A[w, v] for w in 1:n)
  chi_out(v) = d[v] - sum(A[v, w] * d[w] for w in 1:n)
  for v in 1:n
    # R_I [Lemma 3.2, MR1929191]: remove a loopless vertex whose incoming or outgoing
    # paths carry at most d[v] dimensions, shortcutting every path through it; a lone
    # vertex is kept so that the reduced coregular settings are the three settings of
    # [Theorem 1.1, MR1929191]
    if A[v, v] == 0 && n > 1 && (chi_in(v) >= 0 || chi_out(v) >= 0)
      keep = setdiff(1:n, v)
      return A[keep, keep] + A[keep, v] * A[v, keep]', d[keep]
    end
    # R_II [Lemma 3.3, MR1929191]: remove all loops on a vertex of dimension 1
    if A[v, v] > 0 && d[v] == 1
      B = copy(A)
      B[v, v] = 0
      return B, d
    end
    # R_III [Lemma 3.4, MR1929191]: on a vertex of dimension k >= 2 carrying a single
    # loop and, besides the loop, a single incoming (resp. outgoing) arrow from
    # (resp. to) a vertex of dimension 1, remove the loop and thicken that arrow to
    # k parallel arrows
    if A[v, v] == 1 && d[v] >= 2 && (chi_in(v) == -1 || chi_out(v) == -1)
      B = copy(A)
      B[v, v] = 0
      if chi_in(v) == -1
        u = findfirst(w -> w != v && A[w, v] > 0, 1:n)
        B[u, v] = d[v]
      else
        u = findfirst(w -> w != v && A[v, w] > 0, 1:n)
        B[v, u] = d[v]
      end
      return B, d
    end
  end
  return nothing
end

# fully reduce a strongly connected quiver setting, i.e., apply reduction steps until
# the setting is reduced in the sense of [Definition 3.1, MR1929191]
function __bocklandt_reduce(A::Matrix{Int}, d::Vector{Int})
  while (step = __bocklandt_step(A, d)) !== nothing
    A, d = step
  end
  return A, d
end

"""
    bocklandt_reduction(Q::Quiver, d::AbstractVector{Int})

Reduce the quiver setting `(Q, d)` using the reduction steps of
[[MR1929191](https://mathscinet.ams.org/mathscinet/relay-station?mr=1929191)].

The ring of invariants of a quiver setting is the tensor product of those of its
strongly connected components, and vertices of dimension `0` do not contribute, so
these are discarded first, by [Lemma 2.4, MR1929191]. Each component is then
simplified using the three reduction steps of [Section 3, MR1929191], each of which
preserves the ring of invariants up to a polynomial factor:

- ``R_I`` [Lemma 3.2, MR1929191]: a vertex ``v`` without loops with
  ``\\chi(d, e_v) \\geq 0`` or ``\\chi(e_v, d) \\geq 0`` is removed, and every pair of
  arrows ``u \\to v \\to w`` is replaced by an arrow ``u \\to w``;
- ``R_{II}`` [Lemma 3.3, MR1929191]: the loops on a vertex of dimension `1` are
  removed;
- ``R_{III}`` [Lemma 3.4, MR1929191]: the unique loop on a vertex ``v`` of dimension
  ``k \\geq 2`` with ``\\chi(d, e_v) = -1`` (resp. ``\\chi(e_v, d) = -1``) is removed,
  and the unique incoming (resp. outgoing) non-loop arrow is replaced by ``k``
  parallel arrows.

The result, to which no further reduction step applies, is *reduced* in the sense of
[Definition 3.1, MR1929191]; it is returned as the disjoint union of the reduced
components. By [Theorem 3.5, MR1929191] the input setting is coregular if and only if
the reduced setting is, which is what [`is_coregular`](@ref) exploits.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.

# Output

- a dictionary with the reduced quiver `Q` and dimension vector `d`.

# Examples

The setting below is reduced by applying ``R_{III}``, ``R_I``, and ``R_{II}``,
in this order:

```jldoctest
julia> Q = Quiver("1-2, 2-2, 2-1");

julia> setting = bocklandt_reduction(Q, [1, 2]);

julia> setting["Q"]
Quiver with adjacency matrix [0;;]

julia> setting["d"]
1-element Vector{Int64}:
 1
```

A reduced setting is returned unchanged:

```jldoctest
julia> setting = bocklandt_reduction(Quiver("1--2, 2--1"), [1, 1]);

julia> setting["Q"]
Quiver with adjacency matrix [0 2; 2 0]

julia> setting["d"]
2-element Vector{Int64}:
 1
 1
```
"""
function bocklandt_reduction(Q::Quiver, d::AbstractVector{Int})
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("dimension vector must have length $(n_vertices(Q))"))
  all(di >= 0 for di in d) ||
    throw(ArgumentError("dimension vector must be non-negative"))

  # vertices of dimension 0 and arrows between different strongly connected components
  # play no role in the invariant theory [Lemma 2.4, MR1929191]
  A = Matrix{Int}(Q.adjacency)
  vertices = support(d)
  components = strongly_connected_components(Quiver(A[vertices, vertices]))

  reduced = [
    __bocklandt_reduce(A[vertices[c], vertices[c]], Vector{Int}(d[vertices[c]])) for
    c in components
  ]
  return Dict(
    "Q" => reduce(
      disjoint_union, [Quiver(B) for (B, _) in reduced]; init=Quiver(zeros(Int, 0, 0))
    ),
    "d" => reduce(vcat, [e for (_, e) in reduced]; init=Int[]),
  )
end

"""
    is_coregular(Q::Quiver, d::AbstractVector{Int})

Check whether the quiver setting `(Q, d)` is coregular, i.e., whether the ring of
invariants of the `d`-dimensional representation variety of `Q` is a polynomial ring.

Equivalently, this checks whether the affine quotient variety parametrizing
`d`-dimensional semisimple representations of `Q` is smooth, in which case it is an
affine space, by [Theorem 2.1, MR1929191].

By [[Theorem 1.1, MR1929191](https://mathscinet.ams.org/mathscinet/relay-station?mr=1929191)]
this is the case if and only if every strongly connected component of the
[`bocklandt_reduction`](@ref) of `(Q, d)` is one of

- a vertex without loops,
- a vertex with one loop,
- a vertex of dimension `2` with two loops.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.

# Output

- whether the ring of invariants of the setting `(Q, d)` is a polynomial ring.

# Examples

The invariants of pairs of ``2 \\times 2`` matrices form a polynomial ring, but those
of pairs of ``3 \\times 3`` matrices do not, by
[[Procesi](https://mathscinet.ams.org/mathscinet/relay-station?mr=419491)];
the former is the third reduced coregular setting of [Theorem 1.1, MR1929191]:

```jldoctest
julia> is_coregular(jordan_quiver(2), [2])
true

julia> is_coregular(jordan_quiver(2), [3])
false
```

For an acyclic quiver the quotient variety is a point, so the setting is coregular:

```jldoctest
julia> is_coregular(kronecker_quiver(3), [2, 3])
true
```
"""
function is_coregular(Q::Quiver, d::AbstractVector{Int})
  setting = bocklandt_reduction(Q, d)
  A, e = setting["Q"].adjacency, setting["d"]
  return all(
    length(c) == 1 && (A[c[1], c[1]] <= 1 || (A[c[1], c[1]], e[c[1]]) == (2, 2)) for
    c in strongly_connected_components(setting["Q"])
  )
end

########################################################################################
# Cofree quiver settings
########################################################################################

# Everything below implements the classification of cofree quiver settings of
# Bocklandt--Van de Weyer [doi:10.1016/j.jalgebra.2007.08.019]: wedge away vertices
# using their reduction step W, split into prime components, and compare against the
# list of Theorem 1, whose members are recognized by the criteria of Theorems 5, 6, 8
# and 9. Paths and cycles are quasiprimitive throughout: they use every vertex w as a
# source of at most d[w] arrows.

# The number of quasiprimitive cycles through v, counted with arrow multiplicities and
# capped to bound the enumeration; only used when every such cycle passes through v
# exactly once, so that counting closed walks anchored at v is correct.
function __n_quasiprimitive_cycles(A::Matrix{Int}, d::Vector{Int}, v::Int, cap::Int)
  total, budget = Ref(0), copy(d)
  function walk(x::Int, mult::Int)
    (total[] > cap || budget[x] == 0) && return nothing
    budget[x] -= 1
    for y in findall(>(0), A[x, :])
      y == v ? (total[] += mult * A[x, y]) : walk(y, mult * A[x, y])
    end
    budget[x] += 1
    return nothing
  end
  walk(v, 1)
  return total[]
end

# One application of the wedging step W of [doi:10.1016/j.jalgebra.2007.08.019] to a
# vertex v of dimension at least 2 whose unique outgoing (resp. incoming) arrow ends
# (resp. starts) at a vertex of dimension 1: v is removed and its other arrows are
# redirected to that vertex, provided d[v] is at least the number of quasiprimitive
# cycles through v. Wedging preserves cofreeness in both directions [Lemma 3].
function __wedge_step(A::Matrix{Int}, d::Vector{Int})
  n = length(d)
  for v in findall(v -> d[v] >= 2 && A[v, v] == 0, 1:n)
    outs, ins = findall(>(0), A[v, :]), findall(>(0), A[:, v])
    out = length(outs) == 1 && A[v, outs[1]] == 1 && d[outs[1]] == 1
    into = length(ins) == 1 && A[ins[1], v] == 1 && d[ins[1]] == 1
    ((out || into) && __n_quasiprimitive_cycles(A, d, v, d[v]) <= d[v]) || continue
    B = copy(A)
    out ? (B[:, outs[1]] .+= A[:, v]) : (B[ins[1], :] .+= A[v, :])
    keep = setdiff(1:n, v)
    return B[keep, keep], d[keep]
  end
  return nothing
end

# Split a strongly connected quiver setting into its prime components, i.e., the
# summands of its decomposition as an iterated connected sum at vertices of
# dimension 1; a setting is cofree iff its prime components are [Lemma 3].
function __prime_components(A::Matrix{Int}, d::Vector{Int})
  n = length(d)
  for v in findall(==(1), d)
    # the summands at v are the weakly connected components of the quiver minus v,
    # each taken together with v and the arrows between them, and every loop at v
    others = setdiff(1:n, v)
    U = A[others, others]
    pieces = strongly_connected_components(Quiver(U + U'))
    length(pieces) + A[v, v] >= 2 || continue
    out = [(fill(1, 1, 1), [1]) for _ in 1:A[v, v]]
    for piece in pieces
      keep = sort!(vcat(others[piece], v))
      B = A[keep, keep]
      w = findfirst(==(v), keep)
      B[w, w] = 0
      append!(out, __prime_components(B, d[keep]))
    end
    return out
  end
  return [(A, d)]
end

# Decide cofreeness of a prime strongly connected setting by recognizing the members
# of the list of [Theorem 1, doi:10.1016/j.jalgebra.2007.08.019].
function __is_cofree_prime(A::Matrix{Int}, d::Vector{Int})
  n = length(d)
  # a single vertex: no arrows, a cyclic quiver (one loop), any number of loops on a
  # vertex of dimension 1, or the setting Q_2 (two loops on a vertex of dimension 2)
  n == 1 && return A[1, 1] <= 1 || d[1] == 1 || (A[1, 1], d[1]) == (2, 2)

  # (iii) cyclic quiver settings are always cofree [Theorem 5]
  ins, outs = vec(sum(A; dims=1)), vec(sum(A; dims=2))
  all(ins .== 1) && all(outs .== 1) && return true

  # (i) all cycles run through a vertex v of dimension 1 [Theorem 6]: cofree iff
  # d[w] >= #{quasiprimitive paths v -> w} + #{quasiprimitive paths w -> v} - 1 for
  # all other w; the quiver minus v is acyclic, so its adjacency powers count paths
  for v in findall(==(1), d)
    B = copy(A)
    B[v, :] .= 0
    B[:, v] .= 0
    any(!=(0), B^n) && continue
    S = sum(B^k for k in 0:(n - 1))
    return all(
      d[w] >= A[v, :]' * S[:, w] + S[w, :]' * A[:, v] - 1 for w in 1:n if w != v
    )
  end

  # (ii) and (iv) are two cycles sharing a path of s >= 1 vertices: n + 1 arrows, a
  # unique vertex x of out-degree 2, a unique y of in-degree 2, all other degrees 1;
  # the shared path runs from y to x, the two branches lead from x back to y
  sum(outs) == n + 1 || return false
  x, y = findfirst(==(2), outs), findfirst(==(2), ins)
  (isnothing(x) || isnothing(y)) && return false
  shared = [y]
  while shared[end] != x
    length(shared) > n && return false
    push!(shared, findfirst(>(0), A[shared[end], :]))
  end
  function branch(cur::Int)
    b = Int[]
    while cur != y
      (cur == x || cur in shared || cur in b || length(b) > n) && return nothing
      push!(b, cur)
      cur = findfirst(>(0), A[cur, :])
    end
    return b
  end
  targets = findall(>(0), A[x, :])
  b1, b2 = branch(targets[1]), branch(targets[end])
  (isnothing(b1) || isnothing(b2) || length(shared) + length(b1) + length(b2) != n) &&
    return false

  # (ii) a branch is a single vertex of dimension 1: cofree iff the minimal dimension
  # along the other cycle is attained exactly once in the shared path, or not there
  # but exactly once in the other branch [Theorem 8]
  for (c, rest) in ((b1, b2), (b2, b1))
    if length(c) == 1 && d[c[1]] == 1
      m = minimum(d[vcat(shared, rest)])
      return count(==(m), d[shared]) == 1 ||
             (count(==(m), d[shared]) == 0 && count(==(m), d[rest]) == 1)
    end
  end

  # (iv) all branch dimensions at least 2, exactly one shared dimension equal to 2,
  # and the other shared dimensions at least 4 [Theorem 9]
  return all(d[vcat(b1, b2)] .>= 2) &&
         count(==(2), d[shared]) == 1 &&
         all(w -> w == 2 || w >= 4, d[shared])
end

"""
    is_cofree(Q::Quiver, d::AbstractVector{Int})

Check whether the quiver setting `(Q, d)` is cofree, i.e., whether the coordinate ring
of the `d`-dimensional representation variety of `Q` is a graded free module over its
ring of invariants.

By a criterion of Popov this is the case if and only if the setting is coregular (see
[`is_coregular`](@ref)) and its nullcone is equidimensional. The implementation follows
the classification of
[[Bocklandt--Van de Weyer](https://doi.org/10.1016/j.jalgebra.2007.08.019)]:
the setting is cofree if and only if all its strongly connected components are, which
is decided by wedging away vertices (their reduction step ``W``), splitting into prime
components (the summands of the decomposition as an iterated connected sum at vertices
of dimension `1`), and comparing against the list of [Theorem 1, loc. cit.].

Cofreeness is stronger than coregularity: it moreover makes the quotient map from the
representation variety to the affine quotient flat.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.

# Output

- whether the coordinate ring of the setting `(Q, d)` is a graded free module over the
  ring of invariants.

# Examples

Pairs of ``2 \\times 2`` matrices are cofree, pairs of ``3 \\times 3`` matrices are
not even coregular, and cyclic quiver settings are always cofree:

```jldoctest
julia> is_cofree(jordan_quiver(2), [2])
true

julia> is_cofree(jordan_quiver(2), [3])
false

julia> is_cofree(cyclic_quiver(3), [1, 2, 3])
true
```

A coregular setting need not be cofree:

```jldoctest
julia> Q = Quiver("1--2, 2-1");

julia> is_coregular(Q, [2, 2]), is_cofree(Q, [2, 2])
(true, false)

julia> is_coregular(Q, [2, 4]), is_cofree(Q, [2, 4])
(true, true)
```
"""
function is_cofree(Q::Quiver, d::AbstractVector{Int})
  length(d) == n_vertices(Q) ||
    throw(ArgumentError("dimension vector must have length $(n_vertices(Q))"))
  all(di >= 0 for di in d) ||
    throw(ArgumentError("dimension vector must be non-negative"))

  # vertices of dimension 0 do not contribute, and arrows between different strongly
  # connected components only contribute a free matrix factor [Lemma 3]
  A = Matrix{Int}(Q.adjacency)
  vertices = support(d)
  for c in strongly_connected_components(Quiver(A[vertices, vertices]))
    Ac, dc = A[vertices[c], vertices[c]], Vector{Int}(d[vertices[c]])
    # wedge the vertices of dimension at least 2 first, then split into prime
    # components; wedges at vertices of dimension 1 only occur for cyclic quivers,
    # which are cofree anyway [Remark 2]
    while (step = __wedge_step(Ac, dc)) !== nothing
      Ac, dc = step
    end
    all(__is_cofree_prime(B, e) for (B, e) in __prime_components(Ac, dc)) ||
      return false
  end
  return true
end
