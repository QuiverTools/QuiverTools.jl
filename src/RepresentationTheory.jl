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
