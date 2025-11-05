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
  x' * euler_matrix(Q) * y

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

  simples = [unit_vector(n_vertices(Q), i) for i in 1:n_vertices(Q)]
  if interior
    return all(
      simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) < 0,
      simples,
    )
  end
  return all(simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) <= 0, simples)
end
