########################################################################################
# Below lie methods to compute Hodge diamonds translated from the Hodge diamond cutter.
# In turn, these are based on M. Reineke's paper
# "The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli",
# https://doi.org/10.1007/s00222-002-0273-4
########################################################################################

###################################################
# auxiliary functions for hodge_polynomial() below
"""
    solve(A, b)

Solve ``A\\cdot x = b`` for ``A`` upper triangular via back substitution.

This is an internal method only used in the implementation of the Hodge polynomial
and to compute motives.


# Input

- `A::AbstractMatrix`: an upper triangular matrix.
- `b::AbstractVector`: a vector.

# Output

- the solution `x` to the equation.

# Examples

```jldoctest
julia> A = [1 2 3; 0 4 5; 0 0 6];

julia> b = [1, 2, 3];

julia> QuiverTools.solve(A, b)
3-element Vector{Any}:
 -0.25
 -0.125
  0.5
```
"""
function solve(A, b)
  n = length(b)
  x = Vector{Any}(zeros(n))

  x[n] = b[n] / A[n, n]

  for i in (n - 1):-1:1
    x[i] = (b[i] - sum(A[i, j] * x[j] for j in (i + 1):n)) / A[i, i]
  end
  return x
end

"""
Cardinality of general linear group ``\\mathrm{GL}_n(\\mathbb{F}_v)``.
"""
@memoize Dict function CardinalGl(n::Int, q)
  if n == 0
    return 1
  else
    out = q^n - 1
    for i in 1:(n - 1)
      out *= q^n - q^i
    end
    return out
  end
end

"""
Cardinality of representation space ``\\mathrm{R}(Q,d), over \\mathbb{F}_q``.
"""
function CardinalRd(Q::Quiver, d::AbstractVector{Int}, q)
  return q^sum(
    d[i] * d[j] * Q.adjacency[i, j] for i in 1:n_vertices(Q), j in 1:n_vertices(Q); init=0
  )
end

"""
Cardinality of product of general linear groups ``\\mathrm{GL}_{d}(\\mathbb{F}_q)``.
"""
@memoize Dict function CardinalGd(d::AbstractVector{Int}, q)
  return prod(CardinalGl(di, q) for di in d)
end

"""
Entry of the transfer matrix, as per Corollary 6.9 of
[[MR1974891](https://doi.org/10.1007/s00222-002-0273-4)].
"""
function TransferMatrixEntry(Q, e, f, q)
  fe = f - e

  if all(fei >= 0 for fei in fe)
    return q^euler_form(Q, -fe, e) * CardinalRd(Q, fe, q) / CardinalGd(fe, q)
  else
    return 0
  end
end

function Td(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, q)
  # indexing set for the transfer matrix
  I = filter(
    e -> slope(e, theta) > slope(d, theta),
    all_subdimension_vectors(d; nonzero=true, strict=true),
  )
  I = vcat([zero_vector(Q)], I, [d])

  l = length(I)
  T = Matrix{Any}(zeros(l, l))

  for (i, Ii) in enumerate(I)
    for j in i:l  # upper triangular
      T[i, j] = TransferMatrixEntry(Q, Ii, I[j], q)
    end
  end
  return T
end

# auxiliary functions for hodge_polynomial() above
###################################################

"""
    hodge_polynomial(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d))

Return the Hodge polynomial of the moduli space of `theta`-semistable
representations of `Q` with dimension vector `d`.

The algorithm is based on [[]MR1974891](https://doi.org/10.1007/s00222-002-0273-4)],
and the current implementation is translated from the [Hodge diamond cutter]
(https://zenodo.org/doi/10.5281/zenodo.3893509).

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is `canonical_stability(Q, d)`.

# Output

- the Hodge polynomial of the moduli space.

# Examples

The Hodge polynomial of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> d = [2, 3];

julia> theta = [3, -2];

julia> hodge_polynomial(Q, d, theta)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1
```

Cases on a fake wall:
```jldoctest
julia> Q = Quiver([0 1 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]); d = [3, 3, 4, 1];

julia> hodge_polynomial(Q, d)
x^3*y^3 + 3*x^2*y^2 + 3*x*y + 1

julia> Q = Quiver([0 1 0 1 0; 0 0 1 0 2; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]); d = [1, 2, 1, 1, 1];

julia> hodge_polynomial(Q, d)
x^3*y^3 + 4*x^2*y^2 + 4*x*y + 1

"""
function hodge_polynomial(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
)

  # safety checks
  !is_acyclic(Q) && throw(ArgumentError("Q is not acyclic."))
  has_properly_semistables(Q, d, theta) && throw(
    ArgumentError(
      "The quiver moduli problem has properly semistable representations, no description of the Hodge polynomial is known."
    ),
  )

  R, q = polynomial_ring(Singular.QQ, ["q"])
  F = fraction_field(R)

  v = F(q[1]) # worsens performance by ~8%. Necessary?

  T = Td(Q, d, theta, v)

  one_at_the_end = unit_vector(size(T, 1), size(T, 1))

  # @warn "result needs to be a polynomial, otherwise the moduli space is singular."
  solution = solve(T, one_at_the_end)[1] * (1 - v)
  denominator(solution) != 1 && throw(DomainError("Moduli space is singular!"))
  result = numerator(solution)

  S, (x, y) = polynomial_ring(Singular.QQ, ["x", "y"])
  return S(result(x * y))
end

"""
    hodge_polynomial(M::QuiverModuliSpace)

Compute the Hodge polynomial of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Hodge polynomial of the moduli space.

# Examples

The Hodge polynomial of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> hodge_polynomial(M)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1
```
"""
function hodge_polynomial(M::QuiverModuliSpace)
  return hodge_polynomial(M.Q, M.d, M.theta)
end

"""
    hodge_diamond(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d))

Compute the Hodge diamond of the moduli space of
`theta`-semistable representations of `Q` with dimension vector `d`.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is the canonical stability.

# Output

- the Hodge diamond of the moduli space.

# Examples

The Hodge diamond of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> hodge_diamond(Q, [2, 3])
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```

This method correctly handles the moduli spaces being empty or 0-dimensional:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> hodge_diamond(Q, [2, 3], [-3, 2])
0×0 Matrix{Int64}
```
"""
function hodge_diamond(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
)
  g = hodge_polynomial(Q, d, theta)

  # collects the coefficients of the polynomial, converts them to integers
  # and returns them in the diagonal of a matrix.
  return Matrix{Int}(diagonal(Int.(numerator.(collect(Singular.coefficients(g))))))
end
"""
    hodge_diamond(M::QuiverModuliSpace)

Compute the Hodge diamond of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Hodge diamond of the moduli space.

# Examples

The Hodge diamond of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> hodge_diamond(M)
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```
"""
function hodge_diamond(M::QuiverModuliSpace)
  return hodge_diamond(M.Q, M.d, M.theta)
end

"""
    picard_rank(M::QuiverModuliSpace)

Compute the Picard rank of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Picard rank of the moduli space.

# Examples

Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> picard_rank(M)
1
```
"""
function picard_rank(M::QuiverModuliSpace)
  !(is_smooth(M) && is_projective(M)) &&
    throw(ArgumentError("Moduli space is not smooth and projective"))
  return betti_numbers(M)[3]
end

"""
    index(M::QuiverModuliSpace)

Compute the index of the moduli space `M`.

The index of a variety ``X`` is the largest integer which divides
the canonical divisor ``K_X`` in ``Pic(X)``.

This implementation currently only works for the canonical stability.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the index of the moduli space.

# Examples

The 3-Kronecker quiver has index 3:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> index(M)
3
```
The subspace quiver moduli have index 1:
```jldoctest
julia> Q = subspace_quiver(5);

julia> M = QuiverModuliSpace(Q, [1, 1, 1, 1, 1, 2]);

julia> index(M)
1
```
"""
function index(M::QuiverModuliSpace)
  has_properly_semistables(M.Q, M.d, M.theta, M.denom) && throw(
    ArgumentError(
      "The quiver moduli problem has properly semistable representations, no description of the Mukai index is known."
    ),
  )
  !is_amply_stable(M) && throw(
    ArgumentError(
      "The quiver moduli problem is not amply stable, no description of the Mukai index is known."
    ),
  )
  # The index formula of [MR4352662] only holds for a full-support dimension vector, so
  # restrict to the support first: `M(Q, d)` is isomorphic to the moduli space of the
  # support subquiver, but zero entries of `d` would otherwise pollute the gcd (issue #12).
  Qs, ds = support_subquiver(M.Q, M.d)
  return gcd(canonical_stability(Qs, ds))
end

"""
    betti_numbers(M::QuiverModuliSpace)

Compute the Betti numbers of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- a list of Betti numbers of the moduli space, indexed by cohomological
  degree ``0, \\dots, 2\\dim M``.

# Examples

```jldoctest
julia> Q = kronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> betti_numbers(M)
3-element Vector{Int64}:
 1
 0
 1
```

Our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> betti_numbers(M)
13-element Vector{Int64}:
 1
 0
 1
 0
 3
 0
 3
 0
 3
 0
 1
 0
 1
```
"""
function betti_numbers(M::QuiverModuliSpace)
  !is_coprime(M.d, M.theta) && throw(ArgumentError("d and theta are not coprime"))

  N = dimension(M)
  P = poincare_polynomial(M)
  # entry 2k + 1 is the coefficient of L^k in P, i.e., the Betti number
  # in cohomological degree 2k; all odd Betti numbers vanish
  betti = zeros(Int, 2 * N + 1)
  for (c, e) in zip(Singular.coefficients(P), Singular.exponent_vectors(P))
    betti[2 * e[1] + 1] = Int(numerator(c))
  end
  return betti
end

"""
    poincare_polynomial(M::QuiverModuliSpace)

Compute the Poincaré polynomial of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Poincaré polynomial of the moduli space.

# Examples

A Kronecker quiver setup where `M` is the projective line:
```jldoctest
julia> Q = kronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> poincare_polynomial(M)
L + 1
```

The Poincaré polynomial of our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> poincare_polynomial(M)
L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1
```
"""
function poincare_polynomial(M::QuiverModuliSpace)
  !is_coprime(M.d, M.theta) && throw(ArgumentError("d and theta are not coprime"))

  m = motive(M.Q, M.d, M.theta, M.denom)
  v = Singular.transcendence_basis(Singular.parent(m))[1]
  # the stack of stable representations is a G_m-gerbe over the moduli space, so
  # [M] = (L - 1) * [stack motive]
  P = (v - 1) * m

  denominator(P) != 1 && throw(DomainError("must be a polynomial"))
  # returns a polynomial object instead of a FunctionField element.
  return Singular.n_transExt_to_spoly(numerator(P))
end

function power(x, n::Int)
  if n >= 0
    return x^n
  else
    return 1 / x^(-n)
  end
end

function motive(M::QuiverModuliStack)
  M.condition != "stable" && throw(ArgumentError("Motive unknown if not stable"))
  return motive(M.Q, M.d, M.theta)
end

"""
    motive(Q::Quiver, d, theta, denom=sum)

Compute the motive of the moduli stack of `theta`-semistable representations.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is the canonical stability.
- `denom::Function`: a function. Default is the sum.

# Output

- The motive as an element in the function field \\mathbb{Q}(L).

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> motive(Q, [2, 3])
(L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1)//(L - 1)
```
"""
function motive(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum,
)
  K, L = Singular.FunctionField(Singular.QQ, ["L"])
  L = L[1]

  if all(ti == 0 for ti in theta)
    out = power(L, -euler_form(Q, d, d))
    den = 1
    for i in support(d)
      den *= prod(1 - power(L, -nu) for nu in 1:d[i])
    end
    return out / den
  end

  ds = all_destabilizing_subdimension_vectors(d, theta, denom)

  push!(ds, zero_vector(Q), d)
  sort!(ds; by=e -> deglex_key(Q, e))

  T = Matrix{Any}(undef, length(ds), length(ds))
  for (i, j) in Iterators.product(1:length(ds), 1:length(ds))
    if is_subdimension_vector(ds[i], ds[j])
      T[i, j] =
        power(L, euler_form(Q, ds[i] - ds[j], ds[i])) *
        motive(Q, ds[j] - ds[i], zero_vector(Q))
    else
      T[i, j] = 0
    end
  end

  y = [0 for i in 1:length(ds)]
  y[end] = 1
  y = coerce_vector(y)

  # the HN recursion produces the negative of the honest stack motive; negate so
  # this branch agrees with the trivial-stability branch above (see issue #36)
  return -solve(T, y)[1]
end
