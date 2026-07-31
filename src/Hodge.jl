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

########################################################################################
# Intersection cohomology of quiver moduli spaces.
#
# When `d` and `theta` admit properly semistable representations the moduli space is
# singular and its ordinary cohomology is not determined by the Harder--Narasimhan
# recursion. Meinhardt--Reineke
# [[MR4000572](https://mathscinet.ams.org/mathscinet-getitem?mr=4000572)] identify the
# Donaldson--Thomas invariants of the quiver with the intersection cohomology of the
# moduli space, which makes it computable from the motives of the semistable stacks.
#
# With ``\Lambda_\mu`` the monoid of dimension vectors of the slope of ``d``, their
# Lemma in §3.1 and Theorem 3.4 read
#
#     \sum_{e \in \Lambda_\mu} L^{(e, e)/2} [\mathfrak{M}^{ss}_e] t^e
#       = Exp((\sum_{0 \neq e \in \Lambda_\mu} DT_e t^e)/(L^{1/2} - L^{-1/2})),
#     E(IH^*(M^{ss}_d)) = L^{\dim/2} DT_d,   \dim M^{ss}_d = 1 - (d, d),
#
# where ``(-, -)`` is the Euler form and ``Exp`` is the plethystic exponential.
# The motive of the stack ``\mathfrak{M}^{ss}_e`` is what `motive` computes.
#
# Quiver moduli have Hodge structures concentrated on the diagonal, so everything in
# sight is a rational function in the Lefschetz class alone, and the half powers only
# need a square root `w` of it, with ``L^{1/2} = -w`` because ``L^{1/2}`` sits in odd
# degree. The Adams operations are then the substitutions ``w \mapsto w^n``.
#
# For a dimension vector which is primitive in ``\Lambda_\mu`` the plethystic logarithm
# is its own leading term, and the answer is the ordinary Poincaré polynomial again.
########################################################################################

"""
    _mobius(n::Int)

Return the Möbius function ``\\mu(n)``.

This is an internal method, only used in the plethystic logarithm computing
intersection cohomology.

# Examples

```jldoctest
julia> QuiverTools._mobius.(1:10)
10-element Vector{Int64}:
  1
 -1
 -1
  0
 -1
  1
 -1
  0
  0
  1
```
"""
function _mobius(n::Int)
  result = 1
  for p in 2:n            # the range is fixed before the loop divides `n` down
    if n % p == 0
      n = n ÷ p
      n % p == 0 && return 0
      result = -result
    end
  end
  return result
end

"""
    intersection_poincare_polynomial(M::QuiverModuliSpace)

Compute the Poincaré polynomial of the intersection cohomology of the moduli space `M`.

The algorithm is the one of
[[MR4000572](https://mathscinet.ams.org/mathscinet-getitem?mr=4000572)], which identifies
the Donaldson--Thomas invariants of the quiver with the intersection cohomology of
``M^{ss}_\\theta(Q, \\mathbf{d})``. It needs `M.theta` to be generic for the slope of
`M.d`, meaning that the antisymmetrized Euler form vanishes on the dimension vectors of
that slope, and it needs stable representations to exist; both are checked.

The moduli space is smooth exactly when no proper subdimension vector has the slope of
`M.d`, and there intersection cohomology is ordinary cohomology, so this agrees with
[`poincare_polynomial`](@ref).

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Poincaré polynomial of the intersection cohomology of the moduli space.

# Examples

The moduli space for the 3-Kronecker quiver and dimension vector `[2, 2]` is singular,
as `[1, 1]` has the same slope, and it has the intersection cohomology of ``\\mathbb{P}^5``:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 2]);

julia> intersection_poincare_polynomial(M)
L^5 + L^4 + L^3 + L^2 + L + 1
```

Reflection functors identify moduli spaces for different dimension vectors:
```jldoctest
julia> Q = kronecker_quiver(4);

julia> M = QuiverModuliSpace(Q, [3, 3]); N = QuiverModuliSpace(Q, [3, 9]);

julia> intersection_poincare_polynomial(M) == intersection_poincare_polynomial(N)
true
```

In the smooth case this is the ordinary Poincaré polynomial:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> intersection_poincare_polynomial(M)
L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1

julia> intersection_poincare_polynomial(M) == poincare_polynomial(M)
true
```

There is nothing to compute if no representation of dimension vector `M.d` is stable:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(2), [2, 2]);

julia> intersection_poincare_polynomial(M)
ERROR: ArgumentError: there are no stable representations of dimension vector [2, 2]
```
"""
function intersection_poincare_polynomial(M::QuiverModuliSpace)
  M.condition == "semistable" || throw(
    ArgumentError("intersection cohomology is computed for the semistable moduli space")
  )
  Q, d, theta, denom = M.Q, M.d, M.theta, M.denom

  # the dimension vectors of the slope of `d`, i.e. the part of the monoid
  # ``\Lambda_\mu`` below `d`; the zero vector is left out and treated separately
  lattice = filter(
    e -> slope(e, theta, denom) == slope(d, theta, denom),
    all_subdimension_vectors(d; nonzero=true),
  )
  all(euler_form(Q, e, f) == euler_form(Q, f, e) for e in lattice, f in lattice) || throw(
    ArgumentError(
      "the stability parameter is not generic for the slope of $(Vector(d)), " *
      "so intersection cohomology is out of reach",
    ),
  )

  R, ws = polynomial_ring(Singular.QQ, ["w"])
  w = ws[1]
  F = fraction_field(R)
  # the Adams operation ``\psi^n``, and the passage from the Lefschetz class to `w`
  psi(x, n) = F(numerator(x)(w^n))//F(denominator(x)(w^n))
  half(m) =
    F(Singular.n_transExt_to_spoly(numerator(m))(w^2)) //
    F(Singular.n_transExt_to_spoly(denominator(m))(w^2))

  # the generating series, without its constant term
  E = typeof(F(w))
  series = Dict{Vector{Int},E}(
    e => power(F(-w), euler_form(Q, e, e)) * half(motive(Q, e, theta, denom)) for
    e in lattice
  )

  # the ordinary logarithm ``\log(1 + x) = \sum_k (-1)^{k-1}/k x^k``; the kth power is
  # supported on sums of k nonzero dimension vectors, so the sum stops at ``|d|``
  logarithm = Dict{Vector{Int},E}()
  term = series
  for k in 1:sum(d)
    isempty(term) && break
    scale = F((-1)^(k - 1))//F(k)
    for (e, value) in term
      logarithm[e] = get(logarithm, e, zero(F)) + scale * value
    end
    # multiply by `series`, dropping everything that is no longer below `d`
    next = Dict{Vector{Int},E}()
    for (e, value) in term, (f, factor) in series
      is_subdimension_vector(e + f, d) || continue
      next[e + f] = get(next, e + f, zero(F)) + value * factor
    end
    term = next
  end

  # and the plethystic one, ``Log(1 + x) = \sum_n \mu(n)/n \psi^n(\log(1 + x))``; only
  # those `n` with `n * e = d` for some `e` contribute, i.e. the divisors of `gcd(d)`
  total = get(logarithm, Vector{Int}(d), zero(F))
  g = gcd(d)
  for k in 2:g
    (g % k == 0 && _mobius(k) != 0) || continue
    piece = get(logarithm, Vector{Int}(d .÷ k), zero(F))
    total += F(_mobius(k))//F(k) * psi(piece, k)
  end

  # ``DT_d = (L^{1/2} - L^{-1/2}) [Log Q]_{t^d}``, then ``E(IH^*) = L^{\dim/2} DT_d``
  root = F(-w)
  result = power(root, 1 - euler_form(Q, d, d)) * (root - inv(root)) * total
  iszero(result) && throw(
    ArgumentError("there are no stable representations of dimension vector $(Vector(d))")
  )
  isone(denominator(result)) ||
    throw(DomainError("intersection cohomology is not polynomial"))

  # intersection cohomology of these moduli spaces is concentrated in even degree, so the
  # answer has to be a polynomial in `w^2`; that is a real check on the whole computation
  S, Ls = polynomial_ring(Singular.QQ, ["L"])
  L = Ls[1]
  invariant = numerator(result)
  P = zero(S)
  for (c, e) in
      zip(Singular.coefficients(invariant), Singular.exponent_vectors(invariant))
    isodd(e[1]) && throw(DomainError("intersection cohomology in odd degree"))
    P += S(c) * L^(e[1] ÷ 2)
  end
  return P
end

"""
    intersection_betti_numbers(M::QuiverModuliSpace)

Compute the Betti numbers of the intersection cohomology of the moduli space `M`.

See [`intersection_poincare_polynomial`](@ref) for the algorithm and its hypotheses.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- a list of intersection Betti numbers of the moduli space, indexed by cohomological
  degree ``0, \\dots, 2\\dim M``.

# Examples

The singular moduli space for the 3-Kronecker quiver and dimension vector `[2, 2]`:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> intersection_betti_numbers(M)
11-element Vector{Int64}:
 1
 0
 1
 0
 1
 0
 1
 0
 1
 0
 1
```

In the smooth case these are the ordinary Betti numbers:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> intersection_betti_numbers(M) == betti_numbers(M)
true
```
"""
function intersection_betti_numbers(M::QuiverModuliSpace)
  P = intersection_poincare_polynomial(M)
  # entry 2k + 1 is the coefficient of L^k in P, i.e. the intersection Betti number
  # in cohomological degree 2k; all odd ones vanish
  betti = zeros(Int, 2 * dimension(M) + 1)
  for (c, e) in zip(Singular.coefficients(P), Singular.exponent_vectors(P))
    betti[2 * e[1] + 1] = Int(numerator(c))
  end
  return betti
end

"""
    intersection_hodge_diamond(M::QuiverModuliSpace)

Compute the Hodge diamond of the intersection cohomology of the moduli space `M`.

See [`intersection_poincare_polynomial`](@ref) for the algorithm and its hypotheses.
The Hodge structure is of Hodge--Tate type, so the diamond is concentrated on the
diagonal.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Hodge diamond of the intersection cohomology of the moduli space.

# Examples

The singular moduli space for the 3-Kronecker quiver and dimension vector `[2, 2]` has
the intersection cohomology of ``\\mathbb{P}^5``:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> intersection_hodge_diamond(M)
6×6 Matrix{Int64}:
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  1  0  0  0
 0  0  0  1  0  0
 0  0  0  0  1  0
 0  0  0  0  0  1
```

In the smooth case this is the ordinary Hodge diamond:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> intersection_hodge_diamond(M) == hodge_diamond(M)
true
```
"""
function intersection_hodge_diamond(M::QuiverModuliSpace)
  return Matrix{Int}(diagonal(intersection_betti_numbers(M)[1:2:end]))
end
