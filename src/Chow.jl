###############################################################################
# tautological representation of the Chow ring.
# Implements the results of [arXiv:1307.3066](https://doi.org/10.48550/arXiv.1307.3066)
# and [arXiv.2307.01711](https://doi.org/10.48550/arXiv.2307.01711).
###############################################################################

# partial order on the forbidden dimension vectors as defined in
# https://doi.org/10.48550/arXiv.1307.3066
function partial_order(Q::Quiver, f::AbstractVector{Int}, g::AbstractVector{Int})
  if !all(f[i] <= g[i] for i in 1:n_vertices(Q) if is_source(Q, i))
    return false
  elseif !all(f[i] >= g[i] for i in 1:n_vertices(Q) if is_sink(Q, i))
    return false
  elseif !all(f[i] == g[i] for i in 1:n_vertices(Q) if !is_source(Q, i) && !is_sink(Q, i))
    return false
  end
  return true
end

"""
    symmetric_polynomial(vars, degree::Int)

Compute the symmetric polynomial of degree `degree` in the variables `vars`.

# Input

- `vars`: a list of variables.
- `degree`: the degree of the wanted symmetric polynomial.

# Output

- The symmetric polynomial of degree `degree` in the variables `vars`.

# Examples

```julia-repl
julia> using Singular;

julia> R, vars = polynomial_ring(Singular.QQ, ["x", "y", "z"]);

julia> QuiverTools.symmetric_polynomial(vars, 2)
x*y + x*z + y*z
```
"""
function symmetric_polynomial(vars, degree::Int)
  return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

"""
    chow_ring(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d); chi::AbstractVector{Int}=extended_gcd(d)[2])

Compute the Chow ring of the moduli space of `theta`-semistable representations of
`Q` with dimension vector `d`, for a choice of linearization `a`.

This method of the function `chow_ring` also returns the ambient ring ``R``
and the inclusion morphism.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is `canonical_stability(Q, d)`.
- `chi`: a linearization. Default is the extended gcd of `extended_gcd(d)[2]`.

# Output

A tuple containing:
- the Chow ring of the moduli space,
- the polynomial ring above it,
- the inclusion map ``\\iota : A \\to R``.

# Examples

The Chow ring for the projective line has two generators:
```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> CH = chow_ring(M);

julia> QuiverTools.gens(QuiverTools.quotient_ideal(CH))
2-element Vector{Singular.spoly{Singular.n_Q}}:
 x21
 x11^2
```

The Chow ring for our favourite 6-fold has, in this implementation, 16 generators:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> CH = chow_ring(M); I = QuiverTools.quotient_ideal(CH);

julia> length(QuiverTools.gens(I))
16
```
"""
function chow_ring(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d);
  chi::AbstractVector{Int}=extended_gcd(d)[2],
)
  # safety checks
  if !is_coprime(d, theta)
    throw(ArgumentError("d and theta are not coprime"))
  elseif chi' * d != 1
    throw(ArgumentError("``chi`` is not a linearization"))
  end

  # j varies first, then i
  varnames = ["xi$i$j" for i in 1:n_vertices(Q) for j in 1:d[i]]
  R, vars = polynomial_ring(Singular.QQ, varnames)

  # Shorthand to address the variable `xi_{i,j}`.
  function xi(i, j)
    d[i] == 0 && throw(ArgumentError("i is not in the support of d."))
    return vars[sum(d[1:(i - 1)]) + j]
  end

  # This is the naive base that is described in Hans's 2013 paper.
  function base_for_ring()
    bounds = [0:(d[i] - nu) for i in 1:n_vertices(Q) for nu in 1:d[i]]
    lambdas = Iterators.product(bounds...)

    build_elem(lambda) = prod(
      prod(xi(i, nu)^lambda[sum(d[1:(i - 1)]) + nu] for nu in 1:d[i]) for
      i in support(d)
    )
    return map(l -> build_elem(l), lambdas)
  end

  # build the permutation group W
  W = Iterators.product([AbstractAlgebra.SymmetricGroup(d[i]) for i in 1:n_vertices(Q)]...)
  sign(w) = prod(AbstractAlgebra.sign(wi) for wi in w)

  # Action of the symmetric group on R by permutation of the variables.
  permute(f, sigma) =
    f([xi(i, sigma[i][j]) for i in support(d) for j in 1:d[i]]...)

  # The discriminant in the definition of the antisymmetrization.
  delta = 1
  for i in 1:n_vertices(Q)
    d[i] > 1 && (
      delta *= prod(
        xi(i, l) - xi(i, k) for k in 1:(d[i] - 1) for l in (k + 1):d[i]
      )
    )
  end

  antisymmetrize(f) = div(sum(sign(w) * permute(f, w) for w in W), R(delta))

  # All the destabilizing subdimension vectors of `d` with respect to the slope
  # `theta/denom` that are minimal with respect to the total order.
  dest = all_destabilizing_subdimension_vectors(d, theta)
  minimal_forbidden = filter(
    e -> !any(f -> partial_order(Q, f, e), filter(f -> f != e, dest)), dest
  )

  # builds a new forbidden polynomial for the minimal forbidden dimension vector e.
  function new_forbidden(e::AbstractVector{Int})
    out = 1
    for (i, j) in Iterators.product(1:n_vertices(Q), 1:n_vertices(Q))
      for r in 1:e[i], s in (e[j] + 1):d[j]
        out *= (xi(j, s) - xi(i, r))^Q.adjacency[i, j]
      end
    end
    return out
  end
  forbidden_polynomials = [new_forbidden(e) for e in minimal_forbidden]

  varnames2 = ["x$i$j" for i in 1:n_vertices(Q) for j in 1:d[i]]
  A, Avars = polynomial_ring(Singular.QQ, varnames2)

  # Shorthand to address the variables of A `x_{i,j}`.
  function xs(i, j)
    d[i] == 0 && throw(ArgumentError("i is not in the support of d."))
    return Avars[sum(d[1:(i - 1)]) + j]
  end

  targets = [
    [symmetric_polynomial([xi(i, j) for j in 1:d[i]], k) for k in 1:d[i]] for
    i in support(d)
  ]
  targets = reduce(vcat, targets)

  inclusion = AlgebraHomomorphism(A, R, targets)

  anti = [antisymmetrize(f * b) for f in forbidden_polynomials for b in base_for_ring()]
  tautological = [gens(preimage(inclusion, Ideal(R, g)))[1] for g in anti]
  linear = [sum(chi[i] * xs(i, 1) for i in support(d))]

  return (QuotientRing(A, std(Ideal(A, [tautological; linear]))), R, inclusion)
end

"""
    chow_ring(M::QuiverModuliSpace; chi::Union{AbstractVector{Int},UndefInitializer}=undef)

Compute the Chow ring of the moduli space `M` for the given linearization `chi`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `chi::AbstractVector{Int}`: a choice of linearization for the trivial line bundle.
  Default is `extended_gcd(M.d)[2]`.


# Output

- the Chow ring of the moduli space.
"""
function chow_ring(
  M::QuiverModuliSpace; chi::Union{AbstractVector{Int},UndefInitializer}=undef
)
  if !isdefined(M.chow, :chi)
    if (chi isa UndefInitializer)
      setfield!(M.chow, :chi, extended_gcd(M.d)[2])
    else
      setfield!(M.chow, :chi, chi)
    end
    CH, R, inc = chow_ring(M.Q, M.d, M.theta; chi=M.chow.chi)
    setfield!(M.chow, :ring, CH[1])
    setfield!(M.chow, :_R, R)
    setfield!(M.chow, :_inclusion, inc)
  end

  # M.chow.chi was set
  if !(chi isa UndefInitializer) && M.chow.chi != chi
    # reinitializing all the fields
    setfield!(M.chow, :chi, chi)
    CH, R, inc = chow_ring(M.Q, M.d, M.theta; chi=chi)
    setfield!(M.chow, :ring, CH[1])
    setfield!(M.chow, :_R, R)
    setfield!(M.chow, :_inclusion, inc)
    # linearization changed, so these must be reset
    if isdefined(M.chow, :point)
      setfield!(M.chow, :point, undef)
    end
    if isdefined(M.chow, :todd)
      setfield!(M.chow, :todd, undef)
    end
  end
  return M.chow.ring
end

"""
    extended_gcd(x)

Compute the gcd and the Bezout coefficients of a list of integers.

# Input

- `x`: a list of integers.

# Output

A tuple containing:
- the gcd of the integers,
- a choice of Bezout coefficients.

# Examples

```jldoctest
julia> QuiverTools.extended_gcd([2, 3, 4])
2-element Vector{Any}:
 1
  [-1, 1, 0]

julia> QuiverTools.extended_gcd([2, 3])
2-element Vector{Any}:
 1
  [-1, 1]
```
"""
function extended_gcd(x)
  n = length(x)
  if n == 1
    return [x, [1]]
  elseif n == 2
    g, a, b = gcdx(x[1], x[2])
    return [g, [a, b]]
  else
    g, a, b = gcdx(x[1], x[2])
    y = vcat([g], [x[i] for i in 3:n])
    d, c = extended_gcd(y)
    m = vcat([c[1] * a, c[1] * b], [c[i] for i in 2:(n - 1)])
    return [d, m]
  end
end

"""
    chern_class_line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int})

Compute the first Chern class of the line bundle `L(eta)`.

This is given by ``L(eta) = \\bigoplus_{i \\in Q_0} \\det(U_i)^{-eta_i}``.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `eta::AbstractVector{Int]`: a choice of linearization for the trivial line bundle.

# Output

- the first Chern class of the line bundle L(eta) as a polynomial.

# Examples

The line bundles ``\\mathcal{O}(i)`` on the projective line:
```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> l = chern_class_line_bundle(M, [1, -1])
-x11
```

The line bundle corresponding to the canonical stability condition on our favourite
6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chern_class_line_bundle(M, [9, -6])
-3*x21
```
"""
function chern_class_line_bundle(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int},
)
  A = chow_ring(M)
  I = quotient_ideal(A)
  Rvars = gens(base_ring(I))

  chern_class =
    -sum(eta[i] * Rvars[1 + sum(M.d[1:(i - 1)])] for i in support(M.d))

  return coerce_to_quotient(A, chern_class)
end

"""
    chern_character_line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int})

Compute the Chern character of the line bundle `L(eta)`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `eta::AbstractVector{Int}`: a choice of linearization for the trivial line bundle.

# Output

- the Chern character of the line bundle `L(eta)`.

# Examples

Some line bundles on the projective line:
```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> chern_character_line_bundle(M, [1, -1])
-x11 + 1
```

Some Chern characters for our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chern_character_line_bundle(M, [3, -2])
1//720*x21^6 - 1//120*x21^5 + 1//24*x21^4 - 1//6*x21^3 + 1//2*x21^2 - x21 + 1
```
"""
function chern_character_line_bundle(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int},
)
  x = chern_class_line_bundle(M, eta)
  return sum(x^i / factorial(i) for i in 0:dimension(M))
end

"""
    total_chern_class_universal(M::QuiverModuliSpace, i)

Compute the total Chern class of the universal bundle `U_i`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `i`: the universal bundle we want the Chern class of.

# Output

- the total Chern class of the universal bundle ``U_i(\\chi)``.

# Examples

The universal Chern classes on both vertices of our favourite 3-Kronecker quiver:
```jldoctest
julia> Q  = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> total_chern_class_universal(M, 1)
x11 + x12 + 1

julia> total_chern_class_universal(M, 2)
x21 + x22 + x23 + 1
```
"""
function total_chern_class_universal(
  M::QuiverModuliSpace,
  i::Int,
)
  CH = chow_ring(M)
  CHvars = gens(CH)
  return sum(CHvars[sum(M.d[1:(i - 1)]) + r] for r in 1:M.d[i]; init=CH(0)) + CH(1)
end

"""
    point_class(M::QuiverModuliSpace)

Compute the point class of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the point class of the moduli space, as a polynomial in its Chow ring.

# Examples

A projective 7-fold:
```jldoctest
julia> Q = kronecker_quiver(8);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> chow_ring(M; chi=[1, 0]); point_class(M)
x21^7
```

Our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> point_class(M)
x23^2
```
"""
function point_class(
  M::QuiverModuliSpace
)
  if isdefined(M.chow, :point) && M.chow.point != undef
    return M.chow.point
  end

  CH = chow_ring(M)
  num = 1
  den = 1
  N = dimension(M)

  for i in 1:n_vertices(M.Q)
    c = total_chern_class_universal(M, i)
    num *= c^(M.d' * M.Q.adjacency[:, i])
    den *= c^M.d[i]
  end

  quot = div(num, den)
  pt = sum(
    term for term in Singular.terms(quot) if __chow_ring_monomial_grading(M, term) == N;
    init=CH(0),
  )
  setfield!(M.chow, :point, pt)
  return M.chow.point
end

"""
We call the series ``Q(t) = t/(1-e^{-t})`` the Todd generating series.
The function computes the terms of this series up to degree n.
We use this instead of the more conventional notation `Q` to avoid a
clash with the notation for the quiver.
"""
function todd_Q(t, n)
  return sum((-1)^i * (Nemo.bernoulli(i) * t^i) / factorial(i) for i in 0:n)
end

# TODO rewrite this with in-place operations.
"""
Takes an element in a graded ring and discards all homogeneous components
of degree > n
"""
function truncate(f, n)
  return sum(term for term in Singular.terms(f) if Singular.total_degree(term) <= n)
end

"""
    todd_class(M::QuiverModuliSpace)

Compute the Todd class of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Todd class of the moduli space, as a polynomial in its Chow ring.

# Examples

The Todd class of our favourite 3-Kronecker quiver moduli:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> todd_class(M)
-17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 - 77//60*x22*x23 + x23^2 + 5//12*x12 - 3//2*x21 + 9//8*x23 + 1
```
"""
function todd_class(
  M::QuiverModuliSpace
)
  if isdefined(M.chow, :todd) && M.chow.todd != undef
    return M.chow.todd
  end

  N = dimension(M)
  # consider these constructors: https://nemocas.github.io/AbstractAlgebra.jl/latest/mpolynomial/#Polynomial-functions
  A = chow_ring(M)
  R, inclusion = M.chow._R, M.chow._inclusion
  Rvars = gens(R)

  function xi(i, p)
    return Rvars[sum(M.d[1:(i - 1)]) + p]
  end

  num = 1
  den = 1

  for a in arrows(M.Q)
    i, j = a
    for p in 1:M.d[i]
      for q in 1:M.d[j]
        num *= todd_Q(xi(j, q) - xi(i, p), N)
        num = truncate(num, N)
      end
    end
  end

  for i in 1:n_vertices(M.Q)
    for p in 1:M.d[i]
      for q in 1:M.d[i]
        den *= todd_Q(xi(i, q) - xi(i, p), N)
        den = truncate(den, N)
      end
    end
  end

  # this is because Singular does not have a method to get the preimage
  # of a given element, only ideals.
  # In Singular's implementation this does not result in a loss of time anyways...
  num = gens(preimage(inclusion, Ideal(R, num)))[1]
  den = gens(preimage(inclusion, Ideal(R, den)))[1]

  # renormalizing the constant term because it should be 1, but Singular does not keep
  # it fixed.
  num /= constant_coefficient(num)
  den /= constant_coefficient(den)

  num = coerce_to_quotient(A, num)
  den = coerce_to_quotient(A, den)

  setfield!(M.chow, :todd, A(div(A(num), A(den))))
  return M.chow.todd
end

"""
    integral(M::QuiverModuliSpace, f)

Computes the integral of `f` according to the Hirzebruch-Riemann-Roch theorem.

In other words, it computes the Euler characteristic of the vector bundle
whose Chern character is `f`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `f`: the Chern character in to integrate.

# Output

- the integral of `f`.

# Examples

The integral of ``\\mathcal{O}(i)`` on the projective line for some `i`s.

```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> L = chern_character_line_bundle(M, [1, -1]);

julia> [integral(M, L^i) for i in 0:5]
6-element Vector{Singular.n_Q}:
 1
 2
 3
 4
 5
 6
```

Hilbert series for the 3-Kronecker quiver as in our favourite 6-fold:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> L = chern_character_line_bundle(M, [3, -2]);

julia> [integral(M, L^i) for i in 0:5]
6-element Vector{Singular.n_Q}:
 1
 20
 148
 664
 2206
 5999
```
"""
function integral(M::QuiverModuliSpace, f)
  n = dimension(M)
  integ = div(homogeneous_components(M, f * todd_class(M))[n + 1], point_class(M))
  return Singular.constant_coefficient(integ)
end

integral(F::Bundle) = integral(variety(F), chern_character(F))
chi(F::Bundle) = integral(F::Bundle)

"""
Takes a quotient ring R/I and a polynomial f in R and returns the image of f in R/I.
"""
function coerce_to_quotient(R, f)
  I = quotient_ideal(R)
  q = Singular.reduce(f, I)

  B = base_ring(R)
  g = Singular.MPolyBuildCtx(R)
  for (c, e) in zip(Singular.coefficients(q), Singular.exponent_vectors(q))
    Singular.push_term!(g, B(c), e)
  end
  return Singular.finish(g)
end

"""
Takes a ring R and a polynomial in R/I and returns the canonical preimage of f in R.
"""
function pullback_from_quotient(R, f)
  B = base_ring(R)
  g = Singular.MPolyBuildCtx(R)
  for (c, e) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    Singular.push_term!(g, B(c), e)
  end
  return Singular.finish(g)
end

"""
    __chow_ring__monomial_grading(M::QuiverModuliSpace, f)

Compute the "pseudodegree" of the monomial `f` in the Chow ring of the moduli
space `M` passed.

This method is unsafe, as it does not consider the actual degree of the MPolyRingElem
objects passed. Instead, it assumes that the Chow ring passed has variables
``x_{i, j}`` as in the Chow ring paper.
"""
function __chow_ring_monomial_grading(M::QuiverModuliSpace, f)
  return __chow_degrees(M.d)' * collect(Singular.exponent_vectors(f))[1]
end

"""
    __chow_degrees(d)

Compute the vector of degrees for the variables of a Chow ring.

For internal use only.
"""
@memoize Dict function __chow_degrees(d::AbstractVector{Int})
  return vcat([collect(1:di) for di in d if di > 0]...)
end
