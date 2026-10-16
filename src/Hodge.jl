########################################################################################
# Below lie methods to compute the cohomology of quiver moduli spaces from the motive of
# the stack of semistable representations, computed by the Harder--Narasimhan recursion of
# M. Reineke's paper
# "The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli",
# https://doi.org/10.1007/s00222-002-0273-4
########################################################################################

"""
    hodge_polynomial(M::QuiverModuliSpace; verbose::Bool=false)

Compute the Hodge polynomial of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
function hodge_polynomial(M::QuiverModuliSpace; verbose::Bool=false)
  __check_semistable_moduli_space(M)
  return hodge_polynomial(M.Q, M.d, M.theta, M.denom; verbose=verbose)
end

"""
    hodge_polynomial(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum; verbose::Bool=false)

Return the Hodge polynomial of the moduli space of `theta`-semistable
representations of `Q` with dimension vector `d`.

The cohomology of these moduli spaces is of Hodge--Tate type, so this is the Poincaré
polynomial computed by [`poincare_polynomial`](@ref), with `L` replaced by `x * y`; see
there for the hypotheses.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is `canonical_stability(Q, d)`.
- `denom::Function`: the denominator of the slope. Default is `sum`.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
```
"""
function hodge_polynomial(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum;
  verbose::Bool=false,
)
  P = poincare_polynomial(Q, d, theta, denom; verbose=verbose)
  # the cohomology is of Hodge--Tate type, so L^k contributes x^k * y^k
  S, (x, y) = polynomial_ring(Singular.QQ, ["x", "y"])
  terms = zip(Singular.coefficients(P), Singular.exponent_vectors(P))
  return sum((S(c) * (x * y)^e[1] for (c, e) in terms); init=zero(S))
end

"""
    hodge_diamond(M::QuiverModuliSpace; verbose::Bool=false)

Compute the Hodge diamond of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
function hodge_diamond(M::QuiverModuliSpace; verbose::Bool=false)
  __check_semistable_moduli_space(M)
  return hodge_diamond(M.Q, M.d, M.theta, M.denom; verbose=verbose)
end

"""
    hodge_diamond(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum; verbose::Bool=false)

Compute the Hodge diamond of the moduli space of
`theta`-semistable representations of `Q` with dimension vector `d`.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is the canonical stability.
- `denom::Function`: the denominator of the slope. Default is `sum`.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
  denom::Function=sum;
  verbose::Bool=false,
)
  g = hodge_polynomial(Q, d, theta, denom; verbose=verbose)
  iszero(g) && return zeros(Int, 0, 0)

  # h^{p,p} is the coefficient of x^p * y^p; place the coefficients by exponent, as some
  # of them vanish for cohomology with compact support (quivers with cycles)
  h = zeros(Int, Singular.degree(g, 1) + 1)
  for (c, e) in zip(Singular.coefficients(g), Singular.exponent_vectors(g))
    h[e[1] + 1] = Int(numerator(c))
  end
  return Matrix{Int}(diagonal(h))
end

"""
    picard_rank(M::QuiverModuliSpace; verbose::Bool=false)

Compute the Picard rank of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
function picard_rank(M::QuiverModuliSpace; verbose::Bool=false)
  !(is_smooth(M) && is_projective(M)) &&
    throw(ArgumentError("The moduli space is not smooth and projective"))
  # a point, or the empty moduli space, has trivial Picard group
  dimension(M) <= 0 && return 0
  return betti_numbers(M; verbose=verbose)[3]
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
  has_properly_semistable_representation(M.Q, M.d, M.theta, M.denom) && throw(
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
    betti_numbers(M::QuiverModuliSpace; verbose::Bool=false)

Compute the Betti numbers of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

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
function betti_numbers(M::QuiverModuliSpace; verbose::Bool=false)
  # the empty moduli space has no cohomology
  n = dimension(M)
  n == -Inf && return Int[]

  P = poincare_polynomial(M; verbose=verbose)
  # entry 2k + 1 is the coefficient of L^k in P, i.e., the Betti number
  # in cohomological degree 2k; all odd Betti numbers vanish
  betti = zeros(Int, 2 * n + 1)
  for (c, e) in zip(Singular.coefficients(P), Singular.exponent_vectors(P))
    betti[2 * e[1] + 1] = Int(numerator(c))
  end
  return betti
end

"""
    poincare_polynomial(M::QuiverModuliSpace; verbose::Bool=false)

Compute the Poincaré polynomial of the moduli space `M`.

If the quiver has oriented cycles the moduli space is in general not projective, and this
computes the Poincaré polynomial of its cohomology with compact support instead; with
`verbose = true` this is pointed out by a warning.
As the moduli space is smooth, Poincaré duality relates the two: in dimension
``n`` the compactly supported one is ``L^n P(1/L)`` for the ordinary one ``P(L)``.

If there are properly semistable representations the stack is not a ``\\mathbb{G}_m``-gerbe
over the moduli space, and its motive does not give the cohomology. If the moduli space
is nonetheless smooth, its cohomology is its intersection cohomology, which is then
computed by [`intersection_poincare_polynomial`](@ref); this needs `M.theta` to be generic
for the slope of `M.d` and stable representations to exist. Otherwise an error is thrown.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

# Output

- the Poincaré polynomial of the moduli space, with compact support if the quiver has
  oriented cycles.

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

For dimension vector `[2, 2]` there are properly semistable representations, but the
moduli space is smooth, as it is ``\\mathbb{P}^5``, so its intersection cohomology is used:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> poincare_polynomial(M)
L^5 + L^4 + L^3 + L^2 + L + 1
```

With one vertex, two loops and dimension vector `[1]` the moduli space is the affine
plane, whose only compactly supported class sits in degree 4:
```jldoctest; filter = r"└ @ QuiverTools .*"
julia> M = QuiverModuliSpace(loop_quiver(2), [1], [0]);

julia> poincare_polynomial(M; verbose=true)
┌ Warning: The quiver is not acyclic, so this function actually computes the cohomology with compact support.
└ @ QuiverTools src/Hodge.jl
L^2
```
"""
function poincare_polynomial(M::QuiverModuliSpace; verbose::Bool=false)
  __check_semistable_moduli_space(M)
  return poincare_polynomial(M.Q, M.d, M.theta, M.denom; verbose=verbose)
end

# The methods taking `(Q, d, theta, denom)` describe the semistable moduli space, which is
# the stable one only when there are no properly semistable representations.
function __check_semistable_moduli_space(M::QuiverModuliSpace)
  M.condition == "stable" &&
    !semistable_equals_stable(M) &&
    throw(
      ArgumentError(
        "The stable moduli space differs from the semistable one, and its cohomology is not known."
      ),
    )
  return nothing
end

"""
    poincare_polynomial(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d), denom::Function=sum; verbose::Bool=false)

Compute the Poincaré polynomial of the moduli space of `theta`-semistable representations
of `Q` with dimension vector `d`.

As for the method taking a [`QuiverModuliSpace`](@ref), for a quiver with oriented cycles
this is the Poincaré polynomial of the cohomology with compact support, and with properly
semistable representations it is computed from intersection cohomology if the moduli
space is smooth.

# Input

- `Q::Quiver`: a quiver.
- `d::AbstractVector{Int}`: a dimension vector.
- `theta::AbstractVector{Int}`: a stability parameter. Default is `canonical_stability(Q, d)`.
- `denom::Function`: the denominator of the slope. Default is `sum`.
- `verbose::Bool=false`: whether to warn that the result is for cohomology with compact
  support, which happens when the quiver has oriented cycles.

# Output

- the Poincaré polynomial of the moduli space, with compact support if the quiver has
  oriented cycles.

# Examples

The Poincaré polynomial of our favourite 6-fold:
```jldoctest
julia> poincare_polynomial(kronecker_quiver(3), [2, 3])
L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1
```
"""
function poincare_polynomial(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d),
  denom::Function=sum;
  verbose::Bool=false,
)
  verbose &&
    !is_acyclic(Q) &&
    @warn(
      "The quiver is not acyclic, so this function actually computes the cohomology with compact support."
    )
  if has_properly_semistable_representation(Q, d, theta, denom)
    # the stack is then not a G_m-gerbe over the moduli space, so the motive does not give
    # its cohomology; if the moduli space is smooth, its intersection cohomology does
    M = QuiverModuliSpace(Q, d, theta, "semistable", denom)
    is_smooth(M) || throw(
      ArgumentError(
        "The moduli space has properly semistable representations and is singular, so its cohomology is not known; intersection_poincare_polynomial computes its intersection cohomology."
      ),
    )
    return intersection_poincare_polynomial(M)
  end

  m = motive(Q, d, theta, denom)
  v = Singular.transcendence_basis(Singular.parent(m))[1]
  # the stack of stable representations is a G_m-gerbe over the moduli space, so
  # [M] = (L - 1) * [stack motive]
  P = (v - 1) * m

  denominator(P) != 1 && throw(DomainError("Must be a polynomial"))
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
  # the motive below is that of the semistable stack, which is the stable one exactly when
  # there are no properly semistable representations
  if M.condition == "stable"
    has_properly_semistable_representation(M.Q, M.d, M.theta, M.denom) && throw(
      ArgumentError(
        "The motive of the stable stack is unknown when there are properly semistable " *
        "representations",
      ),
    )
  end
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
  # the computation runs in the fraction field of Q[L], which for answers of high degree
  # is much faster than the function field Q(L); the result is converted at the end
  R, (L,) = polynomial_ring(Singular.QQ, ["L"])
  F = fraction_field(R)
  L = F(L)

  all(ti == 0 for ti in theta) && return __to_function_field(__trivial_motive(Q, d, L))

  ds = all_destabilizing_subdimension_vectors(d, theta, denom)

  push!(ds, zero_vector(Q), d)
  sort!(ds; by=e -> deglex_key(Q, e))

  # the trivial-stability motive of a difference ds[j] - ds[i] occurs for many pairs
  trivial = Dict{Vector{Int},typeof(L)}()
  T = Matrix{typeof(L)}(undef, length(ds), length(ds))
  for (i, j) in Iterators.product(1:length(ds), 1:length(ds))
    if is_subdimension_vector(ds[i], ds[j])
      e = Vector(ds[j] - ds[i])
      T[i, j] =
        power(L, euler_form(Q, ds[i] - ds[j], ds[i])) *
        get!(() -> __trivial_motive(Q, e, L), trivial, e)
    else
      T[i, j] = zero(F)
    end
  end

  # the HN recursion produces the negative of the honest stack motive; negate so
  # this branch agrees with the trivial-stability branch above (see issue #36)
  return __to_function_field(-__first_entry_of_solution(T))
end

# The motive of the stack of all representations of dimension vector `d`, as an element
# of the parent of `L`.
function __trivial_motive(Q::Quiver, d::AbstractVector{Int}, L)
  out = power(L, -euler_form(Q, d, d))
  for i in support(d), nu in 1:d[i]
    out = out // (one(L) - power(L, -nu))
  end
  return out
end

# The first entry of the solution `x` of `T * x = e_n`, for `T` upper triangular. Zero
# entries are skipped, since arithmetic with them is not free in a fraction field.
function __first_entry_of_solution(T::AbstractMatrix)
  n = size(T, 1)
  x = Vector{eltype(T)}(undef, n)
  x[n] = inv(T[n, n])
  for i in (n - 1):-1:1
    s = zero(T[i, i])
    for j in (i + 1):n
      iszero(T[i, j]) || (s += T[i, j] * x[j])
    end
    x[i] = -s // T[i, i]
  end
  return x[1]
end

# Convert a rational function in the fraction field of Q[L] to the function field Q(L),
# in which `motive` returns its result.
function __to_function_field(x)
  K, (L,) = Singular.FunctionField(Singular.QQ, ["L"])
  function to_K(p)
    terms = zip(Singular.coefficients(p), Singular.exponent_vectors(p))
    return sum((K(c) * L^e[1] for (c, e) in terms); init=zero(K))
  end
  return to_K(numerator(x)) // to_K(denominator(x))
end

########################################################################################
# Intersection cohomology of quiver moduli spaces.
#
# When `d` and `theta` admit properly semistable representations the moduli space is
# often singular, and in any case its ordinary cohomology is not determined by the
# Harder--Narasimhan recursion. Meinhardt--Reineke
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

Return the Möbius function ``\\mu(n)``, by trial division.

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
  sign, remaining = 1, n
  for divisor in 2:isqrt(n)
    if iszero(remaining % divisor)
      remaining ÷= divisor
      iszero(remaining % divisor) && return 0
      sign = -sign
    end
  end
  # what is left is 1 or the one prime factor above the square root
  return isone(remaining) ? sign : -sign
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

If no proper subdimension vector has the slope of `M.d` there are no properly semistable
representations, so the moduli space is smooth, its intersection cohomology is ordinary
cohomology, and this agrees with [`poincare_polynomial`](@ref). The converse fails: the
moduli space can be smooth although there are properly semistable representations, see
[`is_smooth`](@ref).

The quiver need not be acyclic, but then the moduli space is affine rather than
projective and, as for [`poincare_polynomial`](@ref), the answer is the Poincaré
polynomial for cohomology with compact support.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Poincaré polynomial of the intersection cohomology of the moduli space.

# Examples

For the 3-Kronecker quiver and dimension vector `[2, 2]` there are properly semistable
representations, as `[1, 1]` has the same slope, so the motive does not give the
cohomology. The moduli space is nonetheless smooth, as it is ``\\mathbb{P}^5``, so this
computes its cohomology, and [`poincare_polynomial`](@ref) uses it:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 2]);

julia> intersection_poincare_polynomial(M)
L^5 + L^4 + L^3 + L^2 + L + 1
```

The dimension vector `[3, 3]` for the same quiver gives a singular 10-fold:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [3, 3]);

julia> intersection_poincare_polynomial(M)
L^10 + L^9 + 2*L^8 + 2*L^7 + 2*L^6 + 2*L^5 + 2*L^4 + 2*L^3 + 2*L^2 + L + 1
```

Reflection functors identify moduli spaces for different dimension vectors:
```jldoctest
julia> Q = kronecker_quiver(4);

julia> M = QuiverModuliSpace(Q, [3, 3]); N = QuiverModuliSpace(Q, [3, 9]);

julia> intersection_poincare_polynomial(M) == intersection_poincare_polynomial(N)
true
```

When there are no properly semistable representations this is the ordinary Poincaré
polynomial:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> intersection_poincare_polynomial(M)
L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1

julia> intersection_poincare_polynomial(M) == poincare_polynomial(M)
true
```

The quiver with one vertex and `m` loops and the trivial stability parameter gives the
classical space of matrix invariants, `m`-tuples of operators on a `d`-dimensional vector
space up to simultaneous conjugation. It is affine of dimension ``(m-1)d^2+1``, singular
except for ``d = 1`` or ``m = d = 2``, and its intersection cohomology is worked out in
Theorem 8.2 of [[MR4000572](https://mathscinet.ams.org/mathscinet-getitem?mr=4000572)]:
```jldoctest
julia> intersection_poincare_polynomial(QuiverModuliSpace(loop_quiver(3), [1]))
L^3

julia> intersection_poincare_polynomial(QuiverModuliSpace(loop_quiver(2), [3]))
L^10

julia> intersection_poincare_polynomial(QuiverModuliSpace(loop_quiver(4), [2]))
L^13 + L^11
```

There is nothing to compute if no representation of dimension vector `M.d` is stable:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(2), [2, 2]);

julia> intersection_poincare_polynomial(M)
ERROR: ArgumentError: There are no stable representations of dimension vector [2, 2]
```

Nor is there anything to compute if the stability parameter is not generic for the slope
of `M.d`, which happens on a fake wall:
```jldoctest
julia> Q = Quiver([0 1 1 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]);

julia> intersection_poincare_polynomial(QuiverModuliSpace(Q, [3, 3, 4, 1]))
ERROR: ArgumentError: The stability parameter is not generic for the slope of [3, 3, 4, 1], so intersection cohomology is out of reach
```
"""
function intersection_poincare_polynomial(M::QuiverModuliSpace)
  M.condition == "semistable" || throw(
    ArgumentError("Intersection cohomology is computed for the semistable moduli space")
  )
  Q, d, theta, denom = M.Q, M.d, M.theta, M.denom
  mu = slope(d, theta, denom)

  # the dimension vectors of the slope of `d`, i.e. the part of the monoid
  # ``\Lambda_\mu`` below `d`; the zero vector is left out and treated separately
  lattice = filter(
    e -> slope(e, theta, denom) == mu, all_subdimension_vectors(d; nonzero=true)
  )
  all(euler_form(Q, e, f) == euler_form(Q, f, e) for e in lattice, f in lattice) || throw(
    ArgumentError(
      "The stability parameter is not generic for the slope of $(Vector(d)), " *
      "so intersection cohomology is out of reach",
    ),
  )
  has_stable_representation(Q, d, theta, denom) || throw(
    ArgumentError("There are no stable representations of dimension vector $(Vector(d))")
  )

  R, ws = polynomial_ring(Singular.QQ, ["w"])
  w = ws[1]
  F = fraction_field(R)
  E = Singular.elem_type(F)
  root = F(-w)                     # ``L^{1/2}``, which sits in odd degree
  # substitute `image` for the variable of a univariate rational function
  substitute(num, den, image) = F(num(image))//F(den(image))
  adams(x, n) = substitute(numerator(x), denominator(x), w^n)
  # the motive lives in ``\mathbb{Q}(L)``; move it to the square root, ``L = w^2``
  to_square_root(m) = substitute(
    Singular.n_transExt_to_spoly(numerator(m)),
    Singular.n_transExt_to_spoly(denominator(m)),
    w^2,
  )

  # the generating series, without its constant term
  series = Dict{Vector{Int},E}(
    e => power(root, euler_form(Q, e, e)) * to_square_root(motive(Q, e, theta, denom))
    for e in lattice
  )

  # the ordinary logarithm ``\log(1 + x) = \sum_k (-1)^{k-1}/k x^k``; the kth power is
  # supported on sums of k nonzero dimension vectors, so the sum stops at ``|d|``
  logarithm = Dict{Vector{Int},E}()
  term = series
  for k in 1:sum(d)
    isempty(term) && break
    scale = F((-1)^(k - 1))//k
    for (e, value) in term
      logarithm[e] = get(logarithm, e, zero(F)) + scale * value
    end
    k == sum(d) && break
    # multiply by `series`, dropping everything that is no longer below `d`
    next = Dict{Vector{Int},E}()
    for (e, left) in term, (f, right) in series
      is_subdimension_vector(e + f, d) || continue
      next[e + f] = get(next, e + f, zero(F)) + left * right
    end
    term = next
  end

  # and the plethystic one, ``Log(1 + x) = \sum_n \mu(n)/n \psi^n(\log(1 + x))``; only
  # those `n` with `n * e = d` for some `e` contribute, i.e. the divisors of `gcd(d)`
  target = Vector{Int}(d)
  total = get(logarithm, target, zero(F))
  common = gcd(d)
  for k in 2:common
    iszero(common % k) || continue
    mobius = _mobius(k)
    iszero(mobius) && continue
    total += F(mobius)//k * adams(get(logarithm, target .÷ k, zero(F)), k)
  end

  # ``DT_d = (L^{1/2} - L^{-1/2}) [Log Q]_{t^d}``, then ``E(IH^*) = L^{\dim/2} DT_d``
  result = power(root, 1 - euler_form(Q, d, d)) * (root - inv(root)) * total
  # by [MR4000572] this is the intersection cohomology of a nonempty variety, so it
  # cannot vanish once stable representations exist
  iszero(result) && error(
    "The Donaldson--Thomas invariant of $(Vector(d)) vanishes although stable " *
    "representations exist; this contradicts [MR4000572] and is a bug, please report it",
  )
  isone(denominator(result)) ||
    throw(DomainError("Intersection cohomology is not polynomial"))

  # intersection cohomology of these moduli spaces is concentrated in even degree, so the
  # answer has to be a polynomial in `w^2`; that is a real check on the whole computation
  S, Ls = polynomial_ring(Singular.QQ, ["L"])
  L = Ls[1]
  invariant = numerator(result)
  P = zero(S)
  for (c, e) in
      zip(Singular.coefficients(invariant), Singular.exponent_vectors(invariant))
    isodd(e[1]) && throw(DomainError("Intersection cohomology in odd degree"))
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
  degree ``0, \\dots, 2\\dim M``. The odd ones vanish.

# Examples

For the 3-Kronecker quiver and dimension vector `[2, 2]` there are properly semistable
representations, but the moduli space is smooth: it is ``\\mathbb{P}^5``,
so we get the expected Betti numbers:
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

Intersection cohomology of a projective variety still satisfies Poincaré duality, so the
Betti numbers of a singular moduli space are palindromic just as well:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [2, 4]);

julia> betti = intersection_betti_numbers(M);

julia> betti == reverse(betti)
true

julia> betti[1:2:end]
14-element Vector{Int64}:
 1
 1
 3
 4
 6
 6
 7
 7
 6
 6
 4
 3
 1
 1
```

When there are no properly semistable representations these are the ordinary Betti
numbers:
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

The singular moduli space for the 4-Kronecker quiver and dimension vector `[2, 2]`.
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [2, 2]);

julia> intersection_hodge_diamond(M)
10×10 Matrix{Int64}:
 1  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0  0
 0  0  2  0  0  0  0  0  0  0
 0  0  0  2  0  0  0  0  0  0
 0  0  0  0  2  0  0  0  0  0
 0  0  0  0  0  2  0  0  0  0
 0  0  0  0  0  0  2  0  0  0
 0  0  0  0  0  0  0  2  0  0
 0  0  0  0  0  0  0  0  1  0
 0  0  0  0  0  0  0  0  0  1
```

A case where the canonical stability parameter is on a wall; the moduli space is a 5-fold
with finitely many singular points.
```jldoctest
julia> Q = Quiver("1--2,2--3"); M = QuiverModuliSpace(Q, [2, 2, 2]);

julia> is_smooth(M)
false

julia> codimension_singular_locus(M)
5

julia> intersection_hodge_diamond(M)
6×6 Matrix{Int64}:
 1  0  0  0  0  0
 0  2  0  0  0  0
 0  0  3  0  0  0
 0  0  0  3  0  0
 0  0  0  0  2  0
 0  0  0  0  0  1
```

The moduli space can be smooth even when there are properly semistable representations.
Its intersection cohomology is then its ordinary cohomology, which is how
[`hodge_diamond`](@ref) computes it. For the 3-Kronecker quiver and dimension vector
`[2, 2]` the moduli space is ``\\mathbb{P}^5``:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 2]);

julia> is_smooth(M)
true

julia> intersection_hodge_diamond(M)
6×6 Matrix{Int64}:
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  1  0  0  0
 0  0  0  1  0  0
 0  0  0  0  1  0
 0  0  0  0  0  1

julia> hodge_diamond(M) == intersection_hodge_diamond(M)
true
```

Likewise here the canonical stability parameter is on a wall, yet the moduli space is a
smooth 5-fold:
```jldoctest
julia> Q = Quiver("1-2,1-3,2--3"); d = [2, 2, 2]; M = QuiverModuliSpace(Q, d);

julia> QuiverTools.__is_vgit_chamber_parameter(Q, d, M.theta)
false

julia> is_smooth(M)
true

julia> intersection_hodge_diamond(M)
6×6 Matrix{Int64}:
 1  0  0  0  0  0
 0  2  0  0  0  0
 0  0  3  0  0  0
 0  0  0  3  0  0
 0  0  0  0  2  0
 0  0  0  0  0  1
```

When there are no properly semistable representations this is the ordinary Hodge diamond:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [2, 3]);

julia> intersection_hodge_diamond(M) == hodge_diamond(M)
true
```
"""
function intersection_hodge_diamond(M::QuiverModuliSpace)
  return Matrix{Int}(diagonal(intersection_betti_numbers(M)[1:2:end]))
end
