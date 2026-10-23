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
    __simplify(f)

Force some Gröbner basis simplification of the polynomial `f`.

For internal use only.
"""
__simplify(f) = Singular.reduce(f, Singular.Ideal(parent(f), zero(f)))

"""
    __truncate(M::QuiverModuliSpace, x, n)

Truncate a polynomial in the Chow ring,
discarding terms of degree higher than n.

For internal use only.
"""
function __truncate(M::QuiverModuliSpace, x, n)
  degrees = __chow_degrees(M.d)
  context = Singular.MPolyBuildCtx(parent(x))
  for (c, e) in zip(Singular.coefficients(x), Singular.exponent_vectors(x))
    e' * degrees <= n && Singular.push_term!(context, c, e)
  end
  return Singular.finish(context)
end

"""
    __homogeneous_components(M::QuiverModuliSpace, x; unsafe::Bool=false)
    __homogeneous_components(M::QuiverModuliSpace, x, k::Int)

# Input

- `M::QuiverModuliSpace`: a quiver moduli space.
- `x`: an element of the Chow ring of `M`.
- `k::Int`: a degree.
- `unsafe::Bool=false`: whether to compute the dimension of `M` using
  the faster Euler form instead of `dimension()`. Default is `false`.

Decompose a Chow ring element `x` into its homogeneous components.
If the degree `k` is given, return only the homogeneous component of degree `k`.

For internal use only.
"""
function __homogeneous_components(M::QuiverModuliSpace, x; unsafe::Bool=false)
  unsafe ? (n = 1 - euler_form(M.Q, M.d, M.d)) : (n = dimension(M))

  degrees = __chow_degrees(M.d)
  contexts = [Singular.MPolyBuildCtx(parent(x)) for _ in 0:n]
  for (c, e) in zip(Singular.coefficients(x), Singular.exponent_vectors(x))
    k = e' * degrees
    k <= n && Singular.push_term!(contexts[k + 1], c, e)
  end
  return map(Singular.finish, contexts)
end

function __homogeneous_components(M::QuiverModuliSpace, x, k::Int)
  degrees = __chow_degrees(M.d)
  context = Singular.MPolyBuildCtx(parent(x))
  for (c, e) in zip(Singular.coefficients(x), Singular.exponent_vectors(x))
    e' * degrees == k && Singular.push_term!(context, c, e)
  end
  return Singular.finish(context)
end

"""
    symmetric_polynomial(degree::Int)

Return the symmetric polynomial of degree `degree` in the variables `vars`
as a Julia function.

# Input

- `vars`: a list of variables.
- `degree`: the degree of the wanted symmetric polynomial.

# Output

- The symmetric polynomial of degree `degree` in the variables `vars`.

# Examples

```julia-repl
julia> using Singular;

julia> R, vars = polynomial_ring(Singular.QQ, ["x", "y", "z"]);

julia> f = QuiverTools.symmetric_polynomial(2); f(vars)
x*y + x*z + y*z
```
"""
function symmetric_polynomial(degree::Int)
  f(vars) = sum(prod(e; init=1) for e in IterTools.subsets(vars, degree))
  return f
end
"""
    product_lists(L)

For internal use only.
"""
function product_lists(L)
  length(L) == 1 && return L[1]
  return [p * l for p in product_lists(L[1:(end - 1)]) for l in L[end]]
end

"""
    chow_ring(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d); chi::AbstractVector{Int}=extended_gcd(d)[2], verbose::Bool=false, unsafe::Bool=false)

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
julia> CH = chow_ring(kronecker_quiver(2), [1, 1])[1][1];

julia> QuiverTools.gens(QuiverTools.quotient_ideal(CH))
2-element Vector{Singular.spoly{Singular.n_Q}}:
 x2_1
 x1_1^2
```

The Chow ring for our favourite 6-fold has, in this implementation, 16 generators:
```jldoctest
julia> CH = chow_ring(kronecker_quiver(3), [2, 3])[1][1]; I = QuiverTools.quotient_ideal(CH);

julia> length(QuiverTools.gens(I))
16
```
"""
function chow_ring(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d);
  chi::AbstractVector{Int}=extended_gcd(d)[2],
  verbose::Bool=false,
  unsafe::Bool=false,
)
  chi' * d != 1 && throw(ArgumentError("``chi`` is not a linearization"))
  if !unsafe
    has_properly_semistable_representation(Q, d, theta) &&
      throw(
        ArgumentError(
          "The quiver moduli problem has properly semistable representations, no description of the Chow ring is available."
        ),
      )
    !is_amply_stable(Q, d, theta) && throw(
      ArgumentError(
        "The quiver moduli problem is not amply stable, no description of the Chow ring is available."
      ),
    )
  else
    verbose && @warn "Unsafe computation."
  end

  # j varies first, then i
  varnames = ["xi$(i)_$(j)" for i in 1:n_vertices(Q) for j in 1:d[i]]
  R, vars = Singular.polynomial_ring(Singular.QQ, varnames)

  # the variable `xi_{i,j}` has index var_indices[i] + j
  var_indices = __var_indices(d)

  # Shorthand to address the variable `xi_{i,j}`.
  function xi(i, j)
    d[i] == 0 && throw(ArgumentError("i is not in the support of d."))
    return vars[var_indices[i] + j]
  end

  # build a base of R as an A-module.
  bounds = UnitRange{Int64}[0:(d[i] - nu) for i in 1:n_vertices(Q) for nu in 1:d[i]]
  # the exponent of xi(i, nu) is lambda[var_indices[i] + nu], so `lambda` is the exponent
  # vector of the monomial
  function build_elem(lambda::NTuple)
    context = Singular.MPolyBuildCtx(R)
    Singular.push_term!(context, Singular.QQ(1), collect(lambda))
    return Singular.finish(context)
  end
  base = map(build_elem, Iterators.product(bounds...))

  verbose && @info "base has $(length(base)) elements"

  # the variables xi(i, 1), ..., xi(i, d[i]) of each vertex i in the support of d
  blocks = [(var_indices[i] + 1):(var_indices[i] + d[i]) for i in support(d)]

  # the group W, product of the symmetric groups on the blocks. Each element is stored as
  # its sign and the permutation it induces on the indices of the variables.
  W = [
    (
      prod(sign(Singular.perm(wi)) for wi in w; init=1),
      Int[var_indices[i] + w[i][j] for i in support(d) for j in 1:d[i]],
    ) for
    w in Iterators.product([Combinatorics.permutations(1:d[i]) for i in 1:n_vertices(Q)]...)
  ]

  # Sort the exponents of each block decreasingly, in place, and return the sign of the
  # sorting permutation. Return 0 if a block has a repeated exponent.
  function sort_blocks!(e::Vector{Int})
    s = 1
    for b in blocks
      for k in b, l in (k + 1):last(b)
        e[k] == e[l] && return 0
        e[k] < e[l] && (s = -s)
      end
      sort!(view(e, b); rev=true)
    end
    return s
  end

  # The discriminant in the definition of the antisymmetrization.
  delta = prod(
    xi(i, l) - xi(i, k) for i in 1:n_vertices(Q) for k in 1:(d[i] - 1) for l in (k + 1):d[i];
    init=one(R),
  )

  # The antisymmetrization A(f) of f, divided by delta.
  # If e = b ∘ tau then A(x^e) = sign(tau) A(x^b), and A(x^e) = 0 if a block of e has a
  # repeated exponent. So it suffices to collect the coefficients of f on each W-orbit of
  # exponents, represented by its decreasing rearrangement b, and to expand A(x^b) only for
  # the orbits whose coefficient is nonzero.
  function antisymmetrize(f::Singular.spoly{Singular.n_Q})
    orbits = Dict{Vector{Int},Singular.n_Q}()
    for (c, e) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
      s = sort_blocks!(e)
      s == 0 && continue
      orbits[e] = get(orbits, e, zero(c)) + (s == 1 ? c : -c)
    end
    # distinct orbits give distinct monomials, so the context has no terms to merge
    context = Singular.MPolyBuildCtx(R)
    exponents = Vector{Int}(undef, length(vars))
    for (b, c) in orbits
      iszero(c) && continue
      for (s, indices) in W
        for k in eachindex(indices)
          exponents[k] = b[indices[k]]
        end
        Singular.push_term!(context, s == 1 ? c : -c, exponents)
      end
    end
    # the division is exact, which is faster than just `div`
    return Singular.divexact(Singular.finish(context), delta)
  end

  # All the destabilizing subdimension vectors of `d` with respect to the slope
  # `theta/denom` that are minimal with respect to the total order.
  minimal_forbidden = all_destabilizing_subdimension_vectors(d, theta)
  filter!(
    e -> !any(partial_order(Q, f, e) for f in minimal_forbidden if f != e),
    minimal_forbidden,
  )

  verbose &&
    @info "there are $(length(minimal_forbidden)) minimal forbidden dimension vectors"

  # builds a new forbidden polynomial for the minimal forbidden dimension vector e.
  function new_forbidden(e::AbstractVector{Int})
    out = one(R)
    for (i, j) in Iterators.product(1:n_vertices(Q), 1:n_vertices(Q))
      for r in 1:e[i], s in (e[j] + 1):d[j]
        Q.adjacency[i, j] == 0 && continue
        Singular.mul!(out, out, (xi(j, s) - xi(i, r))^Q.adjacency[i, j])
      end
    end
    return out
  end
  forbidden_polynomials = Singular.spoly{Singular.n_Q}[
    new_forbidden(e) for e in minimal_forbidden
  ]

  varnames2 = ["x$(i)_$(j)" for i in 1:n_vertices(Q) for j in 1:d[i]]
  A, Avars = polynomial_ring(Singular.QQ, varnames2)

  # Shorthand to address the variables of A `x_{i,j}`.
  function xs(i, j)
    d[i] == 0 && throw(ArgumentError("i is not in the support of d."))
    return Avars[var_indices[i] + j]
  end

  symm_polys = [symmetric_polynomial(k) for k in 1:maximum(d)]
  targets = Singular.spoly{Singular.n_Q}[]
  for i in support(d)
    verbose && @info "computing the targets for vertex $(i) out of $(length(support(d)))"
    for k in 1:d[i]
      push!(targets, symm_polys[k]([xi(i, j) for j in 1:d[i]]))
    end
  end
  verbose && @info "there are $(length(targets)) targets"

  inclusion = Singular.AlgebraHomomorphism(A, R, targets)
  verbose && @info "the inclusion map is built"

  anti = Singular.spoly{Singular.n_Q}[]
  verbose && @info "antisymmetrizing the forbidden polynomials, this may take a while..."

  for i in eachindex(forbidden_polynomials)
    verbose && @info "forbidden polynomial $(i) out of $(length(forbidden_polynomials))"
    for b in base
      a = antisymmetrize(forbidden_polynomials[i] * b)
      a != 0 && push!(anti, a)
    end
  end

  verbose && @info "there are $(length(anti)) antisymmetrized forbidden polynomials"

  tautological = Singular.spoly{Singular.n_Q}[
    gens(preimage(inclusion, Singular.Ideal(R, g)))[1] for g in anti
  ]
  verbose && @info "there are $(length(tautological)) tautological polynomials"

  linear = Singular.spoly{Singular.n_Q}[sum(chi[i] * xs(i, 1) for i in support(d))]

  return (
    Singular.QuotientRing(A, std(Singular.Ideal(A, [tautological; linear]))), R, inclusion
  )
end

"""
    chow_ring(M::QuiverModuliSpace; chi::Union{AbstractVector{Int},UndefInitializer}=undef, verbose::Bool=false, unsafe::Bool=false)

Compute the Chow ring of the moduli space `M` for the given linearization `chi`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `chi::AbstractVector{Int}`: a choice of linearization for the trivial line bundle.
  Default is `extended_gcd(M.d)[2]`.


# Output

- the Chow ring of the moduli space.
"""
function chow_ring(
  M::QuiverModuliSpace; chi::Union{AbstractVector{Int},UndefInitializer}=undef,
  verbose::Bool=false,
  unsafe::Bool=false,
)
  # `unsafe` asserts the stable locus is nonempty, so seed the dimension cache with the
  # cheap `1 - <d, d>` before any bundle is built. This stops downstream `dimension(M)`
  # calls from triggering the expensive `has_stable_representation` check (issue #20).
  if unsafe && M.chow._dimension === nothing
    setfield!(M.chow, :_dimension, 1 - euler_form(M.Q, M.d, M.d))
  end

  # gate on `:ring`, not `:chi`: `set_linearization!` sets `chi` without building the
  # ring, so keying off `:chi` would skip the build and leave `ring` undefined (issue #16)
  if !isdefined(M.chow, :ring)
    if !isdefined(M.chow, :chi)
      setfield!(M.chow, :chi, chi isa UndefInitializer ? extended_gcd(M.d)[2] : chi)
    elseif !(chi isa UndefInitializer)
      setfield!(M.chow, :chi, chi)
    end
    CH, R, inc = chow_ring(
      M.Q, M.d, M.theta; chi=M.chow.chi, verbose=verbose, unsafe=unsafe
    )
    setfield!(M.chow, :ring, CH[1])
    setfield!(M.chow, :_R, R)
    setfield!(M.chow, :_inclusion, inc)
  end

  # M.chow.chi was set
  if !(chi isa UndefInitializer) && M.chow.chi != chi
    # reinitializing all the fields
    setfield!(M.chow, :chi, chi)
    CH, R, inc = chow_ring(M.Q, M.d, M.theta; chi=chi, verbose=verbose, unsafe=unsafe)
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
    chern_class_line_bundle(M::QuiverModuliSpace, eta::AbstractVector{Int}; unsafe::Bool=false)

Compute the first Chern class of the line bundle `L(eta)`.

This is given by ``L(eta) = \\bigoplus_{i \\in Q_0} \\det(U_i)^{-eta_i}``.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `eta::AbstractVector{Int]`: a choice of linearization for the trivial line bundle.
- `unsafe::Bool=false`: whether to skip the checks on ample stability and
    existence of properly semistables. Default is `false`.

# Output

- the first Chern class of the line bundle L(eta) as a polynomial.

# Examples

The line bundles ``\\mathcal{O}(i)`` on the projective line:
```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> l = chern_class_line_bundle(M, [1, -1])
-x1_1
```

The line bundle corresponding to the canonical stability condition on our favourite
6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chern_class_line_bundle(M, [9, -6])
-3*x2_1
```
"""
function chern_class_line_bundle(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int};
  unsafe::Bool=false,
)
  A = chow_ring(M; unsafe=unsafe)
  R = base_ring(quotient_ideal(A))
  Rvars = gens(R)
  proj = __projection_to_quotient_ring(A)

  var_indices = __var_indices(M.d)
  chern_class = zero(R)
  for i in support(M.d)
    Singular.add!(chern_class, chern_class, eta[i] * Rvars[var_indices[i] + 1])
  end

  return __simplify(proj(- chern_class))
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
-x1_1 + 1
```

Some Chern characters for our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> chern_character_line_bundle(M, [3, -2])
1//720*x2_1^6 - 1//120*x2_1^5 + 1//24*x2_1^4 - 1//6*x2_1^3 + 1//2*x2_1^2 - x2_1 + 1
```
"""
function chern_character_line_bundle(
  M::QuiverModuliSpace,
  eta::AbstractVector{Int},
)
  x = chern_class_line_bundle(M, eta)
  CH = parent(x)
  ch = one(CH)
  term = one(CH)
  for i in 1:dimension(M)
    term *= x / big(i) # x^i / i!
    ch += term
  end
  return ch
end

"""
    chern_character_tangent_bundle(M::QuiverModuliSpace; unsafe::Bool=false)

Compute the Chern character of the tangent bundle of the quiver moduli space `M`.

Uses the identity
``[\\mathrm{T}_M] = \\sum_{a: i \\to j} [U_i^\\vee \\otimes U_j] - \\sum_i [U_i^\\vee \\otimes U_i] + [\\mathcal{O}_M]``
in the Grothendieck group of `M`.

# Example

The Chern character of the tangent bundle of `P^3`, realized as the moduli space of the
generalized Kronecker quiver with 4 arrows and dimension vector `(1, 1)`:

```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [1, 1]); chow_ring(M);

julia> chern_character_tangent_bundle(M)
-2//3*x1_1^3 + 2*x1_1^2 - 4*x1_1 + 3
```
"""
function chern_character_tangent_bundle(M::QuiverModuliSpace; unsafe::Bool=false)
  CH = chow_ring(M; unsafe=unsafe)
  # the Chern characters of the universal bundles, from their Chern classes
  classes = [
    __homogeneous_components(
      M, total_chern_class_universal(M, i; unsafe=unsafe); unsafe=unsafe
    ) for i in 1:n_vertices(M.Q)
  ]
  ch = [
    __chern_character_from_classes(M, Dict(k - 1 => c[k] for k in eachindex(c)), M.d[i])
    for (i, c) in enumerate(classes)
  ]
  ch_dual = [dual(M, x) for x in ch]
  chT = one(CH)
  chT += sum(ch_dual[i] * ch[j] for (i, j) in arrows(M.Q); init=zero(CH))
  chT -= sum(ch_dual[i] * ch[i] for i in 1:n_vertices(M.Q); init=zero(CH))
  return __simplify(chT)
end

"""
    total_chern_class_universal(M::QuiverModuliSpace, i::Int; unsafe::Bool=false)

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
x1_1 + x1_2 + 1

julia> total_chern_class_universal(M, 2)
x2_1 + x2_2 + x2_3 + 1
```
"""
function total_chern_class_universal(
  M::QuiverModuliSpace,
  i::Int;
  unsafe::Bool=false,
)
  CH = chow_ring(M; unsafe=unsafe)
  CHvars = gens(CH)
  var_indices = __var_indices(M.d)
  return sum(CHvars[var_indices[i] + r] for r in 1:M.d[i]; init=zero(CH)) + one(CH)
end

"""
    segre_class(M::QuiverModuliSpace, c; unsafe::Bool=false)

Compute the total Segre class with total Chern class `c`, i.e. the inverse of `c` in the
Chow ring of `M`.

Writing ``c_j`` for the part of degree ``j`` of `c`, which must have ``c_0 = 1``, the part
``s_k`` of degree ``k`` of the Segre class is given by ``s_0 = 1`` and
```math
s_k = -\\sum_{j=1}^{k} c_j s_{k-j}
```
for ``k \\geq 1``, which is the part of degree ``k`` of ``c \\cdot s = 1``.

# Examples

The Segre class of the universal bundle ``U_1`` on the projective plane:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(3), [1, 1]);

julia> c = total_chern_class_universal(M, 1)
x1_1 + 1

julia> segre_class(M, c)
x1_1^2 - x1_1 + 1
```
"""
function segre_class(M::QuiverModuliSpace, c; unsafe::Bool=false)
  unsafe ? (n = 1 - euler_form(M.Q, M.d, M.d)) : (n = dimension(M))
  cs = __homogeneous_components(M, c; unsafe=unsafe)
  isone(cs[1]) ||
    throw(ArgumentError("the total Chern class must have constant term 1"))
  CH = parent(c)
  s = [one(CH)]
  for k in 1:n
    s_k = -sum(cs[j + 1] * s[k - j + 1] for j in 1:k if !iszero(cs[j + 1]); init=zero(CH))
    push!(s, __simplify(s_k))
  end
  return sum(s)
end

"""
    adams(M::QuiverModuliSpace, x, k)

Compute the Adams operation ``\\Phi^k`` on the Chern character `x`, i.e. multiply its part
of degree ``i`` by ``k^i``.
For internal use only.
"""
function adams(M::QuiverModuliSpace, x, k)
  n = dimension(M)
  # k^n overflows an Int already for k = n = 16
  return [big(k)^i for i in 0:n]' * __homogeneous_components(M, x)
end

"""
    dual(M::QuiverModuliSpace, x)

Compute the Chern character of the dual of a bundle with Chern character `x`, i.e. the
Adams operation ``\\Phi^{-1}`` on `x`.

# Example

On the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); d = [1, 1]; a = [1, 0];

julia> M = QuiverModuliSpace(Q, d); chow_ring(M; chi=a);

julia> dual(M, todd_class(M))
-x2_1 + 1
```
"""
dual(M::QuiverModuliSpace, x) = adams(M, x, -1)

# The Chern classes of a bundle with Chern character `x`, as a dictionary from the degree
# to the class, computed with the Newton identities.
function __chern_classes_from_character(M::QuiverModuliSpace, x)
  CH = parent(x)
  n = dimension(M)
  comps = __homogeneous_components(M, x)
  p = [(CH(-1))^i * CH(factorial(big(i))) * comps[i + 1] for i in 0:n]
  e = [zero(CH) for _ in 1:(n + 1)]
  e[1] = one(CH)
  for i in 1:n
    e[i + 1] = CH(-1//i) * sum(p[j + 1] * e[i - j + 1] for j in 1:i)
    e[i + 1] = __simplify(e[i + 1])
  end
  return Dict(i => e[i + 1] for i in 0:n)
end

# The Chern character of a bundle of rank `r` with Chern classes `e`, a dictionary from the
# degree to the class, computed with the Newton identities.
function __chern_character_from_classes(M::QuiverModuliSpace, e, r)
  CH = parent(e[0])
  n = dimension(M)
  n == 0 && return zero(CH)
  p = vcat([-e[1]], [zero(CH) for _ in 1:(n - 1)])
  for i in 1:(n - 1)
    p[i + 1] = -CH(i + 1) * e[i + 1] - sum(e[j] * p[i - j + 1] for j in 1:i)
  end
  return __simplify(sum(CH((-1)^i//factorial(big(i))) * p[i] for i in 1:n) + r)
end

"""
    __chern_characters_wedge(M::QuiverModuliSpace, x, k)

Compute the Chern characters of the exterior powers, up to degree `k`, of a bundle with
Chern character `x`.
For internal use only.
"""
function __chern_characters_wedge(M::QuiverModuliSpace, x, k)
  CH = parent(x)
  k == 0 && return [one(CH)]
  n = dimension(M)

  # init as zero(CH) for type stability
  wedges = [zero(CH) for _ in 1:(k + 1)]
  wedges[1], wedges[2] = one(CH), x

  for j in 2:k
    wedges[j + 1] =
      CH(1//j) * __truncate(M,
        sum(
          (-one(CH))^(j - i + 1) * wedges[i + 1] * adams(M, x, j - i) for i in 0:(j - 1);
          init=zero(CH),
        ),
        n)
    wedges[j + 1] = __simplify(wedges[j + 1])
  end
  return wedges
end

"""
    __chern_characters_symmetric(M::QuiverModuliSpace, x, r, k)

Compute the Chern characters of the symmetric powers, up to degree `k`, of a bundle of
rank `r` with Chern character `x`.
For internal use only.
"""
function __chern_characters_symmetric(M::QuiverModuliSpace, x, r, k)
  CH = parent(x)
  k == 0 && return [one(CH)]
  n = dimension(M)

  wedges = __chern_characters_wedge(M, x, r)
  # init as zero(CH) for type stability
  syms = [zero(CH) for _ in 1:(k + 1)]
  syms[1], syms[2] = one(CH), x

  for j in 2:k
    syms[j + 1] = __truncate(M,
      sum(
        (-one(CH))^(i + 1) * wedges[i + 1] * syms[j - i + 1] for i in 1:min(j, r);
        init=zero(CH),
      ),
      n)
    syms[j + 1] = __simplify(syms[j + 1])
  end
  return syms
end

"""
    point_class(M::QuiverModuliSpace; unsafe::Bool=false)

Compute the point class of the moduli space `M`.

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.
- `unsafe::Bool=false`: whether to skip the checks on ample stability and
    existence of properly semistables. Default is `false`.

# Output

- the point class of the moduli space, as a polynomial in its Chow ring.

# Examples

A projective 7-fold:
```jldoctest
julia> Q = kronecker_quiver(8);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> chow_ring(M; chi=[1, 0]); point_class(M)
x2_1^7
```

Our favourite 6-fold:
```jldoctest
julia> Q = kronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> point_class(M)
x2_3^2
```

The 7-subspace quiver:
```jldoctest
julia> Q = subspace_quiver(7); d = push!(ones(Int, 7), 2); M = QuiverModuliSpace(Q, d);

julia> point_class(M; unsafe=true)
1//10*x8_1^4
```
"""
function point_class(
  M::QuiverModuliSpace;
  unsafe::Bool=false,
)
  if isdefined(M.chow, :point) && M.chow.point != undef
    return M.chow.point
  end

  CH = chow_ring(M; unsafe=unsafe)
  unsafe ? (N = 1 - euler_form(M.Q, M.d, M.d)) : (N = dimension(M))

  # the point class is the degree N component of the product of c(U_i)^m_i.
  # A negative power of c(U_i) is a power of its Segre class, so we avoid slow divisions.
  m = M.Q.adjacency' * M.d - M.d
  num = one(CH)
  for i in 1:n_vertices(M.Q)
    (m[i] == 0 || M.d[i] == 0) && continue
    c = total_chern_class_universal(M, i; unsafe=unsafe)
    factor = m[i] > 0 ? c : segre_class(M, c; unsafe=unsafe)
    for _ in 1:abs(m[i])
      num = __simplify(__truncate(M, num * factor, N))
    end
  end

  setfield!(M.chow, :point, __homogeneous_components(M, num, N))
  return M.chow.point
end

"""
    todd_class(M::QuiverModuliSpace)

Compute the Todd class of the moduli space `M`.

The Todd class is computed from the Chern character of the tangent bundle, as
```math
\\operatorname{td}(T_M) = \\exp\\Big(\\frac{1}{2}\\operatorname{ch}_1(T_M)
- \\sum_{m \\geq 1} \\frac{B_{2m}}{2m} \\operatorname{ch}_{2m}(T_M)\\Big).
```

# Input

- `M::QuiverModuliSpace`: a moduli space of representations of a quiver.

# Output

- the Todd class of the moduli space, as a polynomial in its Chow ring.

# Examples

The Todd class of our favourite 3-Kronecker quiver moduli:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> todd_class(M)
-17//8*x1_2*x2_1 + x2_1^2 + 823//360*x1_2*x2_2 - 823//1080*x2_2^2 + 553//1080*x2_1*x2_3 - 77//60*x2_2*x2_3 + x2_3^2 + 5//12*x1_2 - 3//2*x2_1 + 9//8*x2_3 + 1
```
"""
function todd_class(
  M::QuiverModuliSpace;
  unsafe::Bool=false,
)
  if isdefined(M.chow, :todd) && M.chow.todd != undef
    return M.chow.todd
  end

  CH = chow_ring(M; unsafe=unsafe)
  N = dimension(M)

  # The Todd class is multiplicative: td(E) = prod_j Q(t_j) over the Chern roots t_j of E,
  # with Q(t) = t/(1 - e^{-t}). Since log Q(t) = t/2 - sum_{m >= 1} B_{2m}/(2m (2m)!) t^{2m}
  # and sum_j t_j^k = k! ch_k(E), this gives
  # log td(E) = ch_1(E)/2 - sum_{m >= 1} B_{2m}/(2m) ch_{2m}(E).
  ch = __homogeneous_components(
    M, chern_character_tangent_bundle(M; unsafe=unsafe); unsafe=unsafe
  ) # ch[k + 1] = ch_k(T_M)
  y = N >= 1 ? CH(1//2) * ch[2] : zero(CH)
  for m in 1:div(N, 2)
    y -= CH(Singular.bernoulli(2m) // (2m)) * ch[2m + 1]
  end

  # y has no constant term and the Chow ring vanishes above degree N,
  # so exp(y) = sum_{k = 0}^N y^k / k!
  td = one(CH)
  term = one(CH)
  for k in 1:N
    term = __simplify(CH(1//k) * term * y)
    td += term
  end
  setfield!(M.chow, :todd, __simplify(td))
  return M.chow.todd
end

"""
    integral(M::QuiverModuliSpace, x; unsafe::Bool=false)

Integrate a Chow class `x` over the fundamental class of `M`,
i.e. return the coefficient of the point class in the top-degree component
of `x`.

# Examples

The class of a hyperplane in `P^3` integrates to its degree:
```jldoctest
julia> M = QuiverModuliSpace(kronecker_quiver(4), [1, 1]); chow_ring(M);

julia> H = chern_class_line_bundle(M, [1, -1]);

julia> integral(M, H^3)
1
```
"""
function integral(M::QuiverModuliSpace, x; unsafe::Bool=false)
  n = dimension(M)
  top = Singular.div(__homogeneous_components(M, x, n), point_class(M; unsafe=unsafe))
  return Singular.constant_coefficient(top)
end

"""
    euler_characteristic(M::QuiverModuliSpace, f)
    euler_characteristic(F::Bundle)

Compute the Euler characteristic of the vector bundle whose Chern character is `f`
(resp. of the bundle `F`), via the Hirzebruch-Riemann-Roch theorem.

# Examples

The Euler characteristic of ``\\mathcal{O}(i)`` on the projective line:

```jldoctest
julia> Q = kronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> L = chern_character_line_bundle(M, [1, -1]);

julia> [euler_characteristic(M, L^i) for i in 0:5]
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

julia> [euler_characteristic(M, L^i) for i in 0:5]
6-element Vector{Singular.n_Q}:
 1
 20
 148
 664
 2206
 5999
```

This method can be used to compute the Euler characteristic of any bundle `F`:

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> U1 = universal_bundle(M, 1);

julia> euler_characteristic(U1)
0
```
"""
euler_characteristic(M::QuiverModuliSpace, f) = integral(M, f * todd_class(M))
euler_characteristic(F::Bundle) = euler_characteristic(variety(F), chern_character(F))
chi(F::Bundle) = euler_characteristic(F)

"""
Takes a quotient ring R/I and returns the projection map from R to R/I.
For internal use only.
"""
function __projection_to_quotient_ring(A)
  I = quotient_ideal(A)
  R = base_ring(I)
  return Singular.AlgebraHomomorphism(R, A, gens(A))
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
  deg = __chow_degrees(M.d)
  exp = first(Singular.AbstractAlgebra.exponent_vectors(f))
  @assert size(deg) == size(exp)
  return exp' * deg
end

"""
    __chow_degrees(d)

Compute the vector of degrees for the variables of a Chow ring.

For internal use only.
"""
function __chow_degrees(d::AbstractVector{Int})
  return vcat([collect(1:di) for di in d if di > 0]...)
end

"""
    __var_indices(d)

Compute the vector whose `i`-th entry is the number of variables of a Chow ring
attached to the vertices before `i`, so that the variable ``x_{i,j}`` has index
`__var_indices(d)[i] + j`.

For internal use only.
"""
function __var_indices(d::AbstractVector{Int})
  return [sum(d[1:(i - 1)]; init=0) for i in eachindex(d)]
end
