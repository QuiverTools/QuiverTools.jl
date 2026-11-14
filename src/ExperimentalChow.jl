###############################################################################
# Experimental parabolic presentation of the Chow ring.
# The reduction, sign conventions, and geometric hypotheses are documented in
# docs/src/experimental-chow.md. This LLM-suggested algorithm needs review.
###############################################################################

"""
    __chow_divided_difference(f::Singular.spoly{Singular.n_Q}, i::Int)

Apply the adjacent divided difference `(f - sᵢ(f))/(xᵢ₊₁ - xᵢ)`.
The variables are ordered as in `gens(parent(f))`.

Expand the quotient of each monomial directly, without polynomial division.
Combine equal exponent vectors before passing them to `MPolyBuildCtx`, which
requires distinct monomials. Coefficients stay in Singular's rational field.

For internal use only.

# Examples

```jldoctest
julia> R, (x, y) = QuiverTools.polynomial_ring(QuiverTools.Singular.QQ, ["x", "y"]);

julia> QuiverTools.__chow_divided_difference(x^2, 1) == -x - y
true

julia> iszero(QuiverTools.__chow_divided_difference(x*y, 1))
true
```
"""
function __chow_divided_difference(f::Singular.spoly{Singular.n_Q}, i::Int)
  j = i + 1
  coefficients = Dict{Vector{Int},Singular.n_Q}()
  for (coefficient, powers) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    a, b = powers[i], powers[j]
    a == b && continue
    signed_coefficient = a > b ? -coefficient : coefficient
    for k in 0:(abs(a - b) - 1)
      exponents = copy(powers)
      if a > b
        exponents[i], exponents[j] = a - 1 - k, b + k
      else
        exponents[i], exponents[j] = a + k, b - 1 - k
      end
      coefficients[exponents] =
        get(coefficients, exponents, zero(Singular.QQ)) + signed_coefficient
    end
  end

  context = Singular.MPolyBuildCtx(parent(f))
  for (powers, coefficient) in coefficients
    iszero(coefficient) || Singular.push_term!(context, coefficient, powers)
  end
  return Singular.finish(context)
end

"""
    __chow_divided_alternant(f::Singular.spoly{Singular.n_Q}, d::AbstractVector{Int})

Divide the alternating sum over the product of the symmetric groups of ranks
`d` by its Vandermonde polynomial. Use the differences `x[j] - x[i]` for `i < j`
within each block.

Evaluate a reduced word for the longest permutation by adjacent divided
differences. Blocks of rank zero or one act as the identity.

For internal use only.
"""
function __chow_divided_alternant(
  f::Singular.spoly{Singular.n_Q}, d::AbstractVector{Int}
)
  offset = 0
  for rank in d
    for top in (rank - 1):-1:1, k in 1:top
      f = __chow_divided_difference(f, offset + k)
      iszero(f) && return f
    end
    offset += rank
  end
  return f
end

"""
    __chow_elementary_expression(f, d::AbstractVector{Int}, A, targets)

Express a block-symmetric polynomial `f` in the elementary symmetric generators
of `A`. The vector `targets` contains their images in `parent(f)`, in vertex
order and then increasing elementary degree within each vertex.

Cancel the lexicographically largest monomial using successive differences of
its exponents. The polynomial `f` must be symmetric in the consecutive blocks
of sizes `d`; otherwise this reduction need not apply.

For internal use only.
"""
function __chow_elementary_expression(
  f::Singular.spoly{Singular.n_Q}, d::AbstractVector{Int}, A, targets
)
  variables = gens(A)
  result = zero(A)
  while !iszero(f)
    coefficient, powers = argmax(
      term -> Tuple(last(term)), zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    )
    elementary, monomial = parent(f)(coefficient), A(coefficient)
    offset = 0
    for rank in d
      for k in 1:rank
        next_power = k == rank ? 0 : powers[offset + k + 1]
        power = powers[offset + k] - next_power
        power >= 0 || error("polynomial is not symmetric in its vertex block")
        elementary *= targets[offset + k]^power
        monomial *= variables[offset + k]^power
      end
      offset += rank
    end
    f -= elementary
    result += monomial
  end
  return result
end

"""
    __chow_rectangle_partitions(rows::Int, columns::Int)

List partitions in a `rows` by `columns` rectangle, padded with trailing zeros.
The empty partition is included, also when either rectangle dimension is zero.

For internal use only.

# Examples

```jldoctest
julia> QuiverTools.__chow_rectangle_partitions(2, 2)
6-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 0]
 [1, 1]
 [2, 0]
 [2, 1]
 [2, 2]
```
"""
function __chow_rectangle_partitions(rows::Int, columns::Int)
  partitions = Vector{Int}[]
  partition = zeros(Int, rows)
  function visit(row, upper)
    if row > rows
      push!(partitions, copy(partition))
      return nothing
    end
    for value in 0:upper
      partition[row] = value
      visit(row + 1, value)
    end
    return nothing
  end
  visit(1, columns)
  return partitions
end

"""
    __chow_parabolic_monomials(d::AbstractVector{Int}, e::AbstractVector{Int})

Iterate over exponent vectors whose parabolic divided alternants give a Schur
basis, up to a common sign. At vertex `i`, set `a = e[i]`, `b = d[i] - e[i]`, and
use exponents `(λ₁+a-1, …, λₐ, b-1, …, 0)` for partitions `λ` in the `a` by `b`
rectangle. The number of vectors is `prod(binomial.(d, e))`.

The Cartesian product over vertices is lazy; only the individual block bases
are stored. The input must satisfy `0 <= e[i] <= d[i]` at each vertex.

For internal use only.
"""
function __chow_parabolic_monomials(d::AbstractVector{Int}, e::AbstractVector{Int})
  blocks = map(eachindex(d)) do i
    a, b = e[i], d[i] - e[i]
    [
      vcat([partition[j] + a - j for j in 1:a], [b - j for j in 1:b])
      for partition in __chow_rectangle_partitions(a, b)
    ]
  end
  return (vcat(parts...) for parts in Iterators.product(blocks...))
end

"""
    __chow_parabolic_alternant(f, monomial, d::AbstractVector{Int}, e::AbstractVector{Int})

Compute the full divided alternant of `f * monomial`, using the invariance of
`f` under permutations within the blocks of ranks `e[i]` and `d[i] - e[i]`.

First apply the parabolic longest-word operator to `monomial`, then multiply
by `f` and apply the remaining divided differences. For a split `(a, b)`, the
remaining reduced word concatenates `(b+j-1, …, j)` for `j = 1:a`; operators
act in reverse word order. See the factorization proof in the experimental
Chow documentation.

For internal use only.
"""
function __chow_parabolic_alternant(
  f::Singular.spoly{Singular.n_Q}, monomial::Singular.spoly{Singular.n_Q},
  d::AbstractVector{Int}, e::AbstractVector{Int},
)
  split = [rank for i in eachindex(d) for rank in (e[i], d[i] - e[i])]
  result = f * __chow_divided_alternant(monomial, split)
  offset = 0
  for i in eachindex(d)
    a, b = e[i], d[i] - e[i]
    for j in a:-1:1, k in j:(b + j - 1)
      result = __chow_divided_difference(result, offset + k)
      iszero(result) && return result
    end
    offset += d[i]
  end
  return result
end

"""
    __chow_truncation_relations(A, weights::AbstractVector{Int}, dimension::Int)

Return the minimal monomial generators of the ideal of weighted degrees above
`dimension`. The positive entries of `weights` are the degrees of `gens(A)`.
A monomial belongs to this list precisely when removing any variable dividing
it lowers its degree to at most `dimension`.

These are valid Chow relations only under the geometric hypotheses of
`__chow_ring_parabolic`.

For internal use only.
"""
function __chow_truncation_relations(A, weights::AbstractVector{Int}, dimension::Int)
  variables = gens(A)
  return [
    prod(variables[k]^powers[k] for k in eachindex(variables); init=one(A))
    for degree in (dimension + 1):(dimension + maximum(weights))
    for powers in _weighted_monomials(weights, degree)
    if
    all(iszero(powers[k]) || degree - weights[k] <= dimension for k in eachindex(weights))
  ]
end

"""
    __chow_graded_dimensions(basis, weights::AbstractVector{Int}, dimension::Int)

Count standard monomials in weighted degrees `0:dimension` for a homogeneous
ideal. The input `basis` must already be a Gröbner basis. The first exponent
vector of each nonzero generator is its leading monomial in the Singular
ring's ordering; standard monomials are those divisible by none of them.

For internal use only.
"""
function __chow_graded_dimensions(basis, weights::AbstractVector{Int}, dimension::Int)
  leading = [first(Singular.exponent_vectors(g)) for g in gens(basis) if !iszero(g)]
  return [
    count(_weighted_monomials(weights, degree)) do powers
      !any(exponents -> all(powers .>= exponents), leading)
    end for degree in 0:dimension
  ]
end

"""
    __chow_parabolic_relations!(relations, seen, f, d, e, A, targets, dimension; schur_degree=nothing)

Append the relations contributed by the forbidden polynomial `f` of vector `e`.
Convert each nonzero divided alternant to elementary symmetric polynomials in
`A`, discarding relations of weighted degree above `dimension`.

With an integer `schur_degree`, use only Schur basis elements of that degree.
With `nothing`, use the whole basis in the relevant degrees. Normalize the
alternants to leading coefficient one before inserting them into `seen`;
stored polynomials must not be mutated afterwards.

For internal use only.
"""
function __chow_parabolic_relations!(
  relations, seen, f::Singular.spoly{Singular.n_Q},
  d::AbstractVector{Int}, e::AbstractVector{Int}, A, targets, dimension::Int;
  schur_degree::Union{Int,Nothing}=nothing,
)
  roots = gens(parent(f))
  forbidden_degree = Singular.total_degree(f)
  root_degree = sum(binomial(rank, 2) for rank in d)
  parabolic_degree = sum(binomial(e[i], 2) + binomial(d[i] - e[i], 2) for i in eachindex(d))
  for powers in __chow_parabolic_monomials(d, e)
    monomial_degree = sum(powers)
    forbidden_degree + monomial_degree - root_degree <= dimension || continue
    isnothing(schur_degree) || monomial_degree - parabolic_degree == schur_degree ||
      continue
    monomial = prod(roots[k]^powers[k] for k in eachindex(roots); init=one(parent(f)))
    relation = __chow_parabolic_alternant(f, monomial, d, e)
    iszero(relation) && continue
    relation *= inv(Singular.leading_coefficient(relation))
    relation in seen && continue
    push!(seen, relation)
    push!(relations, __chow_elementary_expression(relation, d, A, targets))
  end
  return nothing
end

"""
    __chow_ring_parabolic(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, chi::AbstractVector{Int}; adaptive::Bool=true, verbose::Bool=false)

Construct the parabolic presentation of a quiver moduli space's Chow ring.
The moduli space must be nonempty, smooth, projective and fine, with stable
equal to semistable and acyclic support. The public `chow_ring` method checks
stability and the normalization `dot(chi, d) == 1` before calling this routine.

Start with the normalization relation and the dimension truncation. With
`adaptive=true`, add tautological relations in increasing Schur degree until
standard-monomial counts equal independently computed even Betti numbers.
With `adaptive=false`, add every parabolic relation of relevant degree.

Return a named tuple with fields:
- `presentation`: the same `(QuotientRing, root_ring, inclusion)` tuple returned
  by the quiver-level `chow_ring` method.
- `certificate`: `nothing` in full mode; otherwise the stopping `schur_degree`,
  `graded_dimensions`, and independent `reineke_even_betti`.

The transfer calculation, degree bound, and Hilbert comparison are proved in
`docs/src/experimental-chow.md`.

For internal use only.
"""
function __chow_ring_parabolic(
  Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, chi::AbstractVector{Int};
  adaptive::Bool=true, verbose::Bool=false,
)
  Qs, ds = support_subquiver(Q, d)
  is_acyclic(Qs) ||
    throw(ArgumentError("the experimental backend requires an acyclic support"))
  dimension = 1 - euler_form(Q, d, d)
  dimension >= 0 ||
    throw(ArgumentError("the experimental backend requires a nonempty stable moduli space"))

  R, roots = Singular.polynomial_ring(
    Singular.QQ, ["xi$i$j" for i in eachindex(d) for j in 1:d[i]]
  )
  A, variables = Singular.polynomial_ring(
    Singular.QQ, ["x$i$j" for i in eachindex(d) for j in 1:d[i]]
  )
  offsets = cumsum(d) .- d
  root(i, j) = roots[offsets[i] + j]
  targets = [
    symmetric_polynomial(k)([root(i, j) for j in 1:d[i]])
    for i in eachindex(d) for k in 1:d[i]
  ]
  inclusion = Singular.AlgebraHomomorphism(A, R, targets)

  forbidden = all_destabilizing_subdimension_vectors(d, theta)
  filter!(e -> !any(partial_order(Q, f, e) for f in forbidden if f != e), forbidden)
  root_degree = sum(binomial(rank, 2) for rank in d)
  # For a Schur class of degree s, the resulting relation has degree
  # deg(f_e) + deg(Delta_H) + s - deg(Delta_W).
  schur_bounds = map(forbidden) do e
    forbidden_degree = sum(
      Q.adjacency[i, j] * e[i] * (d[j] - e[j]) for i in eachindex(d), j in eachindex(d)
    )
    parabolic_degree = sum(
      binomial(e[i], 2) + binomial(d[i] - e[i], 2) for i in eachindex(d)
    )
    min(sum(e .* (d - e)), dimension + root_degree - forbidden_degree - parabolic_degree)
  end

  weights = [j for rank in d for j in 1:rank]
  relations = __chow_truncation_relations(A, weights, dimension)
  push!(relations, sum(chi[i] * variables[offsets[i] + 1] for i in support(d)))
  polynomials = Dict{Int,Singular.spoly{Singular.n_Q}}()
  seen = Set{Singular.spoly{Singular.n_Q}}()

  expected = Int[]
  if adaptive
    betti = betti_numbers(QuiverModuliSpace(Qs, ds, theta[support(d)]))
    length(betti) == 2dimension + 1 && all(iszero, betti[2:2:end]) ||
      error("the Hilbert comparison requires the projective Hodge--Tate Betti grading")
    expected = betti[1:2:end]
  end

  stages = adaptive ? (0:maximum(schur_bounds; init=0)) : (nothing,)
  for stage in stages
    for (index, e) in enumerate(forbidden)
      schur_bounds[index] < (isnothing(stage) ? 0 : stage) && continue
      f = get!(polynomials, index) do
        prod(
          (root(j, s) - root(i, r))^Q.adjacency[i, j]
          for i in eachindex(d), j in eachindex(d) if Q.adjacency[i, j] > 0
          for r in 1:e[i], s in (e[j] + 1):d[j]; init=one(R)
        )
      end
      __chow_parabolic_relations!(
        relations, seen, f, d, e, A, targets, dimension; schur_degree=stage
      )
    end
    basis = Singular.std(Singular.Ideal(A, relations))
    if !adaptive
      return (
        presentation=(Singular.QuotientRing(A, basis), R, inclusion), certificate=nothing
      )
    end

    actual = __chow_graded_dimensions(basis, weights, dimension)
    all(actual .>= expected) ||
      error("valid relations produced fewer classes than Reineke predicts")
    verbose && @info "Parabolic Chow graded dimensions" schur_degree = stage actual expected
    if actual == expected
      certificate = (
        schur_degree=stage, graded_dimensions=actual, reineke_even_betti=expected
      )
      return (
        presentation=(Singular.QuotientRing(A, basis), R, inclusion),
        certificate=certificate,
      )
    end
  end
  error("full parabolic relations failed the independent Hilbert comparison")
end
