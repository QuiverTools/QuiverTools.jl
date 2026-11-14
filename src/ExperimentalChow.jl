# EXPERIMENTAL, LLM-SUGGESTED implementation. The algebraic reduction and its
# implementation need careful independent review. See docs/src/experimental-chow.md.
# The existing :weyl algorithm remains the default and the comparison reference.

function __check_chow_algorithm(algorithm::Symbol)
  algorithm in (:weyl, :experimental_parabolic) ||
    throw(ArgumentError("unknown Chow algorithm: $algorithm"))
end

function __chow_swap(f, i, j)
  context = Singular.MPolyBuildCtx(parent(f))
  for (c, powers) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    exponents = collect(powers)
    exponents[i], exponents[j] = exponents[j], exponents[i]
    Singular.push_term!(context, c, exponents)
  end
  Singular.finish(context)
end

# (f-s_i(f))/(x_(i+1)-x_i), evaluated monomial by monomial. Combine duplicate
# exponents before MPolyBuildCtx: Singular's builder requires distinct monomials.
function __chow_divided_difference(f, i, j)
  coefficients = Dict{Tuple,Rational{BigInt}}()
  for (c, powers) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    coefficient = BigInt(Singular.numerator(c)) // BigInt(Singular.denominator(c))
    a, b = powers[i], powers[j]
    for k in 0:(abs(a - b) - 1)
      exponents = collect(powers)
      if a > b
        exponents[i], exponents[j] = a - 1 - k, b + k
      else
        exponents[i], exponents[j] = a + k, b - 1 - k
      end
      key = Tuple(exponents)
      coefficients[key] =
        get(coefficients, key, big(0) // big(1)) +
        (a > b ? -coefficient : coefficient)
    end
  end
  context = Singular.MPolyBuildCtx(parent(f))
  for (powers, coefficient) in coefficients
    iszero(coefficient) ||
      Singular.push_term!(context, Singular.QQ(coefficient), collect(powers))
  end
  Singular.finish(context)
end

function __chow_divided_alternant(f, d)
  offset = 0
  for rank in d
    for top in (rank - 1):-1:1, k in 1:top
      f = __chow_divided_difference(f, offset + k, offset + k + 1)
      iszero(f) && return f
    end
    offset += rank
  end
  f
end

# Constructive fundamental theorem of symmetric polynomials, separately in
# each vertex block. A leading monomial determines a product of elementary
# symmetric polynomials, so no elimination/preimage computation is needed.
function __chow_elementary_expression(f, d, A, targets)
  variables = gens(A)
  out = A(0)
  while !iszero(f)
    terms = collect(zip(Singular.coefficients(f), Singular.exponent_vectors(f)))
    coefficient, powers = terms[argmax(Tuple(e) for (_, e) in terms)]
    elementary, answer = parent(f)(coefficient), A(coefficient)
    offset = 0
    for rank in d
      for k in 1:rank
        next = k == rank ? 0 : powers[offset + k + 1]
        power = powers[offset + k] - next
        power >= 0 || error("polynomial is not symmetric in its vertex block")
        elementary *= targets[offset + k]^power
        answer *= variables[offset + k]^power
      end
      offset += rank
    end
    f -= elementary
    out += answer
  end
  out
end

function __chow_rectangle_partitions(rows, columns)
  partitions = Vector{Int}[]
  lambda = zeros(Int, rows)
  function visit(at, upper)
    if at > rows
      push!(partitions, copy(lambda))
      return nothing
    end
    for value in 0:upper
      lambda[at] = value
      visit(at + 1, value)
    end
  end
  visit(1, columns)
  partitions
end

function __chow_parabolic_monomials(d, e)
  blocks = map(eachindex(d)) do i
    a, b = e[i], d[i] - e[i]
    [
      vcat([lambda[j] + a - j for j in 1:a], [b - j for j in 1:b])
      for lambda in __chow_rectangle_partitions(a, b)
    ]
  end
  [reduce(vcat, parts) for parts in Iterators.product(blocks...)]
end

function __chow_parabolic_alternant(f, monomial, d, e)
  split = reduce(vcat, ([e[i], d[i] - e[i]] for i in eachindex(d)))
  value = f * __chow_divided_alternant(monomial, split)
  offset = 0
  for i in eachindex(d)
    a, b = e[i], d[i] - e[i]
    # This is a reduced word for w0(W)w0(H); operators act in reverse order.
    word = [k for j in 1:a for k in (b + j - 1):-1:j]
    for k in reverse(word)
      value = __chow_divided_difference(value, offset + k, offset + k + 1)
      iszero(value) && return value
    end
    offset += d[i]
  end
  value
end

function __chow_graded_dimensions(ideal, weights, dimension)
  leading = [
    collect(first(Singular.exponent_vectors(g))) for g in gens(ideal) if !iszero(g)
  ]
  [
    count(_weighted_monomials(weights, degree)) do powers
      !any(all(powers[k] >= e[k] for k in eachindex(weights)) for e in leading)
    end for degree in 0:dimension
  ]
end

function __chow_ring_parabolic(Q, d, theta, chi; adaptive=true, verbose=false)
  Qs, _ = support_subquiver(Q, d)
  is_acyclic(Qs) ||
    throw(ArgumentError("the experimental backend requires an acyclic support"))
  dimension = 1 - euler_form(Q, d, d)
  dimension >= 0 ||
    throw(ArgumentError("the experimental backend requires a nonempty stable moduli space"))
  n = length(d)
  R, roots = Singular.polynomial_ring(Singular.QQ, ["xi$i$j" for i in 1:n for j in 1:d[i]])
  A, variables = Singular.polynomial_ring(
    Singular.QQ, ["x$i$j" for i in 1:n for j in 1:d[i]]
  )
  offsets = [sum(d[1:(i - 1)]; init=0) for i in 1:n]
  root(i, j) = roots[offsets[i] + j]
  targets = [
    symmetric_polynomial(k)([root(i, j) for j in 1:d[i]]) for i in 1:n for k in 1:d[i]
  ]
  inclusion = Singular.AlgebraHomomorphism(A, R, targets)
  forbidden = all_destabilizing_subdimension_vectors(d, theta)
  filter!(e -> !any(partial_order(Q, f, e) for f in forbidden if f != e), forbidden)
  root_degree = sum(r * (r - 1) ÷ 2 for r in d)
  weights = [j for i in eachindex(d) for j in 1:d[i]]
  relations = Singular.spoly{Singular.n_Q}[
    sum(chi[i] * variables[offsets[i] + 1] for i in support(d))]
  # Valid for a smooth projective moduli space of this dimension. Include the
  # minimal weighted monomial border, so every larger degree vanishes as well.
  for degree in (dimension + 1):(dimension + maximum(weights))
    for powers in _weighted_monomials(weights, degree)
      all(powers[k] == 0 || degree - weights[k] <= dimension for k in eachindex(weights)) ||
        continue
      push!(relations, prod(variables[k]^powers[k] for k in eachindex(variables)))
    end
  end
  expected = Int[]
  if adaptive
    vertices = support(d)
    betti = betti_numbers(QuiverModuliSpace(Qs, d[vertices], theta[vertices]))
    length(betti) == 2dimension + 1 && all(iszero, betti[2:2:end]) ||
      error("the Hilbert certificate requires the projective Hodge--Tate Betti grading")
    expected = betti[1:2:end]
  end
  polynomials = Dict{Int,Any}()
  seen = Set{String}()
  ideal = nothing
  certificate = nothing
  # The transfer can lower degree, so relevant Schur degrees can exceed dim X.
  last_stage = maximum(
    (
      min(sum(e[i] * (d[i] - e[i]) for i in 1:n),
        dimension + root_degree -
        sum(Q.adjacency[i, j] * e[i] * (d[j] - e[j]) for i in 1:n for j in 1:n) -
        sum(binomial(e[i], 2) + binomial(d[i] - e[i], 2) for i in 1:n))
      for e in forbidden
    ); init=0)
  for stage in 0:(adaptive ? last_stage : 0)
    for (index, e) in enumerate(forbidden)
      degree = sum(Q.adjacency[i, j] * e[i] * (d[j] - e[j]) for i in 1:n for j in 1:n)
      parabolic_degree = sum(binomial(e[i], 2) + binomial(d[i] - e[i], 2) for i in 1:n)
      degree + parabolic_degree - root_degree <= dimension || continue
      for powers in __chow_parabolic_monomials(d, e)
        degree + sum(powers) - root_degree <= dimension || continue
        adaptive && sum(powers) - parabolic_degree != stage && continue
        f = get!(polynomials, index) do
          prod(
            (root(j, s) - root(i, r))^Q.adjacency[i, j]
            for i in 1:n for j in 1:n if Q.adjacency[i, j] > 0
            for r in 1:e[i] for s in (e[j] + 1):d[j]; init=R(1))
        end
        monomial = prod(roots[k]^powers[k] for k in eachindex(roots); init=R(1))
        relation = __chow_parabolic_alternant(f, monomial, d, e)
        iszero(relation) && continue
        relation *= inv(Singular.leading_coefficient(relation))
        key = string(relation)
        key in seen && continue
        push!(seen, key)
        push!(relations, __chow_elementary_expression(relation, d, A, targets))
      end
    end
    ideal = Singular.std(Singular.Ideal(A, relations))
    if adaptive
      actual = __chow_graded_dimensions(ideal, weights, dimension)
      all(actual .>= expected) ||
        error("valid relations produced fewer classes than Reineke predicts")
      verbose && @info "Experimental parabolic Chow certificate" stage actual expected
      if actual == expected
        certificate = (
          schur_degree=stage, graded_dimensions=actual, reineke_even_betti=expected
        )
        break
      end
      stage == last_stage &&
        error("full parabolic relations failed the independent Hilbert check")
    end
  end
  (presentation=(Singular.QuotientRing(A, ideal), R, inclusion), certificate=certificate)
end
