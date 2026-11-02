###############################################################################
# Quasimap quantum cohomology of acyclic quiver moduli.
# Zhang--Zhou, arXiv:2511.09875v2, Theorem 2.22.
# See the manual for the finite cut-degree reduction and completion convention.
###############################################################################

"""
    QuasimapQuantumCohomology

A polynomial presentation whose completion in the effective Novikov monoid is
quasimap quantum cohomology. Construct it with [`quasimap_quantum_cohomology`](@ref).

Use [`quasimap_presentation`](@ref), [`quantum_tautological_classes`](@ref),
[`novikov_variables`](@ref), and [`novikov_degrees`](@ref) to access the algebra.
The polynomial presentation retains nilpotents and the Weyl discriminant locus.
Specializing Novikov variables is a polynomial operation, not evaluation of an
arbitrary formal series. No identification with Gromov--Witten theory is implicit;
see [`gw_comparison`](@ref).
"""
struct QuasimapQuantumCohomology
  moduli::QuiverModuliSpace
  ring::Singular.PolyRing{Singular.n_Q}
  _ambient::Singular.PolyRing{Singular.n_Q}
  _relations::Singular.sideal{Singular.spoly{Singular.n_Q}}
  _monoid_relations::Vector{Singular.spoly{Singular.n_Q}}
  _chi::Vector{BigInt}
  _degrees::Vector{Vector{BigInt}}
  _anticanonical_degrees::Vector{BigInt}
end

function show(io::IO, A::QuasimapQuantumCohomology)
  print(io, "Quasimap quantum cohomology presentation with ",
    sum(A.moduli.d), " tautological generators and ", length(A._degrees),
    " Novikov generators\n  Gromov--Witten comparison: ", gw_comparison(A),
    "\n  The Novikov completion is understood.")
end

# Coordinates on Z^vertices / Z*d. The first row also supplies a Bezout
# normalization, and avoids floating point nullspaces or a choice of a vertex
# of dimension one (which need not exist).
function __qm_degree_lattice(d)
  n = length(d)
  U = [BigInt(i == j) for i in 1:n, j in 1:n]
  v = BigInt.(d)
  for i in 2:n
    a, b = v[1], v[i]
    iszero(b) && continue
    g, s, t = gcdx(a, b)
    row1, rowi = copy(U[1, :]), copy(U[i, :])
    U[1, :] = s .* row1 + t .* rowi
    U[i, :] = (-b ÷ g) .* row1 + (a ÷ g) .* rowi
    v[1], v[i] = g, 0
  end
  v[1] == 1 || throw(ArgumentError("the dimension vector must be indivisible"))
  return U
end

# A torus degree is a potential on the expanded vertices, modulo common
# translation. Its generic section uses exactly the arrows of nonnegative
# degree. It is effective iff that support is theta-stable.
function __qm_cut_data(M)
  vertices = [i for i in eachindex(M.d) for _ in 1:M.d[i]]
  theta = BigInt.(M.theta[vertices])
  r = length(vertices)
  edges = [
    (i, j, M.Q.adjacency[vertices[i], vertices[j]])
    for i in 1:r for j in 1:r if M.Q.adjacency[vertices[i], vertices[j]] > 0
  ]
  forbidden = [
    collect(S) for S in subsets(1:r)
                   if !isempty(S) && length(S) < r && sum(theta[S]) > 0
  ]
  crossing = [[(i, j) for (i, j, _) in edges if i in S && !(j in S)]
              for S in forbidden]
  cuts = Vector{Int}[]
  for S in forbidden
    k = [-Int(i in S) for i in 1:r]
    if all(es -> any(k[j] >= k[i] for (i, j) in es), crossing)
      push!(cuts, k)
    end
  end
  return vertices, edges, cuts
end

# A Weyl orbit of binary cuts is determined by the number of selected roots
# in each vertex block. Keep the first cut in each orbit, preserving degree order.
function __qm_cut_orbits(d, cuts)
  offsets = cumsum([0; d])
  indices = unique(eachindex(cuts)) do i
    k = cuts[i]
    Tuple(sum(k[j] for j in (offsets[v] + 1):offsets[v + 1]; init=0)
          for v in eachindex(d))
  end
  return cuts[indices]
end

# Factor an effective monoid element using a finite set of generators. The
# positive stability height makes this a terminating exact integer search.
function __qm_factor(target, height, degrees, heights)::Union{Nothing,Vector{Int}}
  return __qm_factor_tuple(Tuple(target), height, degrees, heights)
end

function __qm_factor_tuple(target::NTuple{N,BigInt}, height, degrees, heights) where {N}
  memo = Dict{Tuple{Int,NTuple{N,BigInt}},Union{Nothing,Vector{Int}}}()
  return __qm_factor_search(target, height, length(degrees), degrees, heights, memo)
end

function __qm_factor_search(v::NTuple{N,BigInt}, h, j, degrees, heights, memo) where {N}
  all(iszero, v) && return zeros(Int, j)
  j == 0 && return nothing
  key = (j, v)
  haskey(memo, key) && return memo[key]
  for power in 0:Int(h ÷ heights[j])
    remainder = ntuple(i -> v[i] - power * degrees[j][i], Val(N))
    answer = __qm_factor_search(remainder, h - power * heights[j], j - 1,
      degrees, heights, memo)
    if !isnothing(answer)
      memo[key] = [answer; power]
      return memo[key]
    end
  end
  memo[key] = nothing
  return nothing
end

function __qm_monoid(M, vertices, cuts, U)
  representatives = [
    BigInt[
      sum(k[j] for j in eachindex(vertices) if vertices[j] == i; init=0)
      for i in eachindex(M.d)
    ] for k in cuts
  ]
  coordinates = [collect((U * e)[2:end]) for e in representatives]
  heights = [-sum(BigInt.(M.theta) .* e) for e in representatives]
  order = sortperm(eachindex(cuts); by=i -> (heights[i], Tuple(coordinates[i])))
  degrees, selected, selected_heights = Vector{BigInt}[], Int[], BigInt[]
  for i in order
    if isnothing(__qm_factor(coordinates[i], heights[i], degrees, selected_heights))
      push!(degrees, coordinates[i])
      push!(selected, i)
      push!(selected_heights, heights[i])
    end
  end
  powers = [
    __qm_factor(c, h, degrees, selected_heights)
    for (c, h) in zip(coordinates, heights)
  ]
  return representatives[selected], degrees, powers
end

# The kernel of the monoid map is essential: distinct products of Novikov
# generators may represent the same integral degree. Laurent coordinates are
# auxiliary to this kernel computation; no Novikov parameter is inverted.
function __qm_monoid_relations(R, qvars, coordinates)
  isempty(qvars) && return typeof(R(0))[]
  rank = length(first(coordinates))
  # Independent columns have no binomial relations.
  C = hcat(coordinates...)
  rankx(C) == length(qvars) && return typeof(R(0))[]
  T, t = polynomial_ring(Singular.QQ,
    ["t$i" for i in 1:rank] ∪ ["u$i" for i in 1:rank])
  Q, q = polynomial_ring(Singular.QQ, ["q$i" for i in eachindex(qvars)])
  images = [
    prod(t[i]^Int(max(c[i], 0)) * t[rank + i]^Int(max(-c[i], 0))
         for i in 1:rank; init=T(1)) for c in coordinates
  ]
  map = Singular.AlgebraHomomorphism(Q, T, images)
  kernel = preimage(map, std(Singular.Ideal(T, [t[i] * t[rank + i] - 1
                                                for i in 1:rank])))
  inclusion = Singular.AlgebraHomomorphism(Q, R, qvars)
  return [inclusion(f) for f in gens(kernel)]
end

function __qm_weyl_data(d, R, roots)
  # A reduced word for the longest element in each symmetric-group block.
  word = Int[]
  offset = 0
  for n in d
    for j in 1:(n - 1), i in j:-1:1
      push!(word, offset + i)
    end
    offset += n
  end
  differences = [(i, roots[i + 1] - roots[i]) for i in word]
  bounds = [0:(n - j) for n in d for j in 1:n]
  basis = [
    prod(roots[i]^e[i] for i in eachindex(roots); init=R(1))
    for e in Iterators.product(bounds...)
  ]
  return differences, basis
end

# The longest divided difference is Alt(f)/prod(x_j-x_i), with i < j in
# each block. Every division is exact; no root or discriminant is inverted.
function __qm_divided_alternant(f, differences)
  R = parent(f)
  for (i, divisor) in differences
    iszero(f) && return f
    context = Singular.MPolyBuildCtx(R)
    for (c, e) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
      e[i], e[i + 1] = e[i + 1], e[i]
      Singular.push_term!(context, c, e)
    end
    f = Singular.divexact(f - Singular.finish(context), divisor)
  end
  return f
end

"""
    quasimap_quantum_cohomology(M::QuiverModuliSpace; chi=nothing, verbose=false)

Compute an exact polynomial presentation of quasimap quantum cohomology over
the effective integral Novikov monoid. Its Novikov completion is the algebra
of Zhang--Zhou, [arXiv:2511.09875v2](https://arxiv.org/abs/2511.09875v2), Theorem 2.22.
This is not a degree truncation or a computation on the regular Bethe locus.

The quiver must be acyclic, the dimension vector nonnegative and indivisible,
the stability parameter orthogonal to it and coprime, and the moduli nonempty.
These hypotheses give a smooth projective fine moduli space. Vertices of
dimension zero are allowed. `chi` is an integral vector with `dot(chi, M.d) == 1`;
its default is a Bezout vector. The normalization is internal to this result
and does not change `M` or its cached Chow ring.

The returned [`QuasimapQuantumCohomology`](@ref) provides a Singular quotient
ring with Chern generators `xij` followed by minimal monoid generators `qk`.
The latter need not be an independent lattice basis; all monoid relations are
included. The polynomial model can have components away from the formal
neighborhood of `q = 0`; its specialization at nonzero parameters is not, in
general, a specialization of the completed algebra.

The algorithm enumerates cuts of the expanded quiver (with `sum(M.d)` vertices),
retains one cut per Weyl orbit, and uses a coinvariant basis of size
`prod(factorial, M.d)`. Each divided alternant is computed by a reduced word
of adjacent divided differences, without enumerating the Weyl group.
The cut enumeration is still exponential and the basis size factorial.
No extra dependencies or user-supplied effective cone are required.

# Examples

The projective plane, including its classical limit:

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(3), [1, 1]));

julia> x = quantum_tautological_classes(A); q = only(novikov_variables(A));

julia> h = x[2][1] - x[1][1]; iszero(quantum_normal_form(A, h^3 - q))
true

julia> gw_comparison(A)
:equal

julia> C = specialize_quantum_parameters(A, [0]);

julia> length(QuiverTools.gens(QuiverTools.Singular.kbase(QuiverTools.quotient_ideal(C))))
3
```
"""
function quasimap_quantum_cohomology(
  M::QuiverModuliSpace; chi=nothing, verbose::Bool=false
)
  any(<(0), M.d) && throw(ArgumentError("the dimension vector must be nonnegative"))
  sum(M.d) > 0 || throw(ArgumentError("the dimension vector must be nonzero"))
  all(>=(0), M.Q.adjacency) ||
    throw(ArgumentError("arrow multiplicities must be nonnegative"))
  is_acyclic(M.Q) || throw(ArgumentError("the quiver must be acyclic"))
  M.denom === sum || throw(ArgumentError("standard King stability is required (denom=sum)"))
  sum(BigInt.(M.d) .* M.theta) == 0 ||
    throw(
      ArgumentError("the stability parameter must be orthogonal to the dimension vector")
    )
  U = __qm_degree_lattice(M.d)
  is_coprime(M) || throw(ArgumentError("the stability parameter must be coprime"))
  is_nonempty(M) || throw(ArgumentError("the moduli space must be nonempty"))
  isnothing(chi) || chi isa AbstractVector ||
    throw(ArgumentError("chi must be an integral vector"))
  normalization = isnothing(chi) ? collect(U[1, :]) : collect(chi)
  all(x -> x isa Integer, normalization) && length(normalization) == length(M.d) &&
  sum(BigInt.(normalization) .* M.d) == 1 ||
    throw(ArgumentError("chi must be an integral vector with dot(chi, M.d) == 1"))
  normalization = BigInt.(normalization)

  vertices, edges, cuts = __qm_cut_data(M)
  cut_count = length(cuts)
  cuts = __qm_cut_orbits(M.d, cuts)
  representatives, coordinates, powers = __qm_monoid(M, vertices, cuts, U)
  r, m = length(vertices), length(coordinates)
  verbose &&
    @info "Effective cut degrees" cuts=cut_count cut_orbits=length(cuts) novikov_generators=m
  R, vars = polynomial_ring(Singular.QQ,
    ["xi$i$j" for i in eachindex(M.d) for j in 1:M.d[i]] ∪
    ["q$i" for i in 1:m])
  roots, qvars = vars[1:r], vars[(r + 1):end]
  differences, basis = __qm_weyl_data(M.d, R, roots)
  alternants = typeof(R(0))[]
  root_coefficients = [M.d[i] + 1 - 2j for i in eachindex(M.d) for j in 1:M.d[i]]
  for (k, exponent) in zip(cuts, powers)
    positive, negative = R(1), R(1)
    for (i, j, multiplicity) in edges
      pairing = k[j] - k[i]
      if pairing > 0
        positive *= (roots[j] - roots[i])^multiplicity
      elseif pairing < 0
        negative *= (roots[j] - roots[i])^multiplicity
      end
    end
    novikov = prod(qvars[i]^exponent[i] for i in 1:m; init=R(1))
    root_sign = isodd(sum(root_coefficients .* k)) ? -1 : 1
    binomial = positive - root_sign * novikov * negative
    for b in basis
      f = __qm_divided_alternant(b * binomial, differences)
      iszero(f) || push!(alternants, f)
    end
  end

  B, bvars = polynomial_ring(Singular.QQ,
    ["x$i$j" for i in eachindex(M.d) for j in 1:M.d[i]] ∪
    ["q$i" for i in 1:m])
  images = [
    symmetric_polynomial(j)(roots[findall(==(i), vertices)])
    for i in eachindex(M.d) for j in 1:M.d[i]
  ]
  inclusion = Singular.AlgebraHomomorphism(B, R, [images; qvars])
  verbose && @info "Descending divided alternants" relations=length(alternants)
  root_ideal = Singular.Ideal(R, isempty(alternants) ? [R(0)] : alternants)
  relations = collect(gens(preimage(inclusion, root_ideal)))
  push!(
    relations,
    sum(
      normalization[i] * bvars[findfirst(==(i), vertices)]
      for i in eachindex(M.d) if M.d[i] > 0; init=B(0)),
  )
  monoid_relations = __qm_monoid_relations(B, bvars[(r + 1):end], coordinates)
  append!(relations, monoid_relations)
  ideal = std(Singular.Ideal(B, relations))
  ring = first(Singular.QuotientRing(B, ideal))
  canonical = BigInt.(canonical_stability(M.Q, M.d))
  anticanonical = [-sum(canonical .* e) for e in representatives]
  return QuasimapQuantumCohomology(M, ring, B, ideal, monoid_relations, normalization,
    representatives, anticanonical)
end

"""
    quasimap_presentation(A::QuasimapQuantumCohomology)

Return the Singular quotient ring defining the polynomial model of `A`.
Its Novikov completion, rather than an implicit localization or truncation,
is quasimap quantum cohomology. The unit in this presentation is the polynomial `1`.
"""
quasimap_presentation(A::QuasimapQuantumCohomology) = A.ring

"""
    quantum_normal_form(A::QuasimapQuantumCohomology, f)

Reduce the presentation element `f` to its exact Gröbner normal form in `A`.
The result has the same parent. Use this before testing a polynomial relation
for equality or zero: Singular quotient expressions are not always reduced
automatically. Elements of other presentations are rejected.
"""
function quantum_normal_form(A::QuasimapQuantumCohomology, f)
  f isa Singular.spoly && parent(f) === A.ring ||
    throw(ArgumentError("f must belong to the quasimap presentation"))
  return __simplify(f)
end

"""
    novikov_variables(A::QuasimapQuantumCohomology)

Return the minimal effective Novikov monoid generators, as elements of the
polynomial presentation. They may satisfy binomial relations. Their order
matches [`novikov_degrees`](@ref).
"""
function novikov_variables(A::QuasimapQuantumCohomology)
  return gens(A.ring)[(sum(A.moduli.d) + 1):end]
end

"""
    novikov_degrees(A::QuasimapQuantumCohomology)

Return integral representatives of the Novikov generator degrees in the
vertex lattice modulo `Z*M.d`, in the same order as [`novikov_variables`](@ref).
Entries at vertices outside the support are zero.
Representatives differing by a multiple of `M.d` describe the same degree.
Returns a copy so callers cannot mutate the presentation's metadata.
"""
novikov_degrees(A::QuasimapQuantumCohomology) = deepcopy(A._degrees)

"""
    quantum_tautological_classes(A::QuasimapQuantumCohomology)

Return a vector of vectors: entry `i` contains the Chern generators of the
normalized universal bundle at vertex `i`. These represent *quantum*
tautological classes; they need not coincide with classical Chern classes,
even when the quasimap and Gromov--Witten products agree.
"""
function quantum_tautological_classes(A::QuasimapQuantumCohomology)
  xs = gens(A.ring)
  return [
    xs[(sum(A.moduli.d[1:(i - 1)]) + 1):sum(A.moduli.d[1:i])]
    for i in eachindex(A.moduli.d)
  ]
end

"""
    gw_comparison(A::QuasimapQuantumCohomology)

Return `:equal` when every nonzero effective quasimap degree has
anticanonical degree at least two, and `:unknown` otherwise. The test is made
on the exact finite generating set of the effective monoid, not just on
degrees of actual curves. In the first case the nonequivariant wall-crossing
corrections `[z*I_beta(z)]_nonnegative` vanish by grading, so the primary
three-point products agree over the completed Novikov ring.

`:unknown` is not evidence of inequality. This sufficient test does not
compute mirror transformations. It also does not assert that high-degree
quantum tautological classes are their classical restrictions, or justify
evaluation of arbitrary formal series at nonzero Novikov parameters.
"""
function gw_comparison(A::QuasimapQuantumCohomology)
  return all(>=(2), A._anticanonical_degrees) ? :equal : :unknown
end

"""
    specialize_quantum_parameters(A::QuasimapQuantumCohomology, values)

Specialize every Novikov generator of the *polynomial presentation* to the
corresponding rational value and return a Singular quotient ring. `values`
must satisfy the Novikov monoid relations. All zeros give the classical
limit; all ones always satisfy the monoid relations.

This does not evaluate formal power series. For a general stability
condition, components of the polynomial model away from the Novikov origin
can survive this operation. No discriminant or arrow weight is inverted.
"""
function specialize_quantum_parameters(A::QuasimapQuantumCohomology, values)
  length(values) == length(A._degrees) ||
    throw(ArgumentError("one rational value is required for each Novikov generator"))
  all(x -> x isa Integer || x isa Rational, values) ||
    throw(ArgumentError("Novikov values must be exact rational numbers"))
  qs = gens(A._ambient)[(sum(A.moduli.d) + 1):end]
  specialization = [qs[i] - A._ambient(values[i]) for i in eachindex(qs)]
  evaluate = Singular.AlgebraHomomorphism(A._ambient, A._ambient,
    [gens(A._ambient)[1:sum(A.moduli.d)]; [A._ambient(v) for v in values]])
  all(f -> iszero(evaluate(f)), A._monoid_relations) ||
    throw(ArgumentError("the values do not satisfy the Novikov monoid relations"))
  relations = [collect(gens(A._relations)); specialization]
  return first(
    Singular.QuotientRing(A._ambient, std(Singular.Ideal(A._ambient, relations)))
  )
end

"""
    quantum_multiplication_matrix(A::QuasimapQuantumCohomology, f, values)

Return `(matrix, basis, ring)` for multiplication by the presentation element
`f`, after [`specialize_quantum_parameters`](@ref) at exact rational `values`.
Columns are the coefficients of `f * basis[j]` in the returned standard
monomial basis. The matrix is a Nemo rational matrix. Nilpotents and repeated
eigenvalues are retained. An infinite-dimensional specialization is rejected.

The result concerns the polynomial specialization, with the same completion
and Gromov--Witten comparison qualifications as the presentation itself.

# Examples

For the projective line at `q = 1`, multiplication by the hyperplane exchanges
the two basis elements:

```jldoctest
julia> A = quasimap_quantum_cohomology(QuiverModuliSpace(kronecker_quiver(2), [1, 1]));

julia> x = quantum_tautological_classes(A); h = x[2][1] - x[1][1];

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [1]);

julia> matrix^2 == one(matrix)
true

julia> matrix, basis, ring = quantum_multiplication_matrix(A, h, [0]);

julia> iszero(matrix^2) && !iszero(matrix)
true
```
"""
function quantum_multiplication_matrix(A::QuasimapQuantumCohomology, f, values)
  f = quantum_normal_form(A, f)
  ring = specialize_quantum_parameters(A, values)
  ideal = quotient_ideal(ring)
  Singular.dimension(ideal) == 0 ||
    throw(ArgumentError("the specialization must be nonzero and finite dimensional"))
  raw_basis = collect(gens(Singular.kbase(ideal)))
  map = Singular.AlgebraHomomorphism(A.ring, ring, gens(ring))
  basis = [ring(b) for b in raw_basis]
  indices = Dict(
    Tuple(first(Singular.exponent_vectors(b))) => i
    for (i, b) in enumerate(basis)
  )
  N = Singular.Nemo
  matrix = N.zero_matrix(N.QQ, length(basis), length(basis))
  for (j, b) in enumerate(basis)
    product = __simplify(map(f) * b)
    for (c, e) in zip(Singular.coefficients(product), Singular.exponent_vectors(product))
      matrix[indices[Tuple(e)], j] = N.QQ(BigInt(numerator(c)), BigInt(denominator(c)))
    end
  end
  return matrix, basis, ring
end
