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

function symmetric_polynomial(degree::Int)
  function f(vars)
    R = parent(first(vars))
    return sum(
      (prod(subset; init=R(1)) for subset in IterTools.subsets(vars, degree));
      init=R(0),
    )
  end

  return f
end

function __permutation_sign(sigma)
  inversions = 0
  for i in 1:(length(sigma) - 1), j in (i + 1):length(sigma)
    inversions += sigma[i] > sigma[j]
  end
  return iseven(inversions) ? 1 : -1
end

function __permutation_sign_product(sigma::Tuple)
  return prod(__permutation_sign(s) for s in sigma; init=1)
end

function __chow_dimension(M::QuiverModuliSpace; unsafe::Bool=false)
  return unsafe ? 1 - euler_form(M.Q, M.d, M.d) : dimension(M)
end

function __truncate_by_degree(x, n::Int)
  R = parent(x)
  return Oscar.simplify(sum((x[i] for i in 0:n); init=R(0)))
end

function __chow_degrees(d::AbstractVector{Int})
  return vcat([collect(1:di) for di in d if di > 0]...)
end

function __resolve_linearization!(M::QuiverModuliSpace, chi)
  if chi isa UndefInitializer
    return linearization(M)
  end

  chi = coerce_vector(chi)
  chi' * M.d != 1 && throw(ArgumentError("$(collect(chi)) is not a linearization."))

  cached = M.linearization_cache[]
  if isnothing(cached) || cached != chi
    M.linearization_cache[] = chi
    M.variety_cache[] = nothing
  end

  return M.linearization_cache[]
end

function __validate_chow_inputs(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int};
  unsafe::Bool=false,
  verbose::Bool=false,
)
  if unsafe
    verbose && @warn "Unsafe computation."
    return nothing
  end

  has_properly_semistables(Q, d, theta) && throw(
    ArgumentError(
      "The quiver moduli problem has properly semistable representations, no description of the Chow ring is available."
    ),
  )
  !is_amply_stable(Q, d, theta) && throw(
    ArgumentError(
      "The quiver moduli problem is not amply stable, no description of the Chow ring is available."
    ),
  )

  return nothing
end

function __presentation_ring(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int};
  chi::AbstractVector{Int},
  verbose::Bool=false,
  unsafe::Bool=false,
)
  chi' * d != 1 && throw(ArgumentError("$(collect(chi)) is not a linearization."))
  __validate_chow_inputs(Q, d, theta; unsafe=unsafe, verbose=verbose)

  root_names = ["xi$i$j" for i in 1:n_vertices(Q) for j in 1:d[i]]
  R, root_vars = Oscar.polynomial_ring(Oscar.QQ, root_names)

  root_offset(i) = sum(d[1:(i - 1)]; init=0)
  xi(i, j) = root_vars[root_offset(i) + j]

  bounds = [0:(d[i] - nu) for i in 1:n_vertices(Q) for nu in 1:d[i]]
  function build_module_element(lambda)
    out = R(1)
    for i in support(d)
      offset = root_offset(i)
      for nu in 1:d[i]
        out *= xi(i, nu)^lambda[offset + nu]
      end
    end
    return out
  end
  module_basis = map(build_module_element, Iterators.product(bounds...))
  verbose && @info "base has $(length(module_basis)) elements"

  permutation_blocks = [Tuple.(collect(permutations(1:d[i]))) for i in 1:n_vertices(Q)]
  W = collect(Iterators.product(permutation_blocks...))
  permuted_indices = Dict{Tuple,Vector{Int}}(
    sigma => reduce(
      vcat,
      [Int[root_offset(i) + sigma[i][j] for j in 1:d[i]] for i in support(d)],
    ) for sigma in W
  )
  permute_vector(exponents, sigma) = exponents[permuted_indices[sigma]]

  function permute_polynomial(f, sigma)
    context = AbstractAlgebra.MPolyBuildCtx(parent(f))
    for (coeff, exponent) in zip(
      AbstractAlgebra.coefficients(f),
      AbstractAlgebra.exponent_vectors(f),
    )
      AbstractAlgebra.push_term!(context, coeff, permute_vector(exponent, sigma))
    end
    return AbstractAlgebra.finish(context)
  end

  delta = prod(
    (
      xi(i, l) - xi(i, k) for i in 1:n_vertices(Q) for k in 1:(d[i] - 1) for
      l in (k + 1):d[i]
    );
    init=R(1),
  )

  function antisymmetrize(f)
    out = R(0)
    for sigma in W
      out += __permutation_sign_product(sigma) * permute_polynomial(f, sigma)
    end
    return AbstractAlgebra.divexact(out, delta)
  end

  minimal_forbidden = all_destabilizing_subdimension_vectors(d, theta)
  filter!(
    e -> !any(partial_order(Q, f, e) for f in minimal_forbidden if f != e),
    minimal_forbidden,
  )
  verbose &&
    @info "there are $(length(minimal_forbidden)) minimal forbidden dimension vectors"

  function forbidden_polynomial(e::AbstractVector{Int})
    out = R(1)
    for i in 1:n_vertices(Q), j in 1:n_vertices(Q)
      for r in 1:e[i], s in (e[j] + 1):d[j]
        out *= (xi(j, s) - xi(i, r))^Q.adjacency[i, j]
      end
    end
    return out
  end
  forbidden = [forbidden_polynomial(e) for e in minimal_forbidden]

  chow_names = ["x$i$j" for i in 1:n_vertices(Q) for j in 1:d[i]]
  A, chow_vars = Oscar.graded_polynomial_ring(Oscar.QQ, chow_names, __chow_degrees(d))
  xs(i, j) = chow_vars[root_offset(i) + j]

  symmetric_generators = [symmetric_polynomial(k) for k in 1:maximum(d)]
  targets = typeof(R(0))[]
  for i in support(d)
    verbose && @info "computing the targets for vertex $(i) out of $(length(support(d)))"
    roots = [xi(i, j) for j in 1:d[i]]
    for k in 1:d[i]
      push!(targets, symmetric_generators[k](roots))
    end
  end
  verbose && @info "there are $(length(targets)) targets"

  inclusion = Oscar.hom(A, R, targets; check=true)
  verbose && @info "the inclusion map is built"

  antisymmetrized = typeof(R(0))[]
  verbose && @info "antisymmetrizing the forbidden polynomials, this may take a while..."
  for (index, polynomial) in enumerate(forbidden)
    verbose && @info "forbidden polynomial $(index) out of $(length(forbidden))"
    for basis_element in module_basis
      candidate = antisymmetrize(polynomial * basis_element)
      candidate != 0 && push!(antisymmetrized, candidate)
    end
  end
  verbose &&
    @info "there are $(length(antisymmetrized)) antisymmetrized forbidden polynomials"

  tautological_relations = typeof(A(0))[]
  for relation in antisymmetrized
    preimage_ideal = Oscar.preimage(inclusion, Oscar.ideal(R, [relation]))
    push!(tautological_relations, Oscar.gens(preimage_ideal)[1])
  end
  verbose && @info "there are $(length(tautological_relations)) tautological polynomials"

  linear_relation = sum((chi[i] * xs(i, 1) for i in support(d)); init=A(0))
  AQ = Oscar.quo(A, Oscar.ideal(A, vcat(tautological_relations, [linear_relation])))[1]
  return AQ
end

function __total_chern_class_universal(ring, d::AbstractVector{Int}, i::Int)
  vars = Oscar.gens(ring)
  offset = sum(d[1:(i - 1)]; init=0)
  return Oscar.simplify(ring(1) + sum((vars[offset + r] for r in 1:d[i]); init=ring(0)))
end

function __point_class(M::QuiverModuliSpace, X; unsafe::Bool=false)
  A = Oscar.chow_ring(X)
  N = __chow_dimension(M; unsafe=unsafe)
  total_chern_classes = [
    __total_chern_class_universal(A, M.d, i) for i in 1:n_vertices(M.Q)
  ]

  num = A(1)
  for i in 1:n_vertices(M.Q)
    c = total_chern_classes[i]
    exponent = sum(M.d[j] * M.Q.adjacency[j, i] for j in 1:n_vertices(M.Q))
    for _ in 1:exponent
      num = __truncate_by_degree(num * c, N)
    end
  end

  for i in 1:n_vertices(M.Q)
    inverse_chern = __truncate_by_degree(inv(total_chern_classes[i]), N)
    for _ in 1:M.d[i]
      num = __truncate_by_degree(num * inverse_chern, N)
    end
  end

  return Oscar.simplify(num[N])
end

function __copy_bundle(F)
  X = parent(F)
  if isdefined(F, :chern)
    return Oscar.abstract_bundle(X, Oscar.rank(F), Oscar.total_chern_class(F))
  end
  return Oscar.abstract_bundle(X, Oscar.chern_character(F))
end

function __attach_universal_weights!(M::QuiverModuliSpace, i::Int, F)
  return set_teleman_weights!(F, weights_universal_bundle(M, i; chi=linearization(M)))
end

function __attach_structure_sheaf_weights!(M::QuiverModuliSpace, F)
  weights = Dict(
    hn_type => [0] for hn_type in all_hn_types(M; unstable=true, ordered=false)
  )
  return set_teleman_weights!(F, weights)
end

function __attach_line_bundle_weights!(M::QuiverModuliSpace, eta::AbstractVector{Int}, F)
  return set_teleman_weights!(F, weights_line_bundle(M, eta))
end

function __attach_canonical_weights!(M::QuiverModuliSpace, F)
  return set_teleman_weights!(F, weights_canonical_bundle(M))
end

function __tautological_bundles(M::QuiverModuliSpace, X)
  bundles = Oscar.AbstractBundle{typeof(X)}[]
  for i in 1:n_vertices(M.Q)
    bundle = Oscar.abstract_bundle(
      X,
      M.d[i],
      __total_chern_class_universal(Oscar.chow_ring(X), M.d, i),
    )
    __attach_universal_weights!(M, i, bundle)
    push!(bundles, bundle)
  end
  return bundles
end

function __tangent_bundle(M::QuiverModuliSpace, bundles)
  X = parent(first(bundles))
  tangent = Oscar.trivial_line_bundle(X)

  for (i, j) in arrows(M.Q)
    tangent += Oscar.dual(bundles[i]) * bundles[j]
  end

  for i in 1:n_vertices(M.Q)
    tangent -= Oscar.dual(bundles[i]) * bundles[i]
  end

  return tangent
end

"""
    chow_ring(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}=canonical_stability(Q, d); chi::Union{AbstractVector{Int},UndefInitializer}=undef, verbose::Bool=false, unsafe::Bool=false)

Construct the quiver moduli space as an `Oscar.AbstractVariety` with its Chow ring,
point class, tangent bundle, and tautological bundles encoded in Oscar's
IntersectionTheory backend.
"""
function chow_ring(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int}=canonical_stability(Q, d);
  chi::Union{AbstractVector{Int},UndefInitializer}=undef,
  verbose::Bool=false,
  unsafe::Bool=false,
)
  return chow_ring(
    QuiverModuliSpace(Q, d, theta, "semistable");
    chi=chi,
    verbose=verbose,
    unsafe=unsafe,
  )
end

"""
    chow_ring(M::QuiverModuliSpace; chi::Union{AbstractVector{Int},UndefInitializer}=undef, verbose::Bool=false, unsafe::Bool=false)

Construct the quiver moduli space `M` as an `Oscar.AbstractVariety`.

The resulting abstract variety carries the tautological Chow presentation,
its point class, tautological bundles, and tangent bundle, so further
intersection-theory computations should use Oscar's API directly.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> X = chow_ring(M; chi=[-1, 1])
AbstractVariety of dim 6

julia> Oscar.point_class(X)
x23^2
```
"""
function chow_ring(
  M::QuiverModuliSpace;
  chi::Union{AbstractVector{Int},UndefInitializer}=undef,
  verbose::Bool=false,
  unsafe::Bool=false,
)
  __validate_chow_inputs(M.Q, M.d, M.theta; unsafe=unsafe, verbose=verbose)
  actual_chi = __resolve_linearization!(M, chi)

  if !isnothing(M.variety_cache[])
    return M.variety_cache[]
  end

  A = __presentation_ring(M.Q, M.d, M.theta; chi=actual_chi, verbose=verbose, unsafe=unsafe)
  X = Oscar.abstract_variety(__chow_dimension(M; unsafe=unsafe), A)

  bundles = __tautological_bundles(M, X)
  Oscar.set_point_class(X, __point_class(M, X; unsafe=unsafe))
  Oscar.set_tautological_bundles(X, bundles)
  Oscar.set_tangent_bundle(X, __tangent_bundle(M, bundles))

  M.variety_cache[] = X
  return X
end

Oscar.abstract_variety(
M::QuiverModuliSpace;
chi::Union{AbstractVector{Int},UndefInitializer}=undef,
verbose::Bool=false,
unsafe::Bool=false
) = chow_ring(M; chi=chi, verbose=verbose, unsafe=unsafe)

"""
    extended_gcd(x)

Compute the gcd and Bezout coefficients of a list of integers.
"""
function extended_gcd(x)
  n = length(x)
  if n == 1
    return [x[1], [1]]
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
