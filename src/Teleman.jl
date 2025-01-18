######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################

export teleman_bound_on_stratum,
  all_teleman_bounds,
  weights_universal_bundle_on_stratum,
  all_weights_universal_bundle,
  weight_irreducible_component_canonical_on_stratum,
  all_weights_irreducible_component_canonical,
  weights_endomorphism_universal_bundle_on_stratum,
  all_weights_endomorphisms_universal_bundle,
  does_rigidity_inequality_hold

# TODO is there a reasoning behind hntype vs all_hn_types? why not hn_type?
"""
    teleman_bound_on_stratum(Q::Quiver, hn_type, theta, denom=sum)

Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type.

# Input

- `Q`: a quiver
- `hn_type`: a Harder-Narasimhan type
- `theta`: a stability parameter.
- `denom`: a denominator for the slope function. Defaults to `sum`.

"""
function teleman_bound_on_stratum(
  Q::Quiver,
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)::Int
  length(hn_type) == 1 &&
    throw(ArgumentError("Weight not defined on the dense stratum"))

  slopes = map(h -> slope(h, theta, denom), hn_type)
  slopes = lcm(denominator.(slopes)) .* slopes
  return sum(
    (slopes[t] - slopes[s]) * euler_form(Q, hn_type[s], hn_type[t]) for
    s in 1:(length(hn_type) - 1) for t in (s + 1):length(hn_type)
  )
end

function teleman_bound_on_stratum(M::QuiverModuli, hn_type::HNType)
  return teleman_bound_on_stratum(M.Q, hn_type, M.theta, M.denom)
end

"""
    all_teleman_bounds(Q::Quiver, d, theta, denom=sum)

Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given ``Q``, ``d``, ``\\theta`` and `denom``.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> all_teleman_bounds(Q, [2, 3], [3, -2])
Dict{HNType{2}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90
```
"""
function all_teleman_bounds(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  #This is only relevant on the unstable locus
  HN = filter(hn_type -> hn_type != [d], all_hn_types(Q, d, theta, denom))
  return Dict(
    hn_type => teleman_bound_on_stratum(Q, hn_type, theta, denom) for hn_type in HN
  )
end

"""
	all_teleman_bounds(M::QuiverModuli)

# Examples

```jldoctest
julia> Q = three_vertex_quiver(1, 2, 3); d = [3, 1, 2]; theta = [5, 3, -9];

julia> M = QuiverModuliSpace(Q, d, theta);

julia> all_teleman_bounds(M)
Dict{HNType{3}, Int64} with 24 entries:
  [[2, 1, 1], [1, 0, 1]]                       => 12
  [[1, 0, 0], [0, 1, 0], [2, 0, 1], [0, 0, 1]] => 306
  [[1, 0, 0], [1, 1, 0], [1, 0, 1], [0, 0, 1]] => 131
  [[2, 0, 0], [1, 0, 1], [0, 1, 1]]            => 64
  [[3, 0, 0], [0, 1, 2]]                       => 150
  [[1, 1, 0], [2, 0, 1], [0, 0, 1]]            => 312
  [[2, 0, 0], [1, 1, 1], [0, 0, 1]]            => 336
  [[2, 0, 0], [1, 1, 0], [0, 0, 2]]            => 242
  [[3, 0, 0], [0, 1, 1], [0, 0, 1]]            => 168
  [[3, 1, 1], [0, 0, 1]]                       => 432
  [[3, 0, 0], [0, 1, 0], [0, 0, 2]]            => 246
  [[0, 1, 0], [3, 0, 2]]                       => 108
  [[0, 1, 0], [2, 0, 1], [1, 0, 1]]            => 76
  [[1, 0, 0], [2, 0, 1], [0, 1, 1]]            => 122
  [[1, 0, 0], [2, 1, 1], [0, 0, 1]]            => 92
  [[2, 0, 0], [0, 1, 0], [1, 0, 2]]            => 312
  [[1, 0, 0], [2, 1, 2]]                       => 18
  [[2, 0, 0], [0, 1, 0], [1, 0, 1], [0, 0, 1]] => 132
  [[1, 0, 0], [1, 1, 1], [1, 0, 1]]            => 68
  ⋮                                            => ⋮
```
"""
function all_teleman_bounds(M::QuiverModuli)
  return all_teleman_bounds(M.Q, M.d, M.theta, M.denom)
end

"""
    weights_universal_bundle_on_stratum(theta, hn_type, denom=sum; chi)

Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type.

"""
function weights_universal_bundle_on_stratum(
  theta::AbstractVector{Int},
  hn_type,
  denom::Function=sum;
  chi::AbstractVector{Int},
)::AbstractVector{Int}
  slopes = map(h -> slope(h, theta, denom), hn_type)
  slopes *= lcm(denominator.(slopes))

  constant_term = sum(slopes[i] * (chi' * hn_type[i]) for i in eachindex(hn_type))

  return -constant_term .+ slopes
end

function weights_universal_bundle_on_stratum(
  M::QuiverModuli,
  hn_type::HNType;
  chi::AbstractVector{Int}=extended_gcd(M.d)[2],
)
  return weights_universal_bundle_on_stratum(M.theta, hn_type, M.denom; chi=chi)
end

"""
    all_weights_universal_bundle(Q::Quiver, d, theta, denom=sum; chi)

Computes the Teleman weights of the universal bundle ``U_i(chi)``
for the linearization ``chi`` on all the non-dense Harder-Narasimhan strata.

"""
function all_weights_universal_bundle(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum;
  chi::AbstractVector{Int}=extended_gcd(d)[2],
)
  HN = filter(hn_type -> hn_type[1] != d, all_hn_types(Q, d, theta, denom))
  return Dict(
    hn_type => weights_universal_bundle_on_stratum(theta, hn_type, denom; chi=chi) for
    hn_type in HN
  )
end

"""
    all_weights_universal_bundle(M::QuiverModuli; chi)

Computes the Teleman weights of the universal bundle ``U_i(chi)``
for the linearization ``chi`` on all the non-dense Harder-Narasimhan strata.

"""
function all_weights_universal_bundle(
  M::QuiverModuli;
  chi::AbstractVector{Int}=extended_gcd(M.d)[2],
)
  return all_weights_universal_bundle(M.Q, M.d, M.theta, M.denom; chi=chi)
end

# TODO test and add safety checks.

"""
    weight_irreducible_component_canonical_on_stratum(Q::Quiver, d, hn_type, theta, denom=sum)

Computes the Teleman weight of the irreducible component of ``\\omega_R|_Z``
on the Harder-Narasimhan stratum `hn_type`.
More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum.

"""
function weight_irreducible_component_canonical_on_stratum(
  Q::Quiver,
  d::AbstractVector{Int},
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)::Int
  kweights = map(di -> slope(di, theta, denom), hn_type)
  kweights = kweights * lcm(denominator.(kweights))

  dd = sum(kweights[m] .* hn_type[m] for m in 1:length(hn_type))
  # The Fano paper shows that under appropriate conditions,
  # the canonical bundle is given by linearizing with minus
  # the canonical stability parameter.
  can = canonical_stability(Q, d)
  can /= gcd(can)
  return can' * dd
end

# TODO test and add safety checks.

"""
    weight_irreducible_component_canonical_on_stratum(M::QuiverModuli, hn_type)

Computes the Teleman weight of the irreducible component of ``\\omega_R|_Z``
on the Harder-Narasimhan stratum `hn_type`.
More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum.

"""
function weight_irreducible_component_canonical_on_stratum(
  M::QuiverModuli,
  hn_type::HNType,
)
  return weight_irreducible_component_canonical_on_stratum(
    M.Q,
    M.d,
    hn_type,
    M.theta,
    M.denom,
  )
end

"""
    all_weights_irreducible_component_canonical(Q::Quiver, d, theta, denom=sum)

Computes the Teleman weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of
``\\mathcal{O}(H)`` on each stratum.

"""
function all_weights_irreducible_component_canonical(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  HN = filter(hn_type -> hn_type != [d], all_hn_types(Q, d, theta))
  return Dict(
    hn_type =>
      weight_irreducible_component_canonical_on_stratum(Q, d, hn_type, theta, denom)
    for hn_type in HN
  )
end

"""
    all_weights_irreducible_component_canonical(M::QuiverModuli)

Computes the Teleman weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of
``\\mathcal{O}(H)`` on each stratum.

"""
function all_weights_irreducible_component_canonical(M::QuiverModuli)
  return all_weights_irreducible_component_canonical(M.Q, M.d, M.theta, M.denom)
end

"""
    weights_endomorphism_universal_bundle_on_stratum(hn_type, theta, denom=sum)

Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type.

"""
function weights_endomorphism_universal_bundle_on_stratum(
  hn_type::HNType,
  theta::AbstractVector{Int},
  denom::Function=sum,
)::AbstractVector{Int}

  # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is
  # slope of first term in the HN type - slope of the last term in the HN type
  kweights = map(di -> slope(di, theta, denom), hn_type)
  kweights = kweights * lcm(denominator.(kweights))
  # return kweights[1] - kweights[end] # this is the largest one
  return [kweights[i] - kweights[j] for i in 1:length(hn_type) for j in 1:length(hn_type)]
end

"""
    weights_endomorphism_universal_bundle_on_stratum(M::QuiverModuli, hn_type)

Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type.

"""
function weights_endomorphism_universal_bundle_on_stratum(
  M::QuiverModuli,
  hn_type::HNType,
)
  return weights_endomorphism_universal_bundle_on_stratum(hn_type, M.theta, M.denom)
end

"""
    all_weights_endomorphisms_universal_bundle(Q::Quiver, d, theta, denom=sum)

Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.

"""
function all_weights_endomorphisms_universal_bundle(
  Q::Quiver,
  d::AbstractVector{Int},
  theta::AbstractVector{Int},
  denom::Function=sum,
)
  HN = filter(hn_type -> hn_type != [d], all_hn_types(Q, d, theta, denom))
  return Dict(
    [hn_type, weights_endomorphism_universal_bundle_on_stratum(hn_type, theta, denom)] for
    hn_type in HN
  )
end

"""
    all_weights_endomorphisms_universal_bundle(Q::Quiver, d, theta, denom=sum)

Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.

"""
function all_weights_endomorphisms_universal_bundle(M::QuiverModuli)
  return all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta, M.denom)
end

"""
    does_rigidity_inequality_hold(M::QuiverModuli)

Checks if the Teleman quantization criterion of
[arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003) holds.

In case the quiver is acyclic, this ensures that the moduli space is infinitesimally
rigid.

# Examples

Our favourite 6-fold is rigid:
```jldoctest
julia> Q = kronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> does_rigidity_inequality_hold(M)
true
```

Some quiver moduli are rigid, but it can't be proved by this criterion:
the following moduli space can be shown to be ``\\mathbb{P}^6``, whose rigidity follows
from the Euler sequence (see
[Example 4.8, arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003)). However,
the Teleman inequality does not hold:
```jldoctest
julia> Q = three_vertex_quiver(1, 6, 1); M = QuiverModuliSpace(Q, [1, 6, 6], [42, 5, -12]);

julia> does_rigidity_inequality_hold(M)
false
```
"""
function does_rigidity_inequality_hold(M::QuiverModuli)
  bounds = all_teleman_bounds(M.Q, M.d, M.theta)
  weights = all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta)
  return all(maximum(weights[hn]) < bounds[hn] for hn in collect(keys(bounds)))
end
