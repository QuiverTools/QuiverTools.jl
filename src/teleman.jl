
######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################

export Teleman_bound_onstratum, all_Teleman_bounds, weights_universal_bundle_onstratum,
	all_weights_universal_bundle, weight_irreducible_component_canonical_on_stratum,
	all_weights_irreducible_component_canonical,
	weights_endomorphism_universal_bundle_on_stratum,
	all_weights_endomorphisms_universal_bundle,
	does_Teleman_inequality_hold

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type."""
function Teleman_bound_onstratum(
	Q::Quiver,
	hntype::Vector{<:AbstractVector{Int}},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)::Int

	if length(hntype) == 1
		throw(ArgumentError("Weight not defined for HN type of length 1."))
	end
	slopes = map(h -> slope(h, theta, denom), hntype)
	slopes = lcm(denominator.(slopes)) .* slopes
	return sum(
			(slopes[t] - slopes[s]) * Euler_form(Q, hntype[s], hntype[t])
			for s in 1:length(hntype)-1
			for t in s+1:length(hntype)
		)
end

function Teleman_bound_onstratum(M::QuiverModuli,
	hntype::Vector{<:AbstractVector{Int}},
	)
	return Teleman_bound_onstratum(M.Q, hntype, M.theta, M.denom)
end


"""
Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given ``Q``, ``d``, ``\\theta`` and `denom``.
	
	
EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> all_Teleman_bounds(Q, [2, 3], [3, -2])
Dict{Vector{AbstractVector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90
```
"""
function all_Teleman_bounds(
	Q::Quiver,
	d::AbstractVector{Int},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)

	#This is only relevant on the unstable locus
	HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
	return Dict([hntype, Teleman_bound_onstratum(Q, hntype, theta, denom)] for hntype in HN)
end
"""

Interface for `all_Teleman_bounds(Q, d, theta)`.

EXAMPLE:
```jldoctest
julia> Q = three_vertex_quiver(1, 2, 3); d = [3, 1, 2]; theta = [5, 3, -9];

julia> all_Teleman_bounds(Q, d, theta)
Dict{Vector{AbstractVector{Int64}}, Int64} with 24 entries:
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
function all_Teleman_bounds(M::QuiverModuli)

	return all_Teleman_bounds(M.Q, M.d, M.theta, M.denom)
end

"""Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type."""
function weights_universal_bundle_onstratum(
	theta::AbstractVector{Int},
	a::AbstractVector{Int},
	hntype,
	denom::Function = sum,
	)::AbstractVector{Int}

	slopes = map(h -> slope(h, theta, denom), hntype)
	slopes *= lcm(denominator.(slopes))

	constant_term = sum(slopes[i] * (a' * hntype[i]) for i in eachindex(hntype))

	return -constant_term .+ slopes
end

function weights_universal_bundle_onstratum(M::QuiverModuli,
	hntype::Vector{<:AbstractVector{Int}},
	a::AbstractVector{Int} = extended_gcd(M.d)[2])

	return weights_universal_bundle_onstratum(M.theta, a, hntype, M.denom)
end

"""Computes the weights of the universal bundle ``U_i(a)`` for the linearization ``a``
on all the non-dense Harder-Narasimhan strata for each 1-PS
corresponding to each HN type."""
function all_weights_universal_bundle(
	Q::Quiver,
	d::AbstractVector{Int},
	theta::AbstractVector{Int},
	a::AbstractVector{Int},
	denom::Function = sum,
	)
	
	HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
	return Dict(
		[hntype, weights_universal_bundle_onstratum(theta, a, hntype, denom)] for
		hntype in HN
	)
end

function all_weights_universal_bundle(M::QuiverModuli,
	a::AbstractVector{Int} = extended_gcd(M.d)[2])

	return all_weights_universal_bundle(M.Q, M.d, M.theta, a, M.denom)
end


"""Computes the weight of the irreducible component of ``\\omega_R|_Z``
on a Harder-Narasimhan stratum for the 1-PS corresponding to each HN type.
More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum."""
function weight_irreducible_component_canonical_on_stratum(
	Q::Quiver,
	d::AbstractVector{Int},
	hntype::Vector{<:AbstractVector{Int}},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)::Int

	kweights = map(di -> slope(di, theta, denom), hntype)
	kweights = kweights * lcm(denominator.(kweights))

	dd = sum(kweights[m] .* hntype[m] for m in 1:length(hntype))
	# The Fano paper shows that under appropriate conditions,
	# the canonical bundle is given by linearizing with minus
	# the canonical stability parameter.
	can = canonical_stability(Q, d)
	can /= gcd(can)
	return can' * dd
end

function weight_irreducible_component_canonical_on_stratum(M::QuiverModuli,
	hntype::Vector{<:AbstractVector{Int}})

	return weight_irreducible_component_canonical_on_stratum(M.Q,
															M.d,
															hntype,
															M.theta,
															M.denom)
end


"""Computes the weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of
``\\mathcal{O}(H)`` on each stratum."""
function all_weights_irreducible_component_canonical(
	Q::Quiver,
	d::AbstractVector{Int},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)

	HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta))
	return Dict(
		[
			hntype,
			weight_irreducible_component_canonical_on_stratum(Q, d, hntype, theta, denom),
		] for hntype in HN
	)
end

function all_weights_irreducible_component_canonical(M::QuiverModuli)

	return all_weights_irreducible_component_canonical(M.Q, M.d, M.theta, M.denom)
end


"""Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type."""
function weights_endomorphism_universal_bundle_on_stratum(
	hntype::AbstractVector{AbstractVector{Int}},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)::AbstractVector{Int}

	# the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is
    # slope of first term in the HN type - slope of the last term in the HN type
	kweights = map(di -> slope(di, theta, denom), hntype)
	kweights = kweights * lcm(denominator.(kweights))
	# return kweights[1] - kweights[end] # this is the largest one
	return [kweights[i] - kweights[j] for i in 1:length(hntype) for j in 1:length(hntype)]
end

function weights_endomorphism_universal_bundle_on_stratum(M::QuiverModuli,
	hntype::Vector{<:AbstractVector{Int}}
	)

	return weights_endomorphism_universal_bundle_on_stratum(hntype, M.theta, M.denom)
end


"""Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type."""
function all_weights_endomorphisms_universal_bundle(
	Q::Quiver,
	d::AbstractVector{Int},
	theta::AbstractVector{Int},
	denom::Function = sum,
	)

	HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
	return Dict(
		[hntype, weights_endomorphism_universal_bundle_on_stratum(hntype, theta, denom)] for
		hntype in HN
	)
end

function all_weights_endomorphisms_universal_bundle(M::QuiverModuli)

	return all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta, M.denom)
end


"""
    does_Teleman_inequality_hold(M::QuiverModuli)

Checks if the Teleman quantization criterion of
[arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003) holds.

In case the quiver is acyclic, this ensures that the moduli space is infinitesimally
rigid.

EXAMPLES:

Our favourite 6-fold is rigid:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> does_Teleman_inequality_hold(M)
true
```

Some quiver moduli are rigid, but it can't be proved by this criterion:
the following moduli space can be shown to be ``\\mathbb{P}^6``, whose rigidity follows 
from the Euler sequence (see
[Example 4.8, arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003)). However,
the Teleman inequality does not hold:
```jldoctest
julia> Q = three_vertex_quiver(1, 6, 1); M = QuiverModuliSpace(Q, [1, 6, 6], [42, 5, -12]);

julia> does_Teleman_inequality_hold(M)
false
```
"""
function does_Teleman_inequality_hold(M::QuiverModuli)
    bounds = all_Teleman_bounds(M.Q, M.d, M.theta)
    weights = all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta)
    return all(maximum(weights[hn]) < bounds[hn] for hn in collect(keys(bounds)))
end