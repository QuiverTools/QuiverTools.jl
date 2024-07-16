
######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################


""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type."""
function Teleman_bound_onstratum(
	Q::Quiver,
	hntype::AbstractVector{AbstractVector{Int}},
	theta::AbstractVector{Int},
	denom::Function = sum,
)::Int
	if length(hntype) == 1
		throw(ArgumentError("Weight not defined for HN type of length 1."))
	end
	slopes = map(h -> slope(h, theta, denom), hntype)
	slopes = lcm(denominator.(slopes)) .* slopes
	return sum(
		(slopes[t] - slopes[s]) * Euler_form(Q, hntype[s], hntype[t]) for
		s ∈ 1:length(hntype)-1 for t ∈ s+1:length(hntype)
	)
end

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given ``Q``, ``d``, ``\\theta`` and `denom``."""
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


"""Computes the weight of the irreducible component of ``\\omega_R|_Z``
on a Harder-Narasimhan stratum for the 1-PS corresponding to each HN type.
More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum."""
function weight_irreducible_component_canonical_on_stratum(
	Q::Quiver,
	d::AbstractVector{Int},
	hntype::AbstractVector{AbstractVector{Int}},
	theta::AbstractVector{Int},
	denom::Function = sum,
)::Int
	kweights = map(di -> slope(di, theta, denom), hntype)
	kweights = kweights * lcm(denominator.(kweights))

	dd = sum(kweights[m] .* hntype[m] for m ∈ 1:length(hntype))
	# The Fano paper shows that under appropriate conditions,
	# the canonical bundle is given by linearizing with minus
	# the canonical stability parameter.
	can = canonical_stability(Q, d)
	can /= gcd(can)
	return can' * dd
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
	return [kweights[i] - kweights[j] for i ∈ 1:length(hntype) for j ∈ 1:length(hntype)]
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

"""
    does_Teleman_inequality_hold(M::QuiverModuli)

Checks if the Teleman quantization criterion of
[arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003) holds.

In case the quiver is acyclic, this ensures that the moduli space is infinitesimally
rigid.
"""
function does_Teleman_inequality_hold(M::QuiverModuli)
    bounds = all_Teleman_bounds(M.Q, M.d, M.theta)
    weights = all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta)
    return all(weights[hn] < bounds[hn] for hn in collect(keys(bounds)))
end