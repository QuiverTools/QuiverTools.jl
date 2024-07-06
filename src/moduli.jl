export QuiverModuli,
	is_nonempty,
	dimension,
	Poincaré_polynomial,
	Betti_numbers,
	is_smooth


struct QuiverModuli
	Q::Quiver
	d::Vector{Int}
	theta::Vector{Int}
	condition::String

end

function show(io::IO, M::QuiverModuli)
	print(
		io,
		"Moduli space of $(M.condition) representations of $(M.Q)
with dimension vector $(M.d) and stability parameter $(M.theta)",
	)
end

function is_nonempty(M::QuiverModuli)
	if M.condition == "stable"
		return has_stables(M.Q, M.d, M.theta)
	elseif M.condition == "semistable"
		return has_semistables(M.Q, M.d, M.theta)
	end
end

function dimension(M::QuiverModuli)
	if M.condition == "stable"
		return 1 - Euler_form(M.Q, M.d, M.d)
	elseif M.condition == "semistable"
		if has_semistables(M.Q, M.d, M.theta)
			return maximum(dimension_of_luna_stratum(tau)
						   for tau in all_luna_types(M.Q, M.d, M.theta))
		end
	else
		return "-∞" # how?
	end
	return "-∞"
end

"""
Checks if the moduli space is infinitesimally rigid by verifying wether
the Teleman quantization criterion of
[arXiv:2311.17003](https://doi.org/10.48550/arXiv.2311.17003) holds.

Only a sufficent criterion is implemented,
so the function may return "not known" even if the moduli space is rigid.
"""
function is_rigid(M::QuiverModuli)
	if is_acyclic(M.Q)
		bounds = all_Teleman_bounds(M.Q, M.d, M.theta)
		weights = all_weights_endomorphisms_universal_bundle(M.Q, M.d, M.theta)
		if all(weights[hn] < bounds[hn] for hn in keys(bounds))
			return true
		end
	end
	return "not known"
end


function dimension_of_luna_stratum(tau)
	throw(NotImplementedError())
end

function Poincaré_polynomial(M::QuiverModuli)
	throw(NotImplementedError())
end

function Betti_numbers(M::QuiverModuli)
	throw(NotImplementedError())
end

function is_smooth(M::QuiverModuli)
	throw(NotImplementedError())
end

# TODO - dimension of Luna stratum
# TODO - Poincaré polynomial
# TODO - Betti numbers
# TODO - is smooth 
# DONE - modify current Picard_rank, Hodge_diamond and so on
#        so that they behave correctly.
# TODO - Chow ring, point class, Todd class, Chern classes,
#        Chern characters, degrees, diagonal.
# TODO - smooth model should return a smooth model and
#        somehow the information for the correspondence?

"""Finds a linearization to construct the universal bundles on the moduli space."""
function linearization(M::QuiverModuli)
	if gcd(M.d) == 1
		throw(NotImplementedError())
	end
	throw(DomainError("Linearization does not exist as gcd(d) = $(gcd(d)) != 1"))
end


"""Returns the Chow ring for the given moduli space.
The optional argument `a` is the linearization used to construct the universal bundles.
If none is given, one is computed automatically."""
function Chow_ring(M::QuiverModuli, a::Vector{Int} = linearization(M), standard::Bool = false)
	return Chow_ring(M.Q, M.d, M.theta, a, standard)
end

"""Returns the point class for the given moduli space.
"""
function point_class(M::QuiverModuli)
	return point_class(M.Q, M.d, M.theta)
end

"""Returns the Todd class for the given moduli space.
"""
function Todd_class(M::QuiverModuli)
	return Todd_class(M.Q, M.d, M.theta)
end

function Hodge_polynomial(M::QuiverModuli)
	return Hodge_polynomial(M.Q, M.d, M.theta)
end

function Picard_rank(M::QuiverModuli)
	return Picard_rank(M.Q, M.d, M.theta)
end

function Hodge_diamond(M::QuiverModuli)
	return Hodge_diamond(M.Q, M.d, M.theta)
end

function index(M::QuiverModuli)
	# right? 
	return gcd(canonical_stability(M.Q, M.d))
end

struct QuiverModuliSpace <: QuiverModuli
	Q::Quiver
	d::AbstractVector{Int}
	theta::AbstractVector{Int}
	condition::String

	function QuiverModuliSpace(Q::Quiver,
        d::AbstractVector{Int},
        theta::AbstractVector{Int},
        condition::String)

        if condition in ["stable", "semistable"] &&
		   length(d) == nvertices(Q) &&
		   length(theta) == nvertices(Q)

			return new(Q, d, theta, condition)
		end
		throw(DomainError("Invalid input"))
	end

    function QuiverModuliSpace(Q::Quiver, d::AbstractVector{Int}, condition::String)
        return QuiverModuliSpace(Q, d, canonical_stability(Q, d), condition)
	end

    function QuiverModuliSpace(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
        return QuiverModuliSpace(Q, d, theta, "semistable")
    end

    function QuiverModuliSpace(Q::Quiver, d::AbstractVector{Int})
        return QuiverModuliSpace(Q, d, "semistable")
    end


end

struct QuiverModuliStack <: QuiverModuli
	Q::Quiver
	d::AbstractVector{Int}
	theta::AbstractVector{Int}
	condition::String

	function QuiverModuliStack(Q::Quiver,
        d::AbstractVector{Int},
        theta::AbstractVector{Int},
        condition::String)

        if condition in ["stable", "semistable"] &&
		   length(d) == nvertices(Q) &&
		   length(theta) == nvertices(Q)

			return new(Q, d, theta, condition)
		end
		throw(DomainError("Invalid input"))
	end

    function QuiverModuliStack(Q::Quiver, d::AbstractVector{Int}, condition::String)
        return QuiverModuliStack(Q, d, canonical_stability(Q, d), condition)
	end

    function QuiverModuliStack(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
        return QuiverModuliStack(Q, d, theta, "semistable")
    end

    function QuiverModuliStack(Q::Quiver, d::AbstractVector{Int})
        return QuiverModuliStack(Q, d, "semistable")
    end

end
