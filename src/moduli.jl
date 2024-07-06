export QuiverModuli, QuiverModuliSpace, QuiverModuliStack

abstract type QuiverModuli end

struct QuiverModuliSpace <: QuiverModuli
    Q::Quiver
    d::AbstractVector{Int}
    theta::AbstractVector{Int}
    condition::String
    denom::Function

    function QuiverModuliSpace(
        Q::Quiver,
        d::AbstractVector{Int},
        theta::AbstractVector{Int} = canonical_stability(Q, d),
        condition::String = "semistable",
        denom::Function = sum,
    )

        if condition in ["stable", "semistable"] &&
           length(d) == nvertices(Q) &&
           length(theta) == nvertices(Q)

            return new(Q, d, theta, condition, denom)
        end
        throw(DomainError("Invalid input"))
    end
end

struct QuiverModuliStack <: QuiverModuli
    Q::Quiver
    d::AbstractVector{Int}
    theta::AbstractVector{Int}
    condition::String

    function QuiverModuliStack(
        Q::Quiver,
        d::AbstractVector{Int},
        theta::AbstractVector{Int} = canonical_stability(Q, d),
        condition::String = "semistable",
        denom::Function = sum,
    )

        if condition in ["stable", "semistable"] &&
           length(d) == nvertices(Q) &&
           length(theta) == nvertices(Q)

            return new(Q, d, theta, condition, denom)
        end
        throw(DomainError("Invalid input"))
    end
end

function show(io::IO, M::QuiverModuli)
    print(
        io,
        "Moduli space of $(M.condition) representations of $(M.Q)
      with dimension vector $(M.d) and stability parameter $(M.theta)",
    )
end


"""
Checks if the moduli space or stack is nonempty.


Examples
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> is_nonempty(M)
true

julia> M = QuiverModuliSpace(Q, [2, 3], [-3, 2]);

julia> is_nonempty(M)
false
```
"""
function is_nonempty(M::QuiverModuli)
    if M.condition == "stable"
        return has_stables(M.Q, M.d, M.theta)
    elseif M.condition == "semistable"
        return has_semistables(M.Q, M.d, M.theta)
    end
end

function is_theta_coprime(M::QuiverModuli)
    return is_coprime(M.d, M.theta)
end


function all_HN_types(M::QuiverModuli; proper::Bool = false, ordered::Bool = true)
    HN = all_HN_types(M.Q, M.d, M.theta, M.denom, ordered = ordered)
    if proper
        return filter(hn_type -> hn_type != [M.d], HN)
    end
end

function is_HN_type(M::QuiverModuli, hn_type::AbstractVector{AbstractVector{Int}})::Bool
    return is_HN_type(M.Q, M.d, M.theta, hn_type, M.denom)
end

function codimension_HN_stratum(
    M::QuiverModuli,
    hn_type::AbstractVector{AbstractVector{Int}},
)
    return codimension_HN_stratum(M.Q, M.d, M.theta, hn_type, M.denom)
end

function codimension_unstable_locus(M::QuiverModuli)
    HN = all_HN_types(M, proper = true)
    return minimum(codimension_HN_stratum(M, hn_type) for hn_type in HN)
end

# TODO add safety checks everywhere in the codebase


function all_Luna_types(M::QuiverModuli; exclude_stable::Bool = false)
    return all_Luna_types(M.Q, M.d, M.theta, M.denom, exclude_stable)
end

function all_Luna_types(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum,
    exclude_stable::Bool = false,
)

    # subdimensions with the same slope as d
    same_slope = filter(
        e ->
            slope(e, theta, denom) == slope(d, theta, denom) &&
                has_stables(Q, e, theta, denom),
        QuiverTools.all_subdimension_vectors(d, nonzero = true),
    )

    Luna_types = []

    bound = sum(d) ÷ minimum(sum(e) for e in same_slope) # the highest possible amount of repetitions for a given stable dimension vector

    for i = 1:bound+1
        for tau in with_replacement_combinations(same_slope, i)
            if sum(tau) == d

                partial = Dict()
                for e in tau
                    if e in keys(partial)
                        partial[e] += 1
                    else
                        partial[e] = 1
                    end
                end

                for e in keys(partial)
                    partial[e] = partitions(partial[e])
                end

                new_Luna = [
                    Dict(zip(keys(partial), p)) for
                    p in IterTools.product(values(partial)...)
                ]
                vcat!(Luna_types, new_Luna)
            end
        end
    end

    if exclude_stable
        return filter(luna -> luna != Dict(d => [1]), Luna_types)
    end
    return Luna_types
end

function is_Luna_type(M::QuiverModuli, tau)
    if sum(M.d) == 0
        return tau == Dict(d => [1])
    end

    if sum(keys(tau)) != M.d
        return false
    end
    if !all(slope(e, M.theta, M.denom) == slope(d, M.theta, M.denom) for e in keys(tau))
        return false
    end

    if !all(has_semistables(M.Q, e, M.theta, M.denom) for e in keys(tau))
        return false
    end
    return true
end

function dimension_of_Luna_stratum(M::QuiverModuli, tau)
    return sum(length(tau[e]) * (1 - Euler_form(M.Q, e, e)) for e in keys(tau))
end


function local_quiver_setting(M::QuiverModuli, tau)
    if !is_Luna_type(M, tau)
        throw(DomainError("Not a Luna type"))
    end

    A = coerce_matrix([
        [generic_ext(M.Q, e, eprime) for eprime in keys(tau) for n in tau[eprime]] for
        e in keys(tau) for m in tau[e]
    ])

    Qloc = Quiver(A)
    dloc = [m for e in keys(tau) for m in tau[e]]

    return Dict("Q" => Qloc, "d" => dloc)
end

function semistable_equals_stable(M::QuiverModuli)

    if is_coprime(M.d, M.theta) || !has_semistables(M.Q, M.d, M.theta, M.denom)
        return true
    end
    return length(all_Luna_types(M, exclude_stable = true)) == 0
end


# TODO what if dim R = 1? or 0?
"""
	is_amply_stable(M::QuiverModuli)

Checks wether the dimension vector ``d`` is amply stable
with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus
in the parameter space is at least ``2``.

Arguments:
- `M::QuiverModuli`: a moduli space of representations of a quiver.

Returns:

- `true` if the moduli space is amply stable, false otherwise.

"""
function is_amply_stable(M::QuiverModuli)
    return codimension_unstable_locus(M) >= 2
end




"""
Solve ``A\\cdot x = b`` for ``A`` upper triangular via back substitution
"""
function solve(A, b)
    n = length(b)
    x = Vector{Any}(zeros(n))

    x[n] = b[n] / A[n, n]

    for i ∈ n-1:-1:1
        x[i] = (b[i] - sum(A[i, j] * x[j] for j ∈ i+1:n)) / A[i, i]
    end
    return x
end

"""
Cardinality of general linear group \$\\mathrm{GL}_n(\\mathbb{F}_v)\$.
"""
@memoize Dict function CardinalGl(n::Int, q)
    if n == 0
        return 1
    else
        return prod(q^n - q^i for i ∈ 0:n-1)
    end
end

"""
Cardinality of representation space \$\\mathrm{R}(Q,d)\$, over \$\\mathbb{F}_q\$.
"""
function CardinalRd(Q::Quiver, d::AbstractVector{Int}, q)
    return q^sum(d[i] * d[j] * Q.adjacency[i, j] for i ∈ 1:nvertices(Q), j ∈ 1:nvertices(Q))
end

"""
Cardinality of product of general linear groups \$\\mathrm{GL}_{d}(\\mathbb{F}_q)\$.
"""
@memoize Dict function CardinalGd(d::AbstractVector{Int}, q)
    return prod(CardinalGl(di, q) for di in d)
end


# TODO document this
"""Entry of the transfer matrix, as per Corollary 6.9"""
function TransferMatrixEntry(Q, e, f, q)
    fe = f - e

    if all(fei >= 0 for fei in fe)
        return q^Euler_form(Q, -fe, e) * CardinalRd(Q, fe, q) / CardinalGd(fe, q)
    else
        return 0
    end
end

function Td(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, q)
    # indexing set for the transfer matrix
    I = filter(
        e -> slope(e, theta) > slope(d, theta),
        all_subdimension_vectors(d, nonzero = true, strict = true),
    )
    I = vcat([zero_vector(nvertices(Q))], I, [d])

    l = length(I)
    T = Matrix{Any}(zeros(l, l))

    for (i, Ii) in enumerate(I)
        for j ∈ i:l  # upper triangular
            T[i, j] = TransferMatrixEntry(Q, Ii, I[j], q)
        end
    end
    return T
end


# TODO DOI below is not open access.
# auxiliary functions for Hodge_polynomial() above
###################################################

"""
Returns the Hodge polynomial of the moduli space of ``\\theta``-semistable
representations of ``Q`` with dimension vector ``d``.

The algorithm is based on [MR1974891](https://doi.org/10.1007/s00222-002-0273-4),
and the current implementation is translated from the [Hodge diamond cutter]
(https://zenodo.org/doi/10.5281/zenodo.3893509).

Examples:

```jldoctest
julia> Q = mKronecker_quiver(3);

julia> d = [2, 3];

julia> theta = [3, -2];

julia> Hodge_polynomial(Q, d, theta)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1
```
"""
function Hodge_polynomial(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})

    # safety checks
    if theta' * d == 0 && !is_coprime(d)
        throw(ArgumentError("d is not coprime"))
    elseif !is_acyclic(Q)
        throw(ArgumentError("Q is not acyclic."))
    end

    R, q = polynomial_ring(Singular.QQ, ["q"])
    F = fraction_field(R)

    v = F(q[1]) # worsens performance by ~8%. Necessary?

    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])

    # @warn "result needs to be a polynomial, otherwise the moduli space is singular."
    solution = solve(T, one_at_the_end)[1] * (1 - v)
    if denominator(solution) != 1
        throw(DomainError("Moduli space is singular!"))
    end
    result = numerator(solution)

    S, (x, y) = polynomial_ring(Singular.QQ, ["x", "y"])
    return result(x * y)
end

function Hodge_polynomial(M::QuiverModuli)
    return Hodge_polynomial(M.Q, M.d, M.theta)
end

"""
Returns the Hodge diamond of the moduli space of
``\\theta``-semistable representations of ``Q`` with dimension vector ``d``.
"""
function Hodge_diamond(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
    g = Hodge_polynomial(Q, d, theta)

    return map(
        ind -> coeff(g, [1, 2], [ind[1] - 1, ind[2] - 1]),
        Iterators.product(1:degree(g, 1)+1, 1:degree(g, 2)+1),
    )
end

function Hodge_diamond(M::QuiverModuli)
    return Hodge_diamond(M.Q, M.d, M.theta)
end

"""
Computes the Picard rank of the moduli space of
``\\theta``-semistable representations of ``Q`` with dimension vector ``d``.
"""
function Picard_rank(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
    # TODO If over the complex numbers this should be h^{1,1},
    # since the moduli space is rational.
    # TODO This should follow from the long exact sequence in cohomology
    # given by the exponential short exact sequence.

    return coeff(Hodge_polynomial(Q, d, theta), 2).num
end

function Picard_rank(M::QuiverModuli)
    return Picard_rank(M.Q, M.d, M.theta)
end

function _Hodge_polynomial_fast(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
)
    # unsafe, curate input!
    # this is about 2% faster than the above, and occupies about 2% less memory.

    R, q = polynomial_ring(AbstractAlgebra.QQ, ["q"])
    F = fraction_field(R)
    v = F(q[1]) # worsens performance by ~8%. Necessary?

    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])

    result = numerator(solve(T, one_at_the_end)[1] * (1 - v))
    return [coeff(result, i) for i = 0:degree(result)] # this is actually all
    # we need for the Hodge diamond because the matrix is diagonal for quiver moduli
end






###############################################################################
# tautological representation of the Chow ring.
# Implements the results of [arXiv:1307.3066](https://doi.org/10.48550/arXiv.1307.3066) and
# [arXiv.2307.01711](https://doi.org/10.48550/arXiv.2307.01711).
###############################################################################

# partial order on the forbidden dimension vectors as in https://doi.org/10.48550/arXiv.1307.3066
function partial_order(Q::Quiver, f::AbstractVector{Int}, g::AbstractVector{Int})
    if !all(f[i] <= g[i] for i ∈ 1:nvertices(Q) if is_source(Q, i))
        return false
    elseif !all(f[i] >= g[i] for i ∈ 1:nvertices(Q) if is_sink(Q, i))
        return false
    elseif !all(f[i] == g[i] for i ∈ 1:nvertices(Q) if !is_source(Q, i) && !is_sink(Q, i))
        return false
    end
    return true
end


"""
Returns the symmetric polynomial of degree ``degree`` in the variables ``vars``.
"""
function symmetric_polynomial(vars, degree::Int)
    return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

function Chow_ring(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    a::AbstractVector{Int},
)
    # TODO cover case d[i] = 0
    # safety checks
    if !is_coprime(d, theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif a' * d != 1
        throw(ArgumentError("a is not a linearization"))
    end

    varnames = ["x$i$j" for i ∈ 1:nvertices(Q) for j ∈ 1:d[i] if d[i] > 0]
    # R, vars = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, varnames)
    R, vars = Singular.polynomial_ring(Singular.QQ, varnames)
    function chi(i, j)
        return vars[sum(d[1:i-1])+j]
    end

    function base_for_ring(name = "naive")
        if name == "naive"
            bounds = [0:(d[i]-nu) for i ∈ 1:nvertices(Q) for nu ∈ 1:d[i]]
            lambdas = Iterators.product(bounds...)

            build_elem(lambda) = prod(
                prod(chi(i, nu)^lambda[sum(d[1:i-1])+nu] for nu ∈ 1:d[i]) for
                i ∈ 1:nvertices(Q)
            )

            return map(l -> build_elem(l), lambdas)
        else
            throw(ArgumentError("unknown base."))
        end
    end

    # build the permutation group W
    W = Iterators.product([AbstractAlgebra.SymmetricGroup(d[i]) for i ∈ 1:nvertices(Q)]...)
    sign(w) = prod(AbstractAlgebra.sign(wi) for wi in w)

    permute(f, sigma) = f([chi(i, sigma[i][j]) for i ∈ 1:nvertices(Q) for j ∈ 1:d[i]]...)

    delta = prod(
        prod(chi(i, l) - chi(i, k) for k ∈ 1:d[i]-1 for l ∈ k+1:d[i]) for
        i ∈ 1:nvertices(Q) if d[i] > 1
    )
    antisymmetrize(f) = sum(sign(w) * permute(f, w) for w in W) / delta

    function all_forbidden(Q, d, theta, denom::Function = sum)
        dest = all_destabilizing_subdimension_vectors(d, theta, denom)
        return filter(
            e -> !any(f -> partial_order(Q, f, e), filter(f -> f != e, dest)),
            dest,
        )
    end

    forbidden_polynomials = [
        prod(
            prod((chi(j, s) - chi(i, r))^Q.adjacency[i, j] for r ∈ 1:e[i], s ∈ e[j]+1:d[j]) for j ∈ 1:nvertices(Q), i ∈ 1:nvertices(Q) if
            Q.adjacency[i, j] > 0 && e[i] > 0 && d[j] > 1
        ) for e in all_forbidden(Q, d, theta)
    ]

    varnames2 = ["x$i$j" for i ∈ 1:nvertices(Q) for j ∈ 1:d[i] if d[i] > 0]
    A, Avars = Singular.polynomial_ring(Singular.QQ, varnames2)

    function xs(i, j)
        return Avars[sum(d[1:i-1])+j]
    end

    targets = [
        [symmetric_polynomial([chi(i, j) for j ∈ 1:d[i]], k) for k ∈ 1:d[i]] for
        i ∈ 1:nvertices(Q)
    ]
    targets = reduce(vcat, targets)

    inclusion = AlgebraHomomorphism(A, R, targets)

    anti = [antisymmetrize(f * b) for f in forbidden_polynomials for b in base_for_ring()]
    tautological = [gens(preimage(inclusion, Ideal(R, g)))[1] for g in anti]
    linear = [sum(a[i] * xs(i, 1) for i ∈ 1:nvertices(Q))]

    return QuotientRing(A, std(Ideal(A, [tautological; linear])))
end

# TODO todd class
# TODO point class
# TODO universal bundle class















# TODO the case for stacks
function dimension(M::QuiverModuliSpace)
    if M.condition == "stable"
        return 1 - Euler_form(M.Q, M.d, M.d)
    elseif M.condition == "semistable"
        if has_semistables(M.Q, M.d, M.theta)
            return maximum(
                dimension_of_luna_stratum(tau) for tau in all_luna_types(M.Q, M.d, M.theta)
            )
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
function Chow_ring(
    M::QuiverModuli,
    a::Vector{Int} = linearization(M),
    standard::Bool = false,
)
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

function index(M::QuiverModuli)
    # right? 
    return gcd(canonical_stability(M.Q, M.d))
end
