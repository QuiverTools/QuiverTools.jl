export QuiverModuli, QuiverModuliSpace, QuiverModuliStack

abstract type QuiverModuli end


# TODO consider this: https://stackoverflow.com/questions/71738970/in-julia-declare-abstractvectorabstractvector
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
    denom::Function

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
	is_nonempty(M::QuiverModuli)

Checks if the quiver moduli is nonempty.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

OUTPUT:
- whether the moduli space is nonempty.

EXAMPLES:
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


"""
	is_coprime(M::QuiverModuli)

Checks if the stability parameter is coprime with the dimension vector,
i.e., if for all subdimension vectors `e` of `d`, \$\\theta\\cdot e \\neq 0\$.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

OUTPUT:
- whether the dimension vector `M.d` is theta-coprime for `M.theta`.

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> is_coprime(M)
true
```
"""
function is_coprime(M::QuiverModuli)
    return is_coprime(M.d, M.theta)
end

"""
	all_HN_types(M::QuiverModuli; unstable::Bool = false, ordered::Bool = true)

Returns all Harder-Narasimhan types of the moduli space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `unstable::Bool = false`: if `true`, returns only Harder-Narasimhan types
corresponding to unstable representations.
- `ordered::Bool = true`: if `true`, returns the Harder-Narasimhan types in
the order introduced by Reineke. TODO link to the paper.

OUTPUT:
- a list of Harder-Narasimhan types for the dimension vector and slope of `M`.

EXAMPLES:

The HN types for a 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> all_HN_types(M)
8-element Vector{Vector{AbstractVector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_HN_types(M, unstable = true)
7-element Vector{Vector{AbstractVector{Int64}}}:
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]
```
"""
function all_HN_types(M::QuiverModuli; unstable::Bool = false, ordered::Bool = true)
    HN = all_HN_types(M.Q, M.d, M.theta, M.denom, ordered)
    if unstable
        return filter(hn_type -> hn_type != [M.d], HN)
    end
	return HN
end

"""
	is_HN_type(M::QuiverModuli, hn_type::AbstractVector{<:AbstractVector{Int}})

Checks if the given sequence of dimension vectors is a valid HN type for
the moduli space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `hn_type::AbstractVector{<:AbstractVector{Int}}`: a sequence of dimension vectors.

OUTPUT:
- whether the given sequence is a valid Harder-Narasimhan type for `M`.

EXAMPLES:

Some HN types for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_HN_type(M, [[2, 3]])
true

julia> is_HN_type(M, [[1, 1], [1, 2]])
true

julia> is_HN_type(M, [[1, 2], [1, 1]])
false
```
"""
function is_HN_type(M::QuiverModuli, hn_type::Vector{<:AbstractVector{Int}})::Bool
    return is_HN_type(M.Q, M.d, hn_type, M.theta, M.denom)
end


"""
	codimension_HN_stratum(M::QuiverModuli, hn_type::Vector{<AbstractVector{Int}})

Computes the codimension of the Harder-Narasimhan stratum corresponding to the given HN type.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `hn_type::Vector{<:AbstractVector{Int}}`: a HN type for `M`.

OUTPUT:
- the codimension of the Harder-Narasimhan stratum corresponding to the given HN type.

EXAMPLES:

Codimensions for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> codimension_HN_stratum(M, [[2, 3]])
0

julia> codimension_HN_stratum(M, [[1, 1], [1, 2]])
"""
function codimension_HN_stratum(
    M::QuiverModuli,
    hn_type::Vector{<:AbstractVector{Int}},
)
    return codimension_HN_stratum(M.Q, hn_type)
end


"""
	codimension_unstable_locus(M::QuiverModuli)

Computes the codimension of the unstable locus in the parameter space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

OUTPUT:
- the codimension of the unstable locus in the parameter space.

EXAMPLES:

Codimensions for the 3-Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> codimension_unstable_locus(M)
3
```
"""
function codimension_unstable_locus(M::QuiverModuli)
    HN = all_HN_types(M, unstable = true)
    return minimum(codimension_HN_stratum(M, hn_type) for hn_type in HN)
end

# TODO add safety checks everywhere in the codebase


# TODO add examples of Luna types
"""
	all_Luna_types(M::QuiverModuli; exclude_stable::Bool = false)

Returns all Luna types of the moduli space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `exclude_stable::Bool = false`: if `true`, excludes the stable Luna type.

OUTPUT:
- a list of Luna types for the dimension vector and slope of `M`.

"""
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

    # treat the zero case separately
    if all(di == 0 for di in d)
        return [Dict(d => [1])]
    end

    # subdimensions with the same slope as d
    same_slope = filter(
        e ->
            slope(e, theta, denom) == slope(d, theta, denom) &&
                has_stables(Q, e, theta, denom),
        QuiverTools.all_subdimension_vectors(d, nonzero = true),
    )

    Luna_types = []

    # the highest possible amount of repetitions for a given stable dimension vector
    bound = sum(d) ÷ minimum(sum(e) for e in same_slope)
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
                    p in Iterators.product(values(partial)...)
                ]
                Luna_types = vcat(Luna_types, new_Luna)
            end
        end
    end

    if exclude_stable
        return filter(luna -> luna != Dict(d => [1]), Luna_types)
    end
    return Luna_types
end

# TODO add examples
"""
	is_Luna_type(M::QuiverModuli, tau)

Checks if the given tau is a valid Luna type for `M`.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`.

OUTPUT:
- whether the given tau is a valid Luna type for `M`.
"""
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

# TODO this should return 0 for the type Dict([0, 0] => [1])??
"""
	dimension_of_Luna_stratum(M::QuiverModuli, tau)

Computes the dimension of the Luna stratum corresponding to the given Luna type in the
moduli space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`.

OUTPUT:
- the dimension of the Luna stratum corresponding to the given Luna type.
"""
function dimension_of_Luna_stratum(M::QuiverModuli, tau)
    return sum(length(tau[e]) * (1 - Euler_form(M.Q, e, e)) for e in keys(tau))
end

"""
	local_quiver_setting(M::QuiverModuli, tau)

Returns the local quiver and dimension vector for the given Luna type.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`.

OUTPUT:
- a dictionary with the local quiver `Q` and dimension vector `d` for the given Luna type.
"""
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


"""
	semistable_equals_stable(M::QuiverModuli)

Checks if stability and semistability are equivalent on the given moduli space.
In other words, checks if there are no properly semistable points in the representation
space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

OUTPUT:
- whether every semistable representation is stable.

EXAMPLES:

If the dimension vector is coprime with the stability parameter, then semistability
and stability are equivalent:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> semistable_equals_stable(M)
true
```

However, this is not necessarily the case:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> semistable_equals_stable(M)
false
```
"""
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

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.

OUTPUT:
- whether the codimension of the unstable locus is at least `2`.

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_amply_stable(M)
true
```
"""
function is_amply_stable(M::QuiverModuli)
    return codimension_unstable_locus(M) >= 2
end




"""
    solve(A, b)

Solve ``A\\cdot x = b`` for ``A`` upper triangular via back substitution.

This is an internal method only used in the implementation of the Hodge polynomial
and to compute motives.


INPUT:
- `A::AbstractMatrix`: an upper triangular matrix.
- `b::AbstractVector`: a vector.

OUTPUT:
- the solution `x` to the equation.

EXAMPLES:

```jldoctest
julia> A = [1 2 3; 0 4 5; 0 0 6];

julia> b = [1, 2, 3];

julia> QuiverTools.solve(A, b)
3-element Vector{Any}:
 -0.25
 -0.125
  0.5
```
"""
function solve(A, b)
    n = length(b)
    x = Vector{Any}(zeros(n))

    x[n] = b[n] / A[n, n]

    for i in n-1:-1:1
        x[i] = (b[i] - sum(A[i, j] * x[j] for j in i+1:n)) / A[i, i]
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
        return prod(q^n - q^i for i in 0:n-1)
    end
end

"""
Cardinality of representation space \$\\mathrm{R}(Q,d)\$, over \$\\mathbb{F}_q\$.
"""
function CardinalRd(Q::Quiver, d::AbstractVector{Int}, q)
    return q^sum(d[i] * d[j] * Q.adjacency[i, j] for i in 1:nvertices(Q), j in 1:nvertices(Q))
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
        for j in i:l  # upper triangular
            T[i, j] = TransferMatrixEntry(Q, Ii, I[j], q)
        end
    end
    return T
end


# TODO DOI below is not open access.
# auxiliary functions for Hodge_polynomial() above
###################################################

"""
    Hodge_polynomial(Q, d, theta)

Returns the Hodge polynomial of the moduli space of ``\\theta``-semistable
representations of ``Q`` with dimension vector ``d``.

The algorithm is based on [MR1974891](https://doi.org/10.1007/s00222-002-0273-4),
and the current implementation is translated from the [Hodge diamond cutter]
(https://zenodo.org/doi/10.5281/zenodo.3893509).

INPUT:
- ``Q``: a quiver.
- ``d``: a dimension vector.
- ``theta``: a stability parameter. Default is the canonical stability.

OUTPUT:
- the Hodge polynomial of the moduli space.

EXAMPLES:

The Hodge polynomial of our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> d = [2, 3];

julia> theta = [3, -2];

julia> Hodge_polynomial(Q, d, theta)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1
```
"""
function Hodge_polynomial(Q::Quiver,
	d::AbstractVector{Int},
	theta::AbstractVector{Int} = canonical_stability(Q, d)
    )

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


"""
    Hodge_polynomial(M::QuiverModuliSpace)

Returns the Hodge polynomial of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the Hodge polynomial of the moduli space.

EXAMPLES:

The Hodge polynomial of our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> Hodge_polynomial(M)
x^6*y^6 + x^5*y^5 + 3*x^4*y^4 + 3*x^3*y^3 + 3*x^2*y^2 + x*y + 1
```
"""
function Hodge_polynomial(M::QuiverModuliSpace)
    return Hodge_polynomial(M.Q, M.d, M.theta)
end

"""
    Hodge_diamond(Q, d, theta)

Returns the Hodge diamond of the moduli space of
``\\theta``-semistable representations of ``Q`` with dimension vector ``d``.

INPUT:
- ``Q``: a quiver.
- ``d``: a dimension vector.
- ``theta``: a stability parameter. Default is the canonical stability.

OUTPUT:
- the Hodge diamond of the moduli space.

EXAMPLES:

The Hodge diamond of our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> Hodge_diamond(Q, [2, 3])
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```
"""
function Hodge_diamond(Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d)
    )::Matrix{Int}
    g = Hodge_polynomial(Q, d, theta)

    # collects the coefficients of the polynomial, converts them to integers
    # and returns them in the diagonal of a matrix.
    return diagonal(Int.(numerator.(collect(Singular.coefficients(g)))))
end
"""
    Hodge_diamond(M::QuiverModuliSpace)

Returns the Hodge diamond of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the Hodge diamond of the moduli space.

EXAMPLES:

The Hodge diamond of our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> Hodge_diamond(M)
7×7 Matrix{Int64}:
 1  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  3  0  0  0  0
 0  0  0  3  0  0  0
 0  0  0  0  3  0  0
 0  0  0  0  0  1  0
 0  0  0  0  0  0  1
```
"""
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


# TODO use betti numbers
function Picard_rank(M::QuiverModuli)
    return Picard_rank(M.Q, M.d, M.theta)
end


"""
    Betti_numbers(M::QuiverModuliSpace)

Returns the Betti numbers of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- a list of Betti numbers of the moduli space.

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> Betti_numbers(M)
3-element Vector{Int64}:
 1
 0
 1
```

Our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> Betti_numbers(M)
13-element Vector{Int64}:
 1
 0
 1
 0
 3
 0
 3
 0
 3
 0
 1
 0
 1
```
 """
function Betti_numbers(M::QuiverModuliSpace)

    if !is_coprime(M.d, M.theta)
        throw(ArgumentError("d and theta are not coprime"))
    end

    N = dimension(M)
    P = Poincare_polynomial(M)
    v = Singular.vars(P)[1]
    # silly way to square the variables but not the coefficients
    f = sum(term * v^Singular.degree(term, 1) for term in Singular.terms(P))

    betti = [f(0)]
    for i in 1:2*N
        f = (f - f(0)) / v
        push!(betti, f(0))
    end
    # returns in Ints
    return Int.(numerator.(betti))
end
"""
    Poincare_polynomial(M::QuiverModuliSpace)

Returns the Poincaré polynomial of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the Poincaré polynomial of the moduli space.

EXAMPLES:

A Kronecker quiver setup where ``M`` is the projective line:
```jldoctest
julia> Q = mKronecker_quiver(2);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> Poincare_polynomial(M)
L + 1
```

The Poincaré polynomial of our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> Poincare_polynomial(M)
L^6 + L^5 + 3*L^4 + 3*L^3 + 3*L^2 + L + 1
```
"""
function Poincare_polynomial(M::QuiverModuliSpace)
    if !is_coprime(M.d, M.theta)
        throw(ArgumentError("d and theta are not coprime"))
    end

    mot = unsafe_motive(M.Q, M.d, M.theta)
    P = (1 - mot["function_field"][2]) * mot["motive"] 

    if denominator(P) != 1
        throw(DomainError("must be a polynomial"))
    end
    # returns a polynomial object instead of a FunctionField element.
    return Singular.n_transExt_to_spoly(numerator(P))
end


# oh my god
function power(x, n::Int)
    if n >= 0
        return x^n
    else
        return 1 / x^(-n)
    end
end

function motive(M::QuiverModuliStack)

    if M.condition != "stable"
        throw(ArgumentError("Motive unknown if not stable"))
    end
    return unsafe_motive(M.Q, M.d, M.theta)["motive"]
end

"""
    unsafe_motive(Q, d, theta, denom)

Returns the motive of the moduli stack of ``\\theta``-semistable representations, and
the function field.

This is an internal method. Use ``motive()`` instead.

INPUT:
- ``Q``: a quiver.
- ``d``: a dimension vector.
- ``theta``: a stability parameter. Default is the canonical stability.
- ``denom``: a function. Default is the sum.

OUTPUT:
- a dictionary with the motive and the function field.

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3);

julia> QuiverTools.unsafe_motive(Q, [2, 3])
Dict{String, Any} with 2 entries:
  "motive"         => (-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)//(L - 1)
  "function_field" => (Function field over Rational field with transcendence ba…

julia> QuiverTools.unsafe_motive(Q, [2, 3])["motive"]
(-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)//(L - 1)
```
"""
function unsafe_motive(Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum
    )


    K, L = Singular.FunctionField(Singular.QQ, ["L"])
    L = L[1]

    if all(ti == 0 for ti in theta)
        out = power(L,(- Euler_form(Q, d, d)))
        div = 1
        for i in 1:nvertices(Q)
            if d[i] > 0
                div *= prod(1 - power(L, -nu) for nu in 1:d[i])
            end
        end
        return out/div
    end

    ds = all_subdimension_vectors(d, nonzero = true, strict = true)
    ds = filter(e -> slope(e, theta, denom) > slope(d, theta, denom), ds)

    push!(ds, zero_vector(nvertices(Q)), d)
    sort!(ds, by = e -> deglex_key(Q, e)) #hopefully

    T = Matrix{Any}(undef, length(ds), length(ds))
    for (i, j) in Iterators.product(1:length(ds), 1:length(ds))
       if is_subdimension_vector(ds[i], ds[j])
           T[i, j] = power(L, Euler_form(Q, ds[i] - ds[j], ds[i])) *
           unsafe_motive(Q, ds[j] - ds[i], zero_vector(nvertices(Q)))
       else
            T[i, j] = 0
       end
    end

    y = [0 for i in 1:length(ds)]
    y[end] = 1
    y = coerce_vector(y)

    return Dict("motive" => solve(T, y)[1],  "function_field" => (K, L))
end

















###############################################################################
# tautological representation of the Chow ring.
# Implements the results of [arXiv:1307.3066](https://doi.org/10.48550/arXiv.1307.3066) and
# [arXiv.2307.01711](https://doi.org/10.48550/arXiv.2307.01711).
###############################################################################

# partial order on the forbidden dimension vectors as in https://doi.org/10.48550/arXiv.1307.3066
function partial_order(Q::Quiver, f::AbstractVector{Int}, g::AbstractVector{Int})
    if !all(f[i] <= g[i] for i in 1:nvertices(Q) if is_source(Q, i))
        return false
    elseif !all(f[i] >= g[i] for i in 1:nvertices(Q) if is_sink(Q, i))
        return false
    elseif !all(f[i] == g[i] for i in 1:nvertices(Q) if !is_source(Q, i) && !is_sink(Q, i))
        return false
    end
    return true
end


"""
	symmetric_polynomial(vars, degree)

Returns the symmetric polynomial of degree ``degree`` in the variables ``vars``.

INPUT:
- ``vars``: a list of variables.
- ``degree``: the degree of the wanted symmetric polynomial.

OUTPUT:
- The symmetric polynomial of degree ``degree`` in the variables ``vars``.

EXAMPLES:

```jldoctest
julia> using Singular;

julia> R, vars = polynomial_ring(Singular.QQ, ["x", "y", "z"]);

julia> QuiverTools.symmetric_polynomial(vars, 2)
x*y + x*z + y*z
```
"""
function symmetric_polynomial(vars, degree::Int)
    return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

function Chow_ring(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int}=canonical_stability(Q, d),
    a::AbstractVector{Int}=extended_gcd(d)[2],
)
    # TODO cover case d[i] = 0
    # safety checks
    if !is_coprime(d, theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif a' * d != 1
        throw(ArgumentError("a is not a linearization"))
    end

    varnames = ["x$i$j" for i in 1:nvertices(Q) for j in 1:d[i] if d[i] > 0]
    R, vars = polynomial_ring(Singular.QQ, varnames)
    function chi(i, j)
        return vars[sum(d[1:i-1])+j]
    end

	"""
	This is the naive base that is described in Hans's 2013 paper.
	"""
    function base_for_ring(name = "naive")
        if name == "naive"
            bounds = [0:(d[i]-nu) for i in 1:nvertices(Q) for nu in 1:d[i]]
            lambdas = Iterators.product(bounds...)

            build_elem(lambda) = prod(
                prod(chi(i, nu)^lambda[sum(d[1:i-1])+nu] for nu in 1:d[i]) for
                i in 1:nvertices(Q)
            )

            return map(l -> build_elem(l), lambdas)
        else
            throw(ArgumentError("unknown base."))
        end
    end

    # build the permutation group W
    W = Iterators.product([AbstractAlgebra.SymmetricGroup(d[i]) for i in 1:nvertices(Q)]...)
    sign(w) = prod(AbstractAlgebra.sign(wi) for wi in w)

    permute(f, sigma) = f([chi(i, sigma[i][j]) for i in 1:nvertices(Q) for j in 1:d[i]]...)

    delta = prod(
        prod(chi(i, l) - chi(i, k) for k in 1:d[i]-1 for l in k+1:d[i]) for
        i in 1:nvertices(Q) if d[i] > 1
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
            prod((chi(j, s) - chi(i, r))^Q.adjacency[i, j]
			for r in 1:e[i], s in e[j]+1:d[j])
				for j in 1:nvertices(Q), i in 1:nvertices(Q) if
            Q.adjacency[i, j] > 0 && e[i] > 0 && d[j] > 1
        ) for e in all_forbidden(Q, d, theta)
    ]

    varnames2 = ["x$i$j" for i in 1:nvertices(Q) for j in 1:d[i] if d[i] > 0]
    A, Avars = polynomial_ring(Singular.QQ, varnames2)

    function xs(i, j)
        return Avars[sum(d[1:i-1])+j]
    end

    targets = [
        [symmetric_polynomial([chi(i, j) for j in 1:d[i]], k) for k in 1:d[i]] for
        i in 1:nvertices(Q)
    ]
    targets = reduce(vcat, targets)

    inclusion = AlgebraHomomorphism(A, R, targets)

    anti = [antisymmetrize(f * b) for f in forbidden_polynomials for b in base_for_ring()]
    tautological = [gens(preimage(inclusion, Ideal(R, g)))[1] for g in anti]
    linear = [sum(a[i] * xs(i, 1) for i in 1:nvertices(Q))]

    return QuotientRing(A, std(Ideal(A, [tautological; linear])))
end


# TODO this should be in a misc.jl file or something
# TODO this should be in Base really...

"""
Computes the gcd and the Bezout coefficients of a list of integers.

This exists for two integers but seemingly not for more than two.
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
		y = vcat([g],  [x[i] for i in 3:n])
		d, c = extended_gcd(y)
		m = vcat([c[1] * a, c[1] * b], [c[i] for i in 1:n - 1])
		return [d, m]
	end
end

# TODO todd class
function todd_class(Q::Quiver,
	d::AbstractVector{Int},
	chi::AbstractVector{Int}=extended_gcd(d)[2]
	)

	"""
	We call the series \$Q(t) = t/(1-e^{-t})\$ the Todd generating series.
	The function computes the terms of this series up to degree n.
	We use this instead of the more conventional notation `Q` to avoid a
	clash with the notation for the quiver.
	"""
	function todd_Q(t, n)
		return sum(
			(-1) ^ i * (bernoulli(i) * t^i) / factorial(i) for i in 0:n
		)
	end

	"""
	Takes an element in a graded ring and discards all homogeneous components
	of degree > n
	"""
	function truncate(f, n)
		return sum([term for term in terms(f) if total_degree(term) <= n])
	end

	throw(NotImplementedError())
    # TODO refactor Chow ring to return ideal, rings, inclusion, variables
    # TODO then use that here
	num = 1
	den = 1

	for a in Q.arrows()
		i, j = a
		for p in 1:d[i]
			for q in 1:d[j]
				num *= todd_Q(chi(j, q) - chi(i, p), sum(d))
				num = truncate(num, sum(d))
			end
		end
	end

	for i in 1:nvertices(Q)
		for p in 1:d[i]
			for q in 1:d[i]
				den *= todd_Q(chi(i, q) - chi(i, p), sum(d))
				den = truncate(den, sum(d))
			end
		end
	end

	return num / den
end

# TODO point class
# TODO universal bundle class











"""
    dimension(M::QuiverModuliStack)

Returns the dimension of the moduli stack.
This differs from the dimension of the moduli space by 1, as we do not quotient out
the stabilizer \$ \\mathbb{G}\$.

INPUT:
- ``M``: a moduli stack of representations of a quiver.

OUTPUT:
- the dimension of the moduli stack.

EXAMPLES:

The dimension of the moduli stack of the 3-Kronecker quiver

```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3]);

julia> dimension(M)
5
```
"""
function dimension(M::QuiverModuliStack)
    if is_nonempty(M)
        return -Euler_form(M.Q, M.d, M.d)        
    end
    return "-∞"
end

"""
    dimension(M::QuiverModuliSpace)

Returns the dimension of the moduli space.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the dimension of the moduli space.

EXAMPLES:

The dimension of the moduli space of the 3-Kronecker quiver
with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> dimension(M)
6
```
"""
function dimension(M::QuiverModuliSpace)
    # the zero representation is semistable, but not stable, for d = 0
    if all(M.d .== 0)
        if M.condition == "semistable"
            return 0
        else
            return "-∞"
        end
    end

    # if the stable locus is nonempty then the dimension is 1 - <d, d>
    if has_stables(M.Q, M.d, M.theta, M.denom)
        return 1 - Euler_form(M.Q, M.d, M.d)
    end

    # if the stable locus is empty, the dimension is the maximum of the dimensions
    # of the Luna strata
    if M.condition == "stable"
        return "-∞"
    elseif M.condition == "semistable"
        if has_semistables(M.Q, M.d, M.theta)
            return maximum(
                dimension_of_Luna_stratum(M, tau) for tau in all_Luna_types(M.Q, M.d, M.theta)
            )
            # TODO what are the Luna strata for d = 0?
            # shouldn't the case d = 0 be handled correctly here?

        end
    end
    # the semistable locus is also empty
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
