export QuiverModuli, QuiverModuliSpace, QuiverModuliStack

export Hodge_diamond, Hodge_polynomial, Picard_rank

export Chow_ring, motive, index, Betti_numbers, Poincare_polynomial,
    is_smooth, is_projective, semisimple_moduli_space, point_class,
    Todd_class, Chern_class_line_bundle, Chern_character_line_bundle,
    total_Chern_class_universal, integral
export all_Luna_types, is_Luna_type, dimension_of_Luna_stratum

export is_nonempty, dimension, is_smooth, semistable_equals_stable,
    codimension_unstable_locus




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
the order introduced by [MR1974891](https://doi.org/10.1007/s00222-002-0273-4).

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

Computes the codimension of the Harder-Narasimhan stratum
corresponding to the given HN type.

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
3
```
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


"""
	all_Luna_types(M::QuiverModuli; exclude_stable::Bool = false)

Returns all Luna types of the moduli space.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `exclude_stable::Bool = false`: if `true`, excludes the stable Luna type.

OUTPUT:
- a list of Luna types for the dimension vector and slope of `M`.

EXAMPLES:

Luna types for a 3-Kronecker quiver:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> all_Luna_types(M)
5-element Vector{Dict{AbstractVector, Vector{Int64}}}:
 Dict([3, 3] => [1])
 Dict([1, 1] => [1], [2, 2] => [1])
 Dict([1, 1] => [3])
 Dict([1, 1] => [2, 1])
 Dict([1, 1] => [1, 1, 1])
```
"""
function all_Luna_types(M::QuiverModuli; exclude_stable::Bool = false)
    return all_Luna_types(M.Q, M.d, M.theta, M.denom, exclude_stable)
end

"""
    all_Luna_types(Q, d, theta, denom, exclude_stable)

Computes all the possible Luna types for the given data.

INPUT:
- `Q`: a quiver.
- `d`: a dimension vector.
- `theta`: a stability parameter. Defaults to the canonical stability.
- `denom`: a function defining the denominator of the slope. Defaults to sum.
- `exclude_stable`: whether to exclude the Luna type of stable representations.

OUTPUT:
A list of Luna types.


EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> all_Luna_types(M)
5-element Vector{Dict{AbstractVector, Vector{Int64}}}:
 Dict([3, 3] => [1])
 Dict([1, 1] => [1], [2, 2] => [1])
 Dict([1, 1] => [3])
 Dict([1, 1] => [2, 1])
 Dict([1, 1] => [1, 1, 1])

julia> X = QuiverModuliSpace(Q, [2, 3]);

julia> all_Luna_types(X)
1-element Vector{Dict{AbstractVector, Vector{Int64}}}:
 Dict([2, 3] => [1])
```
"""
function all_Luna_types(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum,
    exclude_stable::Bool = false,
)::Vector{Dict{AbstractVector, Vector{Int}}}

    d = coerce_vector(d)
    theta = coerce_vector(theta)

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
                    if e in collect(keys(partial))
                        partial[e] += 1
                    else
                        partial[e] = 1
                    end
                end

                for e in keys(partial)
                    partial[e] = partitions(partial[e])
                end

                for p in Iterators.product(values(partial)...)
                    new_Luna_type = Dict(zip(collect(keys(partial)), p))
                    push!(Luna_types, new_Luna_type)
                end
            end
        end
    end

    if exclude_stable
        return filter(luna -> luna != Dict(d => [1]), Luna_types)
    end
    return Luna_types
end

"""
	is_Luna_type(M::QuiverModuli, tau)

Checks if the given tau is a valid Luna type for `M`.

INPUT:
- `M::QuiverModuli`: a moduli space or stack of representations of a quiver.
- `tau::Dict{AbstractVector{Int}, Vector{Int}}`: a Luna type for `M`.

OUTPUT:
- whether the given tau is a valid Luna type for `M`.

EXAMPLES:

Nontrivial Luna types for the 3-Kronecker quiver:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [3, 3]);

julia> l = Dict([1, 1] => [1], [2, 2] => [1]);

julia> is_Luna_type(M, l)
true
```

The zero dimensional case:
```jldoctest
julia> Q = mKronecker_quiver(3); X = QuiverModuliSpace(Q, [0, 0]);

julia> is_Luna_type(X, Dict([0, 0] => [1]))
true
```
"""
function is_Luna_type(M::QuiverModuli, tau)
    if sum(M.d) == 0
        return tau == Dict(M.d => [1])
    end

    ks = collect(keys(tau))
    if sum(ks) != M.d
        return false
    end
    if !all(slope(e, M.theta, M.denom) == slope(M.d, M.theta, M.denom) for e in ks)
        return false
    end

    if !all(has_semistables(M.Q, e, M.theta, M.denom) for e in ks)
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

EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(2); M = QuiverModuliSpace(Q, [2, 2], [1, -1]);

julia> luna = all_Luna_types(M)
2-element Vector{Dict{AbstractVector, Vector{Int64}}}:
 Dict([1, 1] => [2])
 Dict([1, 1] => [1, 1])

julia> [dimension_of_Luna_stratum(M, tau) for tau in luna]
2-element Vector{Int64}:
 1
 2
```
"""
function dimension_of_Luna_stratum(M::QuiverModuli, tau)
    return sum(length(tau[e]) * (1 - Euler_form(M.Q, e, e)) for e in collect(keys(tau)))
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

    ks = collect(keys(tau))
    A = coerce_matrix([
        [generic_ext(M.Q, e, eprime) for eprime in ks for n in tau[eprime]] for
        e in ks for m in tau[e]
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


########################################################################################
# Below lie methods to compute Hodge diamonds translated from the Hodge diamond cutter.
# In turn, these are based on M. Reineke's paper
# "The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli", 
# https://doi.org/10.1007/s00222-002-0273-4
########################################################################################


###################################################
# auxiliary functions for Hodge_polynomial() below

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
    return q^sum(
                d[i] * d[j] * Q.adjacency[i, j]
                for i in 1:nvertices(Q),
                    j in 1:nvertices(Q))
end

"""
Cardinality of product of general linear groups \$\\mathrm{GL}_{d}(\\mathbb{F}_q)\$.
"""
@memoize Dict function CardinalGd(d::AbstractVector{Int}, q)
    return prod(CardinalGl(di, q) for di in d)
end


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
    Picard_rank(M::QuiverModuliSpace)

Returns the Picard rank of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the Picard rank of the moduli space.

EXAMPLES:

Kronecker quiver with dimension vector `[2, 3]`:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> Picard_rank(M)
1
```
"""
function Picard_rank(M::QuiverModuliSpace)
    if !(is_smooth(M) && is_projective(M))
        throw(ArgumentError("Moduli space is not smooth and projective"))
    end
    return Betti_numbers(M)[3]
end


"""
    index(M::QuiverModuliSpace)

Returns the index of the moduli space ``M``.

The index of a variety \$X\$ is the largest which divides the canonical divisor \$K_X\$
in \$Pic(X)\$.

This implementation currently only works for the canonical stability.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the index of the moduli space.

EXAMPLES:

The 3-Kronecker quiver has index 3:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> index(M)
3
```
The subspace quiver moduli have index 1:
```jldoctest
julia> Q = subspace_quiver(5);

julia> M = QuiverModuliSpace(Q, [1, 1, 1, 1, 1, 2]);

julia> index(M)
1
```
"""
function index(M::QuiverModuliSpace)
    if M.theta == canonical_stability(M.Q, M.d) &&
        is_coprime(M.d, M.theta) &&
        is_amply_stable(M)
        return gcd(M.theta)
    end
    throw(NotImplementedError("Only implemented for canonical stability."))
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
    coeff = Int.(numerator.(Singular.coefficients(P)))
    betti = reduce(vcat, [c, 0] for c in coeff[1:end-1])
    push!(betti, coeff[end])

    # if the polynomial did not have degree = N,
    # we add zero coefficients
    if length(betti) < 2*N + 1
        betti = vcat(betti, zeros(2*N + 1 - length(betti)))
    end
    return betti
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

    m = motive(M.Q, M.d, M.theta, M.denom)
    v = Singular.transcendence_basis(Singular.parent(m))[1]
    P = (1 - v) * m 

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
    return motive(M.Q, M.d, M.theta)
end

"""
    motive(Q, d, theta, denom)

Returns the motive of the moduli stack of ``\\theta``-semistable representations.

INPUT:
- ``Q``: a quiver.
- ``d``: a dimension vector.
- ``theta``: a stability parameter. Default is the canonical stability.
- ``denom``: a function. Default is the sum.

OUTPUT:
- The motive as an element in the function field \$\\mathbb{Q}(L)\$.

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3);

julia> motive(Q, [2, 3])
(-L^6 - L^5 - 3*L^4 - 3*L^3 - 3*L^2 - L - 1)//(L - 1)
```
"""
function motive(Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum
    )


    K, L = Singular.FunctionField(Singular.QQ, ["L"])
    L = L[1]

    if all(ti == 0 for ti in theta)
        out = power(L,- Euler_form(Q, d, d))
        den = 1
        for i in 1:nvertices(Q)
            if d[i] > 0
                den *= prod(1 - power(L, -nu) for nu in 1:d[i])
            end
        end
        return out / den
    end

    ds = all_subdimension_vectors(d, nonzero = true, strict = true)
    ds = filter(e -> slope(e, theta, denom) > slope(d, theta, denom), ds)

    push!(ds, zero_vector(nvertices(Q)), d)
    sort!(ds, by = e -> deglex_key(Q, e)) #hopefully

    T = Matrix{Any}(undef, length(ds), length(ds))
    for (i, j) in Iterators.product(1:length(ds), 1:length(ds))
       if is_subdimension_vector(ds[i], ds[j])
           T[i, j] = power(L, Euler_form(Q, ds[i] - ds[j], ds[i])) *
           motive(Q, ds[j] - ds[i], zero_vector(nvertices(Q)))
       else
            T[i, j] = 0
       end
    end

    y = [0 for i in 1:length(ds)]
    y[end] = 1
    y = coerce_vector(y)

    return solve(T, y)[1]
end

###############################################################################
# tautological representation of the Chow ring.
# Implements the results of [arXiv:1307.3066](https://doi.org/10.48550/arXiv.1307.3066)
# and [arXiv.2307.01711](https://doi.org/10.48550/arXiv.2307.01711).
###############################################################################

# partial order on the forbidden dimension vectors as defined in
# https://doi.org/10.48550/arXiv.1307.3066
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

```julia-repl
julia> using Singular;

julia> R, vars = polynomial_ring(Singular.QQ, ["x", "y", "z"]);

julia> QuiverTools.symmetric_polynomial(vars, 2)
x*y + x*z + y*z
```
"""
function symmetric_polynomial(vars, degree::Int)
    return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

# TODO test Chow ring
"""
    Chow_ring(Q, d, theta, a)

Computes the Chow ring of the moduli space of ``\\theta``-semistable representations of
``Q`` with dimension vector ``d``, for a choice of linearization ``a``.

This method of the function Chow_ring also returns the ambient ring \$R\$
and the inclusion morphism.

INPUT:
- ``Q``: a quiver.
- ``d``: a dimension vector.
- ``theta``: a stability parameter. Default is the canonical stability.
- ``a``: a linearization. Default is the extended gcd of ``d``.

OUTPUT:
A tuple containing:
- the Chow ring of the moduli space,
- the polynomial ring above it,
- the inclusion map \$\\iota : A \\to R\$.
"""
@memoize Dict function Chow_ring(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int}=canonical_stability(Q, d),
    a::AbstractVector{Int}=extended_gcd(d)[2],
)
    # safety checks
    if !is_coprime(d, theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif a' * d != 1
        throw(ArgumentError("a is not a linearization"))
    end

    # j varies first, then i
    varnames = ["xi$i$j" for i in 1:nvertices(Q) for j in 1:d[i]]
    R, vars = polynomial_ring(Singular.QQ, varnames)

    # Shorthand to address the variable `xi_{i,j}`.
    function xi(i, j)
        if d[i] == 0
            throw(ArgumentError("i is not in the support of d."))
        end
        return vars[sum(d[1:i-1]) + j]
    end


	# This is the naive base that is described in Hans's 2013 paper.
    function base_for_ring()

        bounds = [0:(d[i]-nu) for i in 1:nvertices(Q) for nu in 1:d[i]]
        lambdas = Iterators.product(bounds...)

        build_elem(lambda) = prod(
            prod(xi(i, nu)^lambda[sum(d[1:i-1])+nu] for nu in 1:d[i])
            for i in 1:nvertices(Q)
                if d[i] > 0
        )
        return map(l -> build_elem(l), lambdas)
    end

    # build the permutation group W
    W = Iterators.product([
                            AbstractAlgebra.SymmetricGroup(d[i])
                            for i in 1:nvertices(Q)
                                ]...)
    sign(w) = prod(AbstractAlgebra.sign(wi) for wi in w)

    # Action of the symmetric group on R by permutation of the variables.
    permute(f, sigma) = f([xi(i, sigma[i][j])
                            for i in 1:nvertices(Q)
                            for j in 1:d[i]
                                if d[i] > 0]...)

    # The discriminant in the definition of the antisymmetrization.
    delta = 1
    for i in 1:nvertices(Q)
        if d[i] > 1
            delta *= prod(xi(i, l) - xi(i, k) for k in 1:d[i]-1 for l in k+1:d[i])
        end
    end

    antisymmetrize(f) = div(sum(sign(w) * permute(f, w) for w in W), R(delta))

    # All the destabilizing subdimension vectors of `d` with respect to the slope
    # `theta/denom` that are minimal with respect to the total order.
    dest = all_destabilizing_subdimension_vectors(d, theta)
    minimal_forbidden = filter(
        e -> !any(f -> partial_order(Q, f, e), filter(f -> f != e, dest)), dest)

    # builds a new forbidden polynomial for the minimal forbidden dimension vector e.
    function new_forbidden(e::AbstractVector{Int})
        out = 1
        for (i, j) in Iterators.product(1:nvertices(Q), 1:nvertices(Q))
                for r in 1:e[i], s in e[j]+1:d[j]
                    out *= (xi(j, s) - xi(i, r))^Q.adjacency[i, j]
                end
            end
        return out
    end
    forbidden_polynomials = [new_forbidden(e) for e in minimal_forbidden]

    varnames2 = ["x$i$j" for i in 1:nvertices(Q) for j in 1:d[i]]
    A, Avars = polynomial_ring(Singular.QQ, varnames2)

    # Shorthand to address the variables of A `x_{i,j}`.
    function xs(i, j)
        if d[i] == 0
            throw(ArgumentError("i is not in the support of d."))
        end
        return Avars[sum(d[1:i-1])+j]
    end

    targets = [
        [symmetric_polynomial([xi(i, j) for j in 1:d[i]], k) for k in 1:d[i]]
                                for i in 1:nvertices(Q)
                                    if d[i] > 0
    ]
    targets = reduce(vcat, targets)

    inclusion = AlgebraHomomorphism(A, R, targets)

    anti = [antisymmetrize(f * b) for f in forbidden_polynomials for b in base_for_ring()]
    tautological = [gens(preimage(inclusion, Ideal(R, g)))[1] for g in anti]
    linear = [sum(a[i] * xs(i, 1) for i in 1:nvertices(Q) if d[i] > 0)]

    return (QuotientRing(A, std(Ideal(A, [tautological; linear]))), R, inclusion)
end

"""
    Chow_ring(M::QuiverModuliSpace, chi)

Computes the Chow ring of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``chi``: a choice of linearization for the trivial line bundle.
    It picks one by default if not provided.


OUTPUT:
- the Chow ring of the moduli space.
"""
function Chow_ring(M::QuiverModuliSpace, chi::AbstractVector{Int}=extended_gcd(M.d)[2])
    return Chow_ring(M.Q, M.d, M.theta, chi)[1]
end


# this should be in a misc.jl file or something
# this should be in Base really...

"""
Computes the gcd and the Bezout coefficients of a list of integers.

INPUT:
- ``x``: a list of integers.

OUTPUT:
A tuple containing:
- the gcd of the integers,
- a choice of Bezout coefficients.

EXAMPLES:

```jldoctest
julia> QuiverTools.extended_gcd([2, 3, 4])
2-element Vector{Any}:
 1
  [-1, 1, 0]
  
julia> QuiverTools.extended_gcd([2, 3])
2-element Vector{Any}:
 1
  [-1, 1]
```
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
		m = vcat([c[1] * a, c[1] * b], [c[i] for i in 2:n - 1])
		return [d, m]
	end
end

"""
    Chern_class_line_bundle(M::QuiverModuliSpace, eta)

Returns the first Chern class of the line bundle L(eta).

This is given by \$L(eta) = \\bigoplus_{i \\in Q_0} \\det(U_i)^{-eta_i}\$.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``eta``: a choice of linearization for the trivial line bundle.

OUTPUT:
- the first Chern class of the line bundle L(eta) as a polynomial.

EXAMPLES:

The line bundles \$\\mathcal{O}(i)\$ on the projective line:
```jldoctest
julia> Q = mKronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> l = Chern_class_line_bundle(M, [1, -1])
-x11
```

The line bundle corresponding to the canonical stability condition on our favourite
6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> Chern_class_line_bundle(M, [9, -6])
-3*x21
```
"""
function Chern_class_line_bundle(M::QuiverModuliSpace,
    eta::AbstractVector{Int})

    A, vars = Chow_ring(M)
    I = quotient_ideal(A)
    Rvars = gens(base_ring(I))

    Chern_class = -sum( eta[i] * Rvars[1 + sum(M.d[1:i-1])] for i in 1:nvertices(M.Q))

    return coerce_to_quotient(A, Chern_class)
end

"""
    Chern_character_line_bundle(M::QuiverModuliSpace, eta)

Returns the Chern character of the line bundle L(eta).

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``eta``: a choice of linearization for the trivial line bundle.

OUTPUT:
- the Chern character of the line bundle L(eta).

EXAMPLES:

Some line bundles on the projective line:
```jldoctest
julia> Q = mKronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> Chern_character_line_bundle(M, [1, -1])
-x11 + 1
```
Some Chern characters for our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> Chern_character_line_bundle(M, [3, -2])
1//720*x21^6 - 1//120*x21^5 + 1//24*x21^4 - 1//6*x21^3 + 1//2*x21^2 - x21 + 1
```
"""
function Chern_character_line_bundle(M::QuiverModuliSpace,
    eta::AbstractVector{Int})

    x = Chern_class_line_bundle(M, eta)
    Chern_character = sum(x^i / factorial(i) for i in 0:dimension(M))

    return Chern_character
end

"""
    total_Chern_class_universal(M::QuiverModuliSpace, i, chi)

Returns the total Chern class of the universal bundle \$U_i(\\chi)\$.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``i``: the universal bundle we want the Chern class of.
- ``chi``: a choice of linearization to construct \$U_i(\\chi)\$.

OUTPUT:
- the total Chern class of the universal bundle \$U_i(\\chi)\$.

EXAMPLES:

The universal Chern classes on both vertices of our favourite 3-Kronecker quiver:
```jldoctest
julia> Q  = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> total_Chern_class_universal(M, 1)
x11 + x12 + 1

julia> total_Chern_class_universal(M, 2)
x21 + x22 + x23 + 1
```
"""
function total_Chern_class_universal(M::QuiverModuliSpace,
    i::Int,
    chi::AbstractVector{Int} = extended_gcd(M.d)[2])

    A, Avars = Chow_ring(M, chi)    
    cUi = sum(
        Avars[sum(M.d[1:i - 1]) + r]
        for r in 1:M.d[i]; init=0
            ) + 1
    return cUi
end


"""
    point_class(M::QuiverModuliSpace)

Returns the point class of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``chi``: a choice of linearization to construct the universal bundles.

OUTPUT:
- the point class of the moduli space, as a polynomial in its Chow ring.

EXAMPLES:

A projective 7-fold:
```jldoctest
julia> Q = mKronecker_quiver(8);

julia> M = QuiverModuliSpace(Q, [1, 1]);

julia> point_class(M, [1, 0])
x21^7
```

Our favourite 6-fold:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> M = QuiverModuliSpace(Q, [2, 3]);

julia> point_class(M)
x23^2
```
"""
@memoize Dict function point_class(M::QuiverModuliSpace,
    chi::AbstractVector{Int} = extended_gcd(M.d)[2],
    )
    num = 1
    den = 1
    N = dimension(M)

    for i in 1:nvertices(M.Q)
        c = total_Chern_class_universal(M, i, chi)
        num *= c^(M.d' * M.Q.adjacency[:, i])
        den *= c^M.d[i]
    end

    quot = div(num, den)
    return sum(term for term in Singular.terms(quot)
                    if __Chow_ring_monomial_grading(M, term) == N; init = 0) 
end

"""
    todd_class(M::QuiverModuliSpace, chi)

Returns the Todd class of the moduli space ``M``.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- ``chi``: a choice of linearization to construct the universal bundles.

OUTPUT:
- the Todd class of the moduli space, as a polynomial in its Chow ring.

EXAMPLES:
The Todd class of our favourite 3-Kronecker quiver moduli:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> Todd_class(M)
-17//8*x12*x21 + x21^2 + 823//360*x12*x22 - 823//1080*x22^2 + 553//1080*x21*x23 - 77//60*x22*x23 + x23^2 + 5//12*x12 - 3//2*x21 + 9//8*x23 + 1
```
"""
@memoize Dict function Todd_class(M::QuiverModuliSpace,
	chi::AbstractVector{Int}=extended_gcd(M.d)[2]
	)

	"""
	We call the series \$Q(t) = t/(1-e^{-t})\$ the Todd generating series.
	The function computes the terms of this series up to degree n.
	We use this instead of the more conventional notation `Q` to avoid a
	clash with the notation for the quiver.
	"""
	function todd_Q(t, n)
		return sum(
			(-1)^i * (Nemo.bernoulli(i) * t^i) / factorial(i) for i in 0:n
		)
	end

	"""
	Takes an element in a graded ring and discards all homogeneous components
	of degree > n
	"""
	function truncate(f, n)
		return sum(
            term for term in Singular.terms(f)
                if Singular.total_degree(term) <= n
                    )
	end

    N = dimension(M)

    A, R, inclusion = Chow_ring(M.Q, M.d, M.theta, chi)
    Rvars = gens(R)

    function xi(i, p)
        return Rvars[sum(M.d[1:i-1]) + p]
    end

	num = 1
	den = 1

	for a in arrows(M.Q)
		i, j = a
		for p in 1:M.d[i]
			for q in 1:M.d[j]
				num *= todd_Q(xi(j, q) - xi(i, p), N)
				num = truncate(num, N)
			end
		end
	end

	for i in 1:nvertices(M.Q)
		for p in 1:M.d[i]
			for q in 1:M.d[i]
				den *= todd_Q(xi(i, q) - xi(i, p), N)
				den = truncate(den, N)
			end
		end
	end

    # this is because Singular does not have a method to get the preimage
    # of a given element, only ideals.
    # In Singular's implementation this does not result in a loss of time anyways...
    num = gens(preimage(inclusion, Ideal(R, num)))[1]
    den = gens(preimage(inclusion, Ideal(R, den)))[1]

    # renormalizing the constant term because Singular is silly like that
    num /= constant_coefficient(num)
    den /= constant_coefficient(den)

    num = coerce_to_quotient(A[1], num)
    den = coerce_to_quotient(A[1], den)

    return div(num, den)
end

"""
    integral(M, f, chi)

Computes the integral of `f` according to the Hirzebruch-Riemann-Roch theorem.

In other words, it computes the Euler characteristic of the vector bundle
whose Chern character is `f`.

INPUT:
- ``M``: a moduli space of representations of a quiver.
- `f`: the Chern character in CH(M) to integrate.
- ``chi``: a choice of linearization to construct the universal bundles.

OUTPUT:
the integral of `f`.

EXAMPLES:

The integral of \$\\mathcal{O}(i)\$ on the projective line for some `i`s.

```jldoctest
julia> Q = mKronecker_quiver(2); M = QuiverModuliSpace(Q, [1, 1]);

julia> L = Chern_character_line_bundle(M, [1, -1]);

julia> [integral(M, L^i) for i in 0:5]
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
```

Hilbert series for the 3-Kronecker quiver as in our favourite 6-fold:

```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> L = Chern_character_line_bundle(M, [3, -2]);

julia> [integral(M, L^i) for i in 0:5]
6-element Vector{Int64}:
    1
   20
  148
  664
 2206
 5999
```
"""
function integral(M::QuiverModuliSpace,
    f,
    chi::AbstractVector{Int}= extended_gcd(M.d)[2],
    )

    N = dimension(M)
    integrand = sum(
                    t
                    for t in collect(Singular.terms(f * Todd_class(M, chi)))
                    if __Chow_ring_monomial_grading(M, t) == N;
                    init = 0
                )

    integ = div(integrand, point_class(M, chi))
    # coercion to Int
    return Int(numerator(Singular.constant_coefficient(integ)))
end


"""
Takes a quotient ring R/I and a polynomial f in R and returns the image of f in R/I.
"""
function coerce_to_quotient(R, f)
    I = quotient_ideal(R)
    q = Singular.reduce(f, I)

    B = base_ring(R)
    g = Singular.MPolyBuildCtx(R)
    for (c, e) = zip(Singular.coefficients(q), Singular.exponent_vectors(q))
        Singular.push_term!(g, B(c), e)
    end
    return Singular.finish(g)
end

"""
Takes a ring R and a polynomial in R/I and returns the canonical preimage of f in R.
"""
function pullback_from_quotient(R, f)
    B = base_ring(R)
    g = Singular.MPolyBuildCtx(R)
    for (c, e) = zip(Singular.coefficients(f), Singular.exponent_vectors(f))
        Singular.push_term!(g, B(c), e)
    end
    return Singular.finish(g)
end


"""
    __Chow_ring__monomial_grading(M, f)

Returns the "pseudodegree" of the monomial `f` in the Chow ring of the moduli
space `M` passed.

This method is unsafe, as it does not consider the actual degree of the MPolyRingElem
objects passed. Instead, it assumes that the Chow ring passed has variables \$x_{i, j}\$
as in the Chow ring paper.
"""
function __Chow_ring_monomial_grading(M::QuiverModuliSpace, f)
    return __Chow_degrees(M.d)' * collect(Singular.exponent_vectors(f))[1]
end


"""
    __Chow_degrees(d)

Returns the vector of degrees for the variables of a Chow ring.

For internal use only.
"""
@memoize Dict function __Chow_degrees(d::AbstractVector{Int})
    return vcat([collect(1:di) for di in d if di > 0]...)
end


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

The dimension of the moduli space of the 3-Kronecker quiver:
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
                dimension_of_Luna_stratum(M, tau)
                for tau in all_Luna_types(M.Q, M.d, M.theta)
                    )
        end
    end
    # the semistable locus is also empty
    return "-∞"
end


"""
    is_smooth(M::QuiverModuliSpace)

Checks if the moduli space is smooth.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- whether the moduli space is smooth.

EXAMPLES:

Setups with `d` `theta`-coprime are smooth:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_smooth(M)
true
```
"""
function is_smooth(M::QuiverModuliSpace)
    if M.condition == "stable"
        return true
    elseif is_coprime(M.d, M.theta)
        return true
    elseif semistable_equals_stable(M)
        return true
    end

    throw(NotImplementedError("Not implemented for properly semistable cases."))
end

"""
    is_smooth(M::QuiverModuliStack)

Checks if the moduli stack is smooth.

This is always trus, as the quotient stack of a smooth variety is smooth.

INPUT:
- ``M``: a moduli stack of representations of a quiver.

OUTPUT:
- true

EXAMPLES:

```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliStack(Q, [2, 3]);

julia> is_smooth(M)
true
```
"""
function is_smooth(M::QuiverModuliStack)
    return true
end

"""
    is_projective(M::QuiverModuli)

Checks if the moduli space is projective.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:

- whether the moduli space is projective.

EXAMPLES:

The moduli space of the 3-Kronecker quiver is projective:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> is_projective(M)
true
```
"""
function is_projective(M::QuiverModuli)
    if is_acyclic(M.Q)
        if M.condition == "semistable"
            return true
        elseif M.condition == "stable"
            return semistable_equals_stable(M)
        end
    end

    SSP = semisimple_moduli_space(M)
    if M.condition == "semistable"
        return dimension(SSP) in [0, "-∞"]
    elseif M.condition == "stable"
        return dimension(SSP) in [1, "-∞"] && semistable_equals_stable(M)
    end
end


"""
    semisimple_moduli_space(M::QuiverModuli)

Returns the moduli space with the zero stability parameter.

Any quiver moduli space is (quasi)projective-over-affine;
this is the affine base.

INPUT:
- ``M``: a moduli space of representations of a quiver.

OUTPUT:
- the moduli space with the zero stability parameter.

EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(3); M = QuiverModuliSpace(Q, [2, 3]);

julia> dimension(semisimple_moduli_space(M))
0
```
"""
function semisimple_moduli_space(M::QuiverModuliSpace)
    return QuiverModuliSpace(M.Q, M.d, zero_vector(nvertices(M.Q)))
end