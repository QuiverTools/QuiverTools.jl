module QuiverTools

using StaticArrays

import Memoization
import AbstractAlgebra
import IterTools
import LinearAlgebraX
import Singular

import Base.show
import Memoization: @memoize
import AbstractAlgebra: fraction_field
import IterTools: subsets
import LinearAlgebraX: rankx
import Singular: polynomial_ring, degree, coeff, AlgebraHomomorphism, preimage, Ideal, QuotientRing, std, gens
import Combinatorics: with_replacement_combinations

export Quiver
export nvertices, narrows, indegree, outdegree, is_acyclic, is_connected, is_sink, is_source
export Euler_form, canonical_stability, is_coprime, slope
export is_Schur_root, generic_ext, generic_hom, canonical_decomposition, in_fundamental_domain
export all_HN_types, has_semistables, has_stables, is_amply_stable
export all_Teleman_bounds, all_weights_endomorphisms_universal_bundle, all_weights_universal_bundle, all_weights_irreducible_component_canonical
export Hodge_diamond, Hodge_polynomial, Picard_rank
export mKronecker_quiver, loop_quiver, subspace_quiver, three_vertex_quiver


"""
A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij})``,
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows ``i \\to j``.

Attributes:

- `adjacency` is the adjacency matrix of the quiver
- `name` is the name of the quiver, defaults to `""`.
"""
mutable struct Quiver
    adjacency::AbstractMatrix{Int}
    name::String

    function Quiver(adjacency::AbstractMatrix{Int}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else
            new(adjacency, name)
        end
    end
    function Quiver(adjacency::AbstractMatrix{Int})
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else
            new(adjacency, "")
        end
    end
end


function show(io::IO, Q::Quiver)
    if Q.name == ""
        print(io, "Quiver with adjacency matrix ")
    else
        print(io, Q.name * ", with adjacency matrix ")
    end
    print(io, Q.adjacency)
end


"""
Returns the (necessarily symmetric) adjacency matrix
of the underlying graph of the quiver.
"""
function underlying_graph(Q::Quiver)
    return Matrix{Int}(Q.adjacency + transpose(Q.adjacency) - diagonal(Q.adjacency))
end

"""
Returns the number of vertices of the quiver.
"""
nvertices(Q::Quiver) = size(Q.adjacency)[1]

"""
Returns the number of arrows of the quiver.
"""
narrows(Q::Quiver) = sum(Q.adjacency)

"""
Checks wether the quiver is acyclic, i.e. has no oriented cycles.
"""
is_acyclic(Q::Quiver) = all(entry == 0 for entry in Q.adjacency^nvertices(Q))

"""
Checks wether the underlying graph of the quiver is connected.

Examples:
```jldoctest
julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0]);

julia> is_connected(Q)
true

julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2]);

julia> is_connected(Q)
false

julia> # The 4-Kronecker quiver:

julia> Q = mKronecker_quiver(4);

julia> is_connected(Q)
true

julia> # The 4-loop quiver:

julia> Q = loop_quiver(4);

julia> is_connected(Q)
true

julia> # The 4-subspace quiver:

julia> Q = subspace_quiver(4);

julia> is_connected(Q)
true

julia> # The A10 quiver:

julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0] );

julia> is_connected(A10)
true

julia> # The A10 quiver without one arrow:

julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0] );

julia> is_connected(A10)
false
```
"""
function is_connected(Q::Quiver)
    paths = underlying_graph(Q)
    for i in 2:nvertices(Q)-1
        paths += paths * underlying_graph(Q)
    end
    for i in 1:nvertices(Q), j in 1:nvertices(Q)
        if i != j && paths[i, j] == 0 && paths[j, i] == 0
            return false
        end
    end
    return true
end

"""
Returns the number of incoming arrows to the vertex ``j``.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(4);

julia> indegree(Q, 1)
0

julia> indegree(Q, 2)
4
```
"""
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j])

"""
Returns the number of outgoing arrows from the vertex ``i``.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(4);

julia> outdegree(Q, 1)
4

julia> outdegree(Q, 2)
0
```
"""
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :])

"""
Checks if the vertex ``i`` is a source, i.e., a vertex with no incoming arrows.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(4);

julia> is_source(Q, 1)
true

julia> is_source(Q, 2)
false
```
"""
is_source(Q::Quiver, i::Int) = indegree(Q, i) == 0

"""
Checks if the vertex ``j`` is a sink, i.e., a vertex with no outgoing arrows.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(4);

julia> is_sink(Q, 1)
false

julia> is_sink(Q, 2)
true
```
"""
is_sink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

function arrows(Q::Quiver)
    n = nvertices(Q)
    return reduce(vcat, [[i, j] for k in 1:Q.adjacency[i, j]] for i in 1:n for j in 1:n if Q.adjacency[i, j] > 0)
end

identity_matrix(n::Int) = map(ind -> ind[1] == ind[2] ? 1 : 0, Iterators.product(1:n, 1:n))

function diagonal(m::AbstractMatrix{Int})
    n = size(m)[1]
    return map(ind -> ind[1] == ind[2] ? m[ind...] : 0, Iterators.product(1:n, 1:n))
end

"""
Returns the Euler matrix of the quiver.

The Euler matrix of a quiver ``Q`` is defined as 
```math
E = I - A,
```
where ``A`` is the adjacency matrix of ``Q`` and ``I`` is the identity matrix of the same size as ``A``.
"""
@memoize Dict Euler_matrix(Q::Quiver) = coerce_matrix(identity_matrix(nvertices(Q)) - Q.adjacency)

"""
Computes the Euler form of the quiver for vectors ``x`` and ``y``.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where ``E`` is the Euler matrix of the quiver.
"""
Euler_form(Q::Quiver, x::AbstractVector{Int}, y::AbstractVector{Int}) = x' * Euler_matrix(Q) * y

"""
The canonical stability parameter for the couple ``(Q,d)`` is given by ``<d,- > - < - ,d>``
"""
function canonical_stability(Q::Quiver, d::AbstractVector{Int})::AbstractVector{Int}
    return coerce_vector(-(-transpose(Euler_matrix(Q)) + Euler_matrix(Q)) * d)

"""Checks wether the given dimension vector ``d`` is ``\\theta``-coprime for
the stability parameter ``\\theta``."""
function is_coprime(d::AbstractVector{Int}, theta::AbstractVector{Int})
    return all(e -> theta' * e != 0, all_proper_subdimension_vectors(d))
end

"""Checks if the gcd of all the entries of d is ``1``."""
function is_coprime(d::AbstractVector{Int})
    return gcd(d) == 1
end

# TODO should this return a function?
"""
Returns the slope of the dimension vector ``d`` with respect to the stability parameter ``\\theta``
and a choice of a denominator function.
"""
function slope(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    return (theta' * d) // denom(d)
end
"""
Returns the subdimension vectors of ``d`` with a strictly larger slope than ``d``.
"""
@memoize Dict function all_destabilizing_subdimension_vectors(d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    return filter(e -> slope(e, theta, denom) > slope(d, theta, denom), all_nonzero_subdimension_vectors(d))
end

"""
Returns the list of all sequences ``(d^1,...,d^l)`` which sum to ``d`` such that ``\\mu(d^1) > ... > \\mu(d^l).``

Examples:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> QuiverTools.all_slope_decreasing_sequences(Q, d, theta)
8-element Vector{Vector{AbstractVector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]
```
"""
function all_slope_decreasing_sequences(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denominator::Function=sum)::Vector{Vector{AbstractVector{Int}}}
    coerce_vector!(d)
    coerce_vector!(theta)
    # List all subdimension vectors e of bigger slope than d.
    subdimensions = filter(e -> slope(e, theta, denominator) > slope(d, theta, denominator), all_nonzero_subdimension_vectors(d))

    # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    subdimensions = sort(subdimensions, by=e -> slope(e, theta, denominator))
    # The slope decreasing sequences which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

    # I will rewrite this as functional programming later

    function slope_filter(e, fstar)
        return slope(e, theta, denominator) > slope(fstar[1], theta, denominator)
    end

    function subdimensions_filter(e)
        return filter(fstar -> slope_filter(e, fstar), all_slope_decreasing_sequences(Q, d - e, theta, denominator))
    end

    allSlopeDecreasing = vcat(
        map(e -> map(fstar -> [e, fstar...], subdimensions_filter(e)), subdimensions)...,
    )

    # allSlopeDecreasing = []
    # for e in subdimensions
    #     for fstar in filter(fstar -> slope_filter(e, fstar), all_slope_decreasing_sequences(Q, d - e, theta, denominator))
    #         push!(allSlopeDecreasing, [e, fstar...])
    #     end
    # end
    # Add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
    return [[coerce_vector(d)], allSlopeDecreasing...]
end


"""Checks if there is a ``\\theta``-semistable representation of dimension vector ``d``.

Examples:
```jldoctest
julia> A2 = mKronecker_quiver(1); theta = [1,-1];

julia> has_semistables(A2, [1,1], theta)
true

julia> has_semistables(A2, [2,2], theta)
true

julia> has_semistables(A2, [1,2], theta)
false

julia> has_semistables(A2, [0,0], theta)
true

julia> # The 3-Kronecker quiver:

julia> K3 = mKronecker_quiver(3); theta = [3,-2];

julia> has_semistables(K3, [2,3], theta)
true

julia> has_semistables(K3, [1,4], theta)
false
```
"""
@memoize Dict function has_semistables(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, denom)
        subdimensionsBiggerSlope = filter(e -> slope(e, theta, denom) > slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensionsBiggerSlope)
    end
end

"""Checks if Q has a ``theta``-stable representation of dimension vector ``d``.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> has_stables(Q, d, theta)
true

julia> Q = mKronecker_quiver(2); d = [2,2]; theta = [1,-1];

julia> has_stables(Q, d, theta)
false

julia> has_semistables(Q, d, theta)
true
```
"""
@memoize Dict function has_stables(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, denom)
        subdimensions_bigger_or_equal_slope = filter(e -> slope(e, theta, denom) >= slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensions_bigger_or_equal_slope)
    end
end

"""
Checks if ``d`` is a Schur root for ``Q``.

By a lemma of Schofield (See Lemma 4.2 of [arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147)),
this is equivalent to the existence of a stable representation of dimension vector ``d``
for the canonical stability parameter.
	
Examples:

```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3];

julia> is_Schur_root(Q, d)
true
```
"""
is_Schur_root(Q::Quiver, d::AbstractVector{Int}) = has_stables(Q, d, canonical_stability(Q, d))

function is_real_root(Q, d)
    return Euler_form(Q, d, d) == 1
end

function is_imaginary_root(Q, d)
    return Euler_form(Q, d, d) <= 0
end

function is_isotropic_root(Q, d)
    return Euler_form(Q, d, d) == 0
end


"""Checks if ``e`` is a generic subdimension vector of ``d``.

A dimension vector ``e`` is called a generic subdimension vector of ``d`` if a generic representation
of dimension vector ``d`` possesses a subrepresentation of dimension vector ``e``.

By a result of Schofield (see Thm. 5.3 of [arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147)),
``e`` is a generic subdimension vector of ``d`` if and only if
```math
<e',d-e> \\geq 0
```
for all generic subdimension vectors ``e'`` of ``e``.
"""
@memoize Dict function is_generic_subdimension_vector(Q::Quiver, e::AbstractVector{Int}, d::AbstractVector{Int})
    if e == d || all(ei == 0 for ei in e)
        return true
    end
    # considering subdimension vectors that violate the numerical condition
    Euler_matrix_temp = Euler_matrix(Q) * (d - e) #to speed up computation of <eprime,d-e>
    subdimensions = filter(eprime -> eprime' * Euler_matrix_temp < 0, all_subdimension_vectors(e))
    # none of the subdimension vectors violating the condition should be generic
    return !any(eprime -> is_generic_subdimension_vector(Q, eprime, e), subdimensions)
end

"""
Returns the list of all generic subdimension vectors of ``d``.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3];

julia> QuiverTools.all_generic_subdimension_vectors(Q, d)
7-element Vector{AbstractVector{Int64}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 2]
 [0, 3]
 [1, 3]
 [2, 3]
```
"""
function all_generic_subdimension_vectors(Q::Quiver, d::AbstractVector{Int})::Vector{AbstractVector{Int}}
    return filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end

"""
Returns a list of all the Harder Narasimhan types of representations of ``Q``
with dimension vector ``d``, with respect to the slope function theta/denom.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_HN_types(Q, d, theta)
8-element Vector{Vector{AbstractVector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_HN_types(Q, [3,0], [0,0])
1-element Vector{Vector{AbstractVector{Int64}}}:
 [[3, 0]]

julia> Q = three_vertex_quiver(1,4,1); d = [4,1,4];

julia> theta = canonical_stability(Q, d);

julia> all_HN_types(Q, d, theta)
106-element Vector{Vector{AbstractVector{Int64}}}:
 [[4, 1, 4]]
 [[4, 1, 3], [0, 0, 1]]
 [[4, 0, 3], [0, 1, 1]]
 [[4, 0, 3], [0, 1, 0], [0, 0, 1]]
 [[3, 1, 2], [1, 0, 2]]
 [[3, 1, 2], [1, 0, 1], [0, 0, 1]]
 [[3, 0, 2], [1, 1, 2]]
 [[3, 0, 2], [0, 1, 0], [1, 0, 2]]
 [[3, 0, 2], [1, 0, 1], [0, 1, 1]]
 [[3, 0, 2], [1, 1, 1], [0, 0, 1]]
 ⋮
 [[3, 0, 0], [1, 1, 2], [0, 0, 2]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 4]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 3], [0, 0, 1]]
 [[3, 0, 0], [0, 1, 0], [1, 0, 2], [0, 0, 2]]
 [[3, 0, 0], [1, 0, 1], [0, 1, 1], [0, 0, 2]]
 [[3, 0, 0], [1, 1, 1], [0, 0, 3]]
 [[3, 0, 0], [1, 1, 0], [0, 0, 4]]
 [[4, 0, 0], [0, 1, 1], [0, 0, 3]]
 [[4, 0, 0], [0, 1, 0], [0, 0, 4]]
```
"""
@memoize Dict function all_HN_types(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum; ordered=true)

    coerce_vector!(d)
    if all(di == 0 for di in d)
        return [[d]]
    end
    # We consider just proper subdimension vectors which admit a semistable representation and for which μ(e) > μ(d)
    # Note that we also eliminate d by the following
    subdimensions = filter(e -> has_semistables(Q, e, theta, denom), all_destabilizing_subdimension_vectors(d, theta, denom))

    # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    if ordered
        subdimensions = sort(subdimensions, by=e -> slope(e, theta, denom))
    end

    # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN type of f = d-e and μ(e) > μ(f^1) holds.

    alltypes =
        Vector{AbstractVector{Int}}[
            [e, efstar...] for e in subdimensions
            for efstar in filter(
                fstar -> slope(e, theta, denom) > slope(fstar[1], theta, denom),
                all_HN_types(Q, d - e, theta, denom, ordered=ordered)
            )
        ]

    # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
    if has_semistables(Q, d, theta, denom)
        pushfirst!(alltypes, [d])
    end
    return alltypes
end

"""
Returns the codimension of the given HN stratum.

Examples:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> HN = all_HN_types(Q, d, theta)
8-element Vector{Vector{AbstractVector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> [QuiverTools.codimension_HN_stratum(Q, stratum) for stratum in HN]
8-element Vector{Int64}:
  0
  3
  4
 10
  8
  9
 12
 18
```
"""
function codimension_HN_stratum(Q::Quiver, stratum::AbstractVector{AbstractVector{Int}})
    if length(stratum) == 1
        return 0
    else
        return -sum(Euler_form(Q, stratum[i], stratum[j]) for i in 1:length(stratum)-1 for j in i+1:length(stratum))
    end
end

"""
Checks wether the dimension vector ``d`` is amply stable with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus in the parameter space is at least ``2``.
"""
function is_amply_stable(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    # We say that representations of a given dimension vector d are amply stable (for any notion of stability) if the codimension of the semistable locus is at least 2.
    # We verify this by computing the codimension of each HN stratum.
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
    return all(stratum -> codimension_HN_stratum(Q, stratum) >= 2, HN)
end


######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################


""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type."""
function Teleman_bound_onstratum(Q::Quiver, hntype::AbstractVector{AbstractVector{Int}}, theta::AbstractVector{Int}, denom::Function=sum)::Int
    if length(hntype) == 1
        throw(ArgumentError("Weight not defined for HN type of length 1."))
    end
    slopes = map(h -> slope(h, theta, denom), hntype)
    slopes = lcm(denominator.(slopes)) .* slopes
    return sum((slopes[t] - slopes[s]) * Euler_form(Q, hntype[s], hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype))
end

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given ``Q``, ``d``, ``\\theta`` and `denom``."""
function all_Teleman_bounds(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
    return Dict([hntype, Teleman_bound_onstratum(Q, hntype, theta, denom)] for hntype in HN)
end

"""Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type."""
function weights_universal_bundle_onstratum(theta::AbstractVector{Int}, a::AbstractVector{Int}, hntype, denom::Function=sum)::AbstractVector{Int}
    slopes = map(h -> slope(h, theta, denom), hntype)
    slopes *= lcm(denominator.(slopes))

    constant_term = sum(slopes[i] * (a' * hntype[i]) for i in eachindex(hntype))

    return -constant_term .+ slopes
end
"""Computes the weights of the universal bundle ``U_i(a)`` for the linearization ``a``
on all the non-dense Harder-Narasimhan strata for each 1-PS corresponding to each HN type."""
function all_weights_universal_bundle(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, a::AbstractVector{Int}, denom::Function=sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
    return Dict([hntype, weights_universal_bundle_onstratum(theta, a, hntype, denom)] for hntype in HN)
end


"""Computes the weight of the irreducible component of ``\\omega_R|_Z``
on a Harder-Narasimhan stratum for the 1-PS corresponding to each HN type.
More explicitly, if ``\\omega_X = \\mathcal{O}(rH)``, this returns the weight of
the pullback of O(H) on the given stratum."""
function weight_irreducible_component_canonical_on_stratum(Q::Quiver, d::AbstractVector{Int}, hntype::AbstractVector{AbstractVector{Int}}, theta::AbstractVector{Int}, denom::Function=sum)::Int
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

"""Computes the weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of ``\\mathcal{O}(H)`` on each stratum."""
function all_weights_irreducible_component_canonical(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta))
    return Dict([hntype, weight_irreducible_component_canonical_on_stratum(Q, d, hntype, theta, denom)] for hntype in HN)
end

"""Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type."""
function weights_endomorphism_universal_bundle_on_stratum(hntype::AbstractVector{AbstractVector{Int}}, theta::AbstractVector{Int}, denom::Function=sum)::AbstractVector{Int}
    # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is slope of first term in the HN type - slope of the last term in the HN type
    kweights = map(di -> slope(di, theta, denom), hntype)
    kweights = kweights * lcm(denominator.(kweights))
    # return kweights[1] - kweights[end] # this is the largest one
    return [kweights[i] - kweights[j] for i in 1:length(hntype) for j in 1:length(hntype)]
end

"""Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type."""
function all_weights_endomorphisms_universal_bundle(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, denom::Function=sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
    return Dict([hntype, weights_endomorphism_universal_bundle_on_stratum(hntype, theta, denom)] for hntype in HN)
end


#####################################################
# Canonical decomposition
#####################################################

# TODO check if the paper of Schofield is open access
"""
Computes the dimension of the ``\\mathrm{Ext}^1`` group between generic representations
of dimension vectors ``a`` and ``b``.

According to Thm. 5.4 in Schofield's [*General representations of quivers*](https://doi.org/10.1112/plms/s3-65.1.46),
we have

```math
ext(a,b) = max\\{- \\langle c, b\\rangle~~|~~c~\\text{is a generic subdimension vector of } a\\}.
```
"""
function generic_ext(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
    return maximum(-Euler_form(Q, c, b) for c in all_generic_subdimension_vectors(Q, a))
end

"""
Computes the dimension of the ``\\mathrm{Hom}`` group between generic representations
of dimension vectors ``a`` and ``b``.
"""
function generic_hom(Q::Quiver, a::AbstractVector{Int}, b::AbstractVector{Int})
    return Euler_form(Q, a, b) + generic_ext(Q, a, b)
end


# TODO add references
# TODO implement Derksen-Weyman?
# TODO add examples
"""
Computes the canonical decomposition of the dimension vector ``d`` for the given quiver ``Q``.

If ``\\beta_1, \\dots, \\beta_{\\ell}`` is a sequence
of Schur roots such that, for all ``i \\neq j``, one has

```math
\\mathrm{ext}(\\beta_i, \\beta_j) = \\mathrm{ext}(\\beta_j, \\beta_i) = 0,
```

then the general representation of dimension ``\\sum_i \\beta_i`` is
isomorphic to the direct sum of irreducible representations
of dimension vectors ``\\beta_i``.

Such a decomposition is called the canonical decomposition.
"""
function canonical_decomposition(Q::Quiver, d::AbstractVector{Int})
    # if is_Schur_root(Q, d)
    #     return [d]
    # end
    generic_subdimensions = filter(e -> e != d, all_generic_subdimension_vectors(Q, d))
    for e in generic_subdimensions
        if d - e in generic_subdimensions && generic_ext(Q, e, d - e) == 0 && generic_ext(Q, d - e, e) == 0
            return vcat(canonical_decomposition(Q, e), canonical_decomposition(Q, d - e))
        end
    end
    return [d] # if nothing above worked then d is a Schur root.
end

# """
# Checks wether the stability parameter theta is on a wall with respect to the wall-and-chamber decomposition for the dimension vector d.
# The wall and chamber decomposition is described in Section 2.2, MR4352662
# """
# function is_on_a_fake_wall(d::AbstractVector{Int}, theta::AbstractVector{Int}) 
#     return any(e -> e'*theta == 0, all_proper_subdimension_vectors(d))
# end

# TODO this is not relevant for anyone who is not me right now.
"""
Given an HN type ``(d^1,...,d^l)`` for the quiver Q,
returns the upper triangular matrix whose ``i,j``-th entry is ``\\mathrm{ext}(d^i,d^j)``.

The sum of all the entries is the codimension of the HN stratum;
the sum of all the rectangles starting on the "up-diagonal" (where the 1s go in a Jordan form)
and going all the way to the entry ``1,\\ell`` is at least ``1``.

The Teleman inequality is satisfied for this stratum
iif one of these rectangles sums to ``2`` or more.
"""
function extension_matrix(Q::Quiver, hntype::AbstractVector{AbstractVector{Int}})
    n = length(hntype)

    if n <= 1
        throw(ArgumentError("HN type must have length at least 2, this makes no sense for the dense stratum"))
    else
        M = zeros(Int, n, n)
        for i in 1:n-1, j in i+1:n
            M[i, j] = -Euler_form(Q, hntype[i], hntype[j])
        end
        return M
    end
end

# This can be modified once we decide how to properly represent a Luna type. Dict? Vector? Set? custom type?
function Luna_type_from_vector(vec::Vector{AbstractVector{Int}})
    Luna_type = Dict{AbstractVector{Int},Int}()
    for entry in vec
        if haskey(Luna_type, entry)
            Luna_type[entry] += 1
        else
            Luna_type[entry] = 1
        end
    end
    return Luna_type
end

function all_Luna_types(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}; denom::Function=sum)
    same_slope = filter(e -> slope(e, theta, denom) == slope(d, theta, denom) && has_stables(Q, e, theta, denom), QuiverTools.all_nonzero_subdimension_vectors(d))
    Luna_types = []
    bound = sum(d) ÷ minimum(sum(e) for e in same_slope) # the highest possible amount of repetitions for a given stable dimension vector
    for i in 1:bound+1
        for tau in with_replacement_combinations(same_slope, i)
            if sum(tau) == d
                push!(Luna_types, Luna_type_from_vector(tau))
            end
        end
    end
    return Luna_types
end


"""
Checks if the dimension vector ``d`` is in the fundamental domain of the quiver ``Q``.

The fundamental domain is the cone of dimension vectors in ``\\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i,

```math
(s_i, d) := \\langle d, s_i\\rangle + \\langle s_i, d\\rangle  \\leq 0,
```

where ``s_i`` is the dimension vector with all entries set to ``0`` and the i-th
set to ``1``.
"""
function in_fundamental_domain(Q::Quiver, d::AbstractVector{Int}; interior::Bool=false)
    # https://arxiv.org/abs/2209.14791 uses a strict inequality,
    # while https://arxiv.org/abs/2310.15927 uses a non-strict.
    # here we set it to non-strict by default.

    simples = [unit_vector(nvertices(Q), i) for i in 1:nvertices(Q)]
    if interior
        return all(simple -> Euler_form(Q, d, simple) + Euler_form(Q, simple, d) < 0, simples)
    end
    return all(simple -> Euler_form(Q, d, simple) + Euler_form(Q, simple, d) <= 0, simples)
end

##############################################################
# walls and chambers bases for stability parameters

function all_Schurian_decompositions(Q::Quiver, d::AbstractVector{Int})
    # the end of the recursion is necessarily at a simple root
    # TODO multisets instead of lists.
    if all(di == 0 for di in d)
        return [[]]
    elseif sum(d) == 1
        return [[d]]
    end
    Schur_subroots = filter(e -> is_Schur_root(Q, e), all_nonzero_subdimension_vectors(d))
    # @info d, Schur_subroots

    out = []
    for e in Schur_subroots, fstar in all_Schurian_decompositions(Q, d - e)
        push!(out, [e, fstar...])
    end
    return out
end

function has_semistables(Q::Quiver, d::AbstractVector{Int})
    # if there are less summands than vertices, there will always be a stability parameter
    # which multiplies all the summands (and hence d) to zero.
    allow_stability = []

    allSchurian = all_Schurian_decompositions(Q, d)
    for candidate in allSchurian
        if length(candidate) < nvertices(Q)
            push!(allow_stability, candidate)
        else
            # if not, one has to check.
            if LinearAlgebraX.rankx(hcat(candidate...)) < nvertices(Q)
                push!(allow_stability, candidate)
            end
        end
    end
    return allow_stability
end

"""
Checks if the stability parameter ``\\theta`` belongs to the cone of parameters admitting
stable representations of dimension vector ``d``.
Assumes that the dimension vector ``d`` is Schurian (for now).
"""
function in_stable_cone(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, strict::Bool=false)
    if !is_Schur_root(Q, d)
        throw(ArgumentError("d is not Schurian"))
    end
    if strict
        return all(e -> theta' * e < 0, all_generic_subdimension_vectors(Q, d))
    end
    return all(e -> theta' * e <= 0, all_generic_subdimension_vectors(Q, d))
end

# how to find inner walls now? These are given by a finite subset of the special subdimension vectors of d.
# which ones? and how to find them?

######################################################################################
# Below lie methods to compute Hodge diamonds translated from the Hodge diamond cutter.
# In turn, these are based on M. Reineke's paper
# "The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli", 
# https://doi.org/10.1007/s00222-002-0273-4
######################################################################################


###################################################
# auxiliary functions for Hodge_polynomial() below

"""
Solve ``A\\cdot x = b`` for ``A`` upper triangular via back substitution
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
    I = filter(e -> slope(e, theta) > slope(d, theta), all_proper_subdimension_vectors(d))
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
Returns the Hodge polynomial of the moduli space of ``\\theta``-semistable
representations of ``Q`` with dimension vector ``d``.

The algorithm is based on [MR1974891](https://doi.org/10.1007/s00222-002-0273-4),
and the current implementation is translated from the [Hodge diamond cutter]
(https://zenodo.org/doi/10.5281/zenodo.3893509).
"""
function Hodge_polynomial(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})

    # safety checks
    if theta' * d == 0 && !is_coprime(d)
        throw(ArgumentError("d is not coprime"))
    elseif theta' * d != 0 && gcd(theta' * d, sum(d)) != 1
        throw(ArgumentError("d is not coprime in the sense of Definition 6.3 of MR1974891."))
    elseif !is_acyclic(Q)
        throw(ArgumentError("Q is not acyclic."))
    end

    R, q = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, "q") # writing ["q"] throws bugs. No idea why.
    F = AbstractAlgebra.fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?

    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])

    # @warn "result needs to be a polynomial, otherwise the moduli space is singular."
    solution = solve(T, one_at_the_end)[1] * (1 - v)
    if denominator(solution) != 1
        throw(DomainError("Moduli space is singular!"))
    end
    result = numerator(solution)

    S, (x, y) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x", "y"])
    return result(x * y)
end

"""
Returns the Hodge diamond of the moduli space of
``\\theta``-semistable representations of ``Q`` with dimension vector ``d``.
"""
function Hodge_diamond(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
    g = Hodge_polynomial(Q, d, theta)
    return map(ind -> coeff(g, [ind[1] - 1, ind[2] - 1]).num, Iterators.product(1:degree(g, 1)+1, 1:degree(g, 2)+1))
end

"""
Computes the Picard rank of the moduli space of
``\\theta``-semistable representations of ``Q`` with dimension vector ``d``.
"""
function Picard_rank(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
    # TODO If over the complex numbers this should be h^{1,1}, since the moduli space is rational.
    # TODO This should follow from the long exact sequence in cohomology given by the exponential short exact sequence.

    return coeff(Hodge_polynomial(Q, d, theta), 2).num
end

function _Hodge_polynomial_fast(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int})
    # unsafe, curate input!
    # this is about 2% faster than the above, and occupies about 2% less memory.

    R, q = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, "q")
    F = AbstractAlgebra.fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?

    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])

    result = numerator(solve(T, one_at_the_end)[1] * (1 - v))
    # return [coeff(result, i) for i in 0:degree(result)] # this is actually all we need for the Hodge diamond because the matrix is diagonal for quiver moduli
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
Returns the symmetric polynomial of degree ``degree`` in the variables ``vars``.
"""
function symmetric_polynomial(vars, degree::Int)
    return sum(prod(e) for e in IterTools.subsets(vars, degree))
end

function Chow_ring(Q::Quiver, d::AbstractVector{Int}, theta::AbstractVector{Int}, a::AbstractVector{Int})
    # TODO cover case d[i] = 0
    # safety checks
    if !is_coprime(d, theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif a' * d != 1
        throw(ArgumentError("a is not a linearization"))
    end

    varnames = ["x$i$j" for i in 1:nvertices(Q) for j in 1:d[i] if d[i] > 0]
    # R, vars = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, varnames)
    R, vars = Singular.polynomial_ring(Singular.QQ, varnames)
    function chi(i, j)
        return vars[sum(d[1:i-1])+j]
    end

    function base_for_ring(name="naive")
        if name == "naive"
            bounds = [0:(d[i]-nu) for i in 1:nvertices(Q) for nu in 1:d[i]]
            lambdas = Iterators.product(bounds...)

            build_elem(lambda) = prod(prod(chi(i, nu)^lambda[sum(d[1:i-1])+nu] for nu in 1:d[i]) for i in 1:nvertices(Q))

            return map(l -> build_elem(l), lambdas)
        else
            throw(ArgumentError("unknown base."))
        end
    end

    # build the permutation group W
    W = Iterators.product([AbstractAlgebra.SymmetricGroup(d[i]) for i in 1:nvertices(Q)]...)
    sign(w) = prod(AbstractAlgebra.sign(wi) for wi in w)

    permute(f, sigma) = f([chi(i, sigma[i][j]) for i in 1:nvertices(Q) for j in 1:d[i]]...)

    delta = prod(prod(chi(i, l) - chi(i, k) for k in 1:d[i]-1 for l in k+1:d[i]) for i in 1:nvertices(Q) if d[i] > 1)
    antisymmetrize(f) = sum(sign(w) * permute(f, w) for w in W) / delta

    function all_forbidden(Q, d, theta, denom::Function=sum)
        dest = all_destabilizing_subdimension_vectors(d, theta, denom)
        return filter(e -> !any(f -> partial_order(Q, f, e), filter(f -> f != e, dest)), dest)
    end

    forbidden_polynomials = [prod(prod((chi(j, s) - chi(i, r))^Q.adjacency[i, j] for r in 1:e[i], s in e[j]+1:d[j]) for j in 1:nvertices(Q), i in 1:nvertices(Q) if Q.adjacency[i, j] > 0 && e[i] > 0 && d[j] > 1) for e in all_forbidden(Q, d, theta)]

    varnames2 = ["x$i$j" for i in 1:nvertices(Q) for j in 1:d[i] if d[i] > 0]
    A, Avars = Singular.polynomial_ring(Singular.QQ, varnames2)

    function xs(i, j)
        return Avars[sum(d[1:i-1])+j]
    end

    targets = [[symmetric_polynomial([chi(i, j) for j in 1:d[i]], k) for k in 1:d[i]] for i in 1:nvertices(Q)]
    targets = reduce(vcat, targets)

    inclusion = AlgebraHomomorphism(A, R, targets)

    anti = [antisymmetrize(f * b) for f in forbidden_polynomials for b in base_for_ring()]
    tautological = [gens(preimage(inclusion, Ideal(R, g)))[1] for g in anti]
    linear = [sum(a[i] * xs(i, 1) for i in 1:nvertices(Q))]

    return QuotientRing(A, std(Ideal(A, [tautological; linear])))
end

# TODO todd class
# TODO point class
# TODO universal bundle class



"""
Computes the quiver ``\\hat{Q}`` defined in [arXiv:0706.4306](https://doi.org/10.48550/arXiv.0706.4306)
for the given quiver ``Q`` and a given subdimension vector ``e``.
"""
@memoize Dict function smooth_model_quiver(Q::Quiver, e::AbstractVector{Int})
    n = nvertices(Q)
    A = zeros(Int, n + 1, n + 1)
    A[2:n+1, 2:n+1] = Q.adjacency
    A[1, 2:n+1] = e
    return Quiver(A, "Smooth model quiver of " * Q.name * " for the dimension vector " * e)
end

# TODO think of a better way to do this?
@memoize Dict function smooth_model_root(d::AbstractVector{Int}; inclusion::Bool=false)
    if inclusion
        return coerce_vector([0, d...])
    end
    return coerce_vector([1, d...])
end



#################
# Technical tools
#################

@memoize Dict function zero_vector(n::Int)
    return coerce_vector((zeros(Int, n)))
end

function thin_dimension_vector(Q::Quiver)
    return coerce_vector(ones(Int, nvertices(Q)))
end

@memoize Dict function all_subdimension_vectors(d::AbstractVector{Int})::Array{AbstractVector{Int}}
    coerce_vector!(d)
    return coerce_vector.(collect(Iterators.product(map(di -> 0:di, d)...)))
end
@memoize Dict function all_nonzero_subdimension_vectors(d::AbstractVector{Int})
    return filter(e -> !all(ei == 0 for ei in e), all_subdimension_vectors(d))
end

@memoize Dict function all_proper_subdimension_vectors(d::AbstractVector{Int})
    return filter(e -> any(ei != 0 for ei in e) && e != d, all_subdimension_vectors(d))
end

@memoize Dict function unit_vector(n::Int, i::Int)
    v = zeros(Int, n)
    v[i] = 1
    return coerce_vector(v)
end

function unit_vector(Q::Quiver, i::Int)
    return unit_vector(nvertices(Q), i)
end


function coerce_vector(v::AbstractVector)
    return SVector{length(v)}(v)
end
function coerce_vector(v::Tuple)
    return SVector{length(v)}(v)
end

function coerce_vector!(x)
    x = coerce_vector(x)
end

function coerce_matrix(m::AbstractMatrix)
    return SMatrix{size(m)...}(m)
end


#######################################################
# Include all the submodules
#######################################################

include("constructors.jl")

######################
# end of QuiverTools
######################
end