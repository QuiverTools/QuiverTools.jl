module QuiverTools

using StaticArrays

using Memoization: Memoization
using IterTools: IterTools
using LinearAlgebraX: LinearAlgebraX
using Singular: Singular
using AbstractAlgebra: AbstractAlgebra
using Nemo: Nemo

import Base.show
import Memoization: @memoize
import IterTools: subsets
import LinearAlgebraX: rankx
import Singular:
    polynomial_ring,
    degree,
    coeff,
    constant_coefficient,
    AlgebraHomomorphism,
    preimage,
    Ideal,
    quotient_ideal,
    QuotientRing,
	fraction_field,
    std,
    gens,
    base_ring
import Combinatorics: with_replacement_combinations, partitions

export Quiver
export nvertices, narrows, arrows, indegree, outdegree,
    is_acyclic, is_connected, is_sink, is_source
export Euler_form, canonical_stability, is_coprime, slope
export is_Schur_root,
    generic_ext, generic_hom, canonical_decomposition, in_fundamental_domain
export all_HN_types, is_HN_type, has_semistables, has_stables, is_amply_stable
export all_Teleman_bounds,
    all_weights_endomorphisms_universal_bundle,
    all_weights_universal_bundle,
    all_weights_irreducible_component_canonical
export Hodge_diamond, Hodge_polynomial, Picard_rank
export QuiverModuli, is_nonempty, dimension, Poincare_polynomial, Betti_numbers,
    is_smooth, semistable_equals_stable, codimension_unstable_locus

# TODO add doctests to every method in this file.

"""
A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij})``,
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows ``i \\to j``.

Attributes:

- `adjacency` is the adjacency matrix of the quiver
- `name` is the name of the quiver, defaults to `""`.
"""
struct Quiver
    adjacency::AbstractMatrix{Int}
    name::String

    """
        Quiver(adjacency, name = "")

    Constructs a quiver starting from its adjacency matrix, and an optional name.

    EXAMPLES:
    ```jldoctest
    julia> m = [0 1; 2 0];

    julia> Quiver(m, "my quiver")
    my quiver, with adjacency matrix [0 1; 2 0]
    ```
    """
    function Quiver(adjacency::AbstractMatrix{Int}, name::String = "")
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else
            new(adjacency, name)
        end
    end

    """
        Quiver(arrows)

    Constructs a quiver based on its arrows encoded in a string.

    the string `arrows` must be of the form
    
    ```i---j,k-...-s```

    where `i`, `j` and all the vertices are positive integers.
    The amount of characters between `i` and  `j` is then the number of arrows i -> j.

    EXAMPLE:
    ```jldoctest
    julia> Q = Quiver("1--2,1---3,2----3")
    Quiver with adjacency matrix [0 2 3; 0 0 4; 0 0 0]
    ```
    """
    function Quiver(arrows::String)
        pairs = split(arrows, ",")
        pairs = map(p -> split(p, "-"), pairs)
        pairs = map( p -> [parse(Int, p[1]), parse(Int, p[end]), length(p) - 1], pairs)

        n = maximum(maximum(pair[1:2]) for pair in pairs)
        A = zeros(Int, n, n)
        for pair in pairs
            A[pair[1], pair[2]] = pair[3]
        end
        return Quiver(A)
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

function deglex_key(Q::Quiver, e::AbstractVector{Int})::Int
    b = maximum(e) + 1
    n = nvertices(Q)

    return (
            sum(e[i] * b^(n - i) for i in 1:length(e))
            + sum(e) * b^n
        )
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

EXAMPLES:
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

EXAMPLES:
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

EXAMPLES:
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

EXAMPLES:
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

EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(4);

julia> is_sink(Q, 1)
false

julia> is_sink(Q, 2)
true
```
"""
is_sink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

"""
    arrows(Q::Quiver)

Returns a list of all arrows of the quiver ``Q``.

INPUT:
- `Q`: a quiver

OUTPUT:
- a list of all arrows of the quiver ``Q``.

EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> arrows(Q)
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 2]
 [1, 2]
```
"""
function arrows(Q::Quiver)
    n = nvertices(Q)
    return reduce(
        vcat,
        [[i, j] for k in 1:Q.adjacency[i, j]] for i in 1:n for
        j in 1:n if Q.adjacency[i, j] > 0
    )
end


# this is the wheel reinvention department.
# I don't want to load the whole LinearAlgebra package just for this.
"""
Returns the identity matrix of size ``n``.
"""
identity_matrix(n::Int) = map(ind -> ind[1] == ind[2] ? 1 : 0, Iterators.product(1:n, 1:n))

function diagonal(m::AbstractMatrix{Int})
    n = size(m)[1]
    return map(ind -> ind[1] == ind[2] ? m[ind...] : 0, Iterators.product(1:n, 1:n))
end
function diagonal(v::AbstractVector)
    n = length(v)
    return map(ind -> ind[1] == ind[2] ? v[ind[1]] : 0, Iterators.product(1:n, 1:n))
end
"""
Returns the Euler matrix of the quiver.

The Euler matrix of a quiver ``Q`` is defined as 
```math
E = I - A,
```
where ``A`` is the adjacency matrix of ``Q`` and ``I``
is the identity matrix of the same size as ``A``.
"""
@memoize Dict Euler_matrix(Q::Quiver) =
    coerce_matrix(identity_matrix(nvertices(Q)) - Q.adjacency)

"""
Computes the Euler form of the quiver for vectors ``x`` and ``y``.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where ``E`` is the Euler matrix of the quiver.
"""
Euler_form(Q::Quiver, x::AbstractVector{Int}, y::AbstractVector{Int}) =
    x' * Euler_matrix(Q) * y

"""
The canonical stability parameter for the couple ``(Q,d)`` is given by ``<d,-> - <-,d>``
"""
function canonical_stability(Q::Quiver, d::AbstractVector{Int})::AbstractVector{Int}
    return coerce_vector(-(-transpose(Euler_matrix(Q)) + Euler_matrix(Q)) * d)
end

"""Checks wether the given dimension vector ``d`` is ``\\theta``-coprime for
the stability parameter ``\\theta``."""
function is_coprime(d::AbstractVector{Int}, theta::AbstractVector{Int})
    return all(
        e -> theta' * e != 0,
        all_subdimension_vectors(d, nonzero = true, strict = true),
    )
end

"""Checks if the gcd of all the entries of d is ``1``."""
function is_coprime(d::AbstractVector{Int})
    return gcd(d) == 1
end


"""
Returns the slope of the dimension vector ``d``
with respect to the stability parameter ``\\theta``
and a choice of a denominator function.
"""
function slope(d::AbstractVector{Int},
    theta::AbstractVector{Int},
    denom::Function = sum)

    return (theta' * d) // denom(d)
end
"""
Returns the subdimension vectors of ``d`` with a strictly larger slope than ``d``.
"""
@memoize Dict function all_destabilizing_subdimension_vectors(
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    denom::Function = sum,
)
    return filter(
        e -> slope(e, theta, denom) > slope(d, theta, denom),
        all_subdimension_vectors(d, nonzero = true),
    )
end

"""
Returns the list of all sequences ``(d^1,...,d^l)`` which sum to ``d``
such that ``\\mu(d^1) > ... > \\mu(d^l).``

EXAMPLES:
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
function all_slope_decreasing_sequences(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    denominator::Function = sum,
)::Vector{Vector{AbstractVector{Int}}}

    d = coerce_vector(d)
    theta = coerce_vector(theta)
    # List all subdimension vectors e of bigger slope than d.
    subdimensions = filter(
        e -> slope(e, theta, denominator) > slope(d, theta, denominator),
        all_subdimension_vectors(d, nonzero = true),
    )

    # We sort the subdimension vectors by slope because that will return the list of
    # all HN types in ascending order with respect to the partial order from
    # Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    subdimensions = sort(subdimensions, by = e -> slope(e, theta, denominator))
    # The slope decreasing sequences which are not of the form (d)
    # are given by (e,f^1,...,f^s) where e is a proper subdimension vector
    # such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of
    # f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

    function slope_filter(e, fstar)
        return slope(e, theta, denominator) > slope(fstar[1], theta, denominator)
    end

    function subdimensions_filter(e)
        return filter(
            fstar -> slope_filter(e, fstar),
            all_slope_decreasing_sequences(Q, d - e, theta, denominator),
        )
    end

    allSlopeDecreasing = vcat(
        map(e -> map(fstar -> [e, fstar...], subdimensions_filter(e)), subdimensions)...,
    )


    # Add d again, at the beginning, because it is smallest
    # with respect to the partial order from Def. 3.6
    return [[coerce_vector(d)], allSlopeDecreasing...]
end


"""Checks if there is a ``\\theta``-semistable representation of dimension vector ``d``.

EXAMPLES:
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
```
The 3-Kronecker quiver:

```jldoctest
julia> K3 = mKronecker_quiver(3); theta = [3,-2];

julia> has_semistables(K3, [2,3], theta)
true

julia> has_semistables(K3, [1,4], theta)
false
```
"""
@memoize Dict function has_semistables(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum,
)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, denom)
        subdimensionsBiggerSlope = filter(
            e -> slope(e, theta, denom) > slope_d,
            all_subdimension_vectors(d, nonzero = true, strict = true),
        )
        # to have semistable representations, none of the vectors above must be
        # a generic subdimension vector.
        return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensionsBiggerSlope)
    end
end

"""Checks if Q has a ``theta``-stable representation of dimension vector ``d``.

EXAMPLES:
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

The zero dimension vector has no stables:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [0,0];

julia> has_stables(Q, d)
false
```
"""
@memoize Dict function has_stables(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum,
)
    if all(di == 0 for di in d)
        return false
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, denom)
        subdimensions_bigger_or_equal_slope = filter(
            e -> slope(e, theta, denom) >= slope_d,
            all_subdimension_vectors(d, nonzero = true, strict = true),
        )
        # to have semistable representations,
        # none of the vectors above must be generic subdimension vectors.
        return all(
            e -> !is_generic_subdimension_vector(Q, e, d),
            subdimensions_bigger_or_equal_slope,
        )
    end
end

"""
Checks if ``d`` is a Schur root for ``Q``.

By [Lemma 4.2, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147),
this is equivalent to the existence of a stable representation of dimension vector ``d``
for the canonical stability parameter.
	
EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(3); d = [2,3];

julia> is_Schur_root(Q, d)
true
```
"""
is_Schur_root(Q::Quiver, d::AbstractVector{Int}) =
    has_stables(Q, d, canonical_stability(Q, d))

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

A dimension vector ``e`` is called a generic subdimension vector of ``d``
if a generic representation of dimension vector ``d`` possesses a subrepresentation
of dimension vector ``e``.

By [Theorem 5.3, arXiv:0802.2147](https://doi.org/10.48550/arXiv.0802.2147),
``e`` is a generic subdimension vector of ``d`` if and only if
```math
<e',d-e> \\geq 0
```
for all generic subdimension vectors ``e'`` of ``e``.
"""
@memoize Dict function is_generic_subdimension_vector(
    Q::Quiver,
    e::AbstractVector{Int},
    d::AbstractVector{Int},
)
    if e == d || all(ei == 0 for ei in e)
        return true
    end
    # considering subdimension vectors that violate the numerical condition
    Euler_matrix_temp = Euler_matrix(Q) * (d - e) #to speed up computation of <eprime,d-e>
    subdimensions =
        filter(eprime -> eprime' * Euler_matrix_temp < 0, all_subdimension_vectors(e))
    # none of the subdimension vectors violating the condition should be generic
    return !any(eprime -> is_generic_subdimension_vector(Q, eprime, e), subdimensions)
end

"""
Returns the list of all generic subdimension vectors of ``d``.

EXAMPLES:
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
function all_generic_subdimension_vectors(
    Q::Quiver,
    d::AbstractVector{Int},
    )::Vector{AbstractVector{Int}}

    return filter(e -> is_generic_subdimension_vector(Q, e, d),
                    all_subdimension_vectors(d)
                    )
end

"""
Returns a list of all the Harder Narasimhan types of representations of ``Q``
with dimension vector ``d``, with respect to the slope function theta/denom.

EXAMPLES:
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
@memoize Dict function all_HN_types(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    denom::Function = sum,
    ordered = true,
    )

    d = coerce_vector(d)
    if all(di == 0 for di in d)
        return [[d]]
    end
    # We consider just proper subdimension vectors which admit a semistable
    # representation and for which μ(e) > μ(d)
    # Note that we also eliminate d by the following
    subdimensions = filter(
        e -> has_semistables(Q, e, theta, denom),
        all_destabilizing_subdimension_vectors(d, theta, denom),
    )

    # We sort the subdimension vectors by slope because that will return the list of
    # all HN types in ascending order with respect to the partial order from
    # Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    if ordered
        subdimensions = sort(subdimensions, by = e -> slope(e, theta, denom))
    end

    # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a
    # proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN
    # type of f = d-e and μ(e) > μ(f^1) holds.

    alltypes = Vector{AbstractVector{Int}}[
        [e, efstar...] for e in subdimensions for efstar in filter(
            fstar -> slope(e, theta, denom) > slope(fstar[1], theta, denom),
            all_HN_types(Q, d - e, theta, denom, ordered),
        )
    ]

    # Possibly add d again, at the beginning, because it is smallest
    # with respect to the partial order from Def. 3.6
    if has_semistables(Q, d, theta, denom)
        pushfirst!(alltypes, [d])
    end
    return alltypes
end

"""
	is_hn_type(Q, d, dstar, theta, denom)

Checks if the given ordered list of subdimension vectors ``dstar`` is an HN type
for the datum ``(Q, d)`` and the slope stability given by ``(theta, denom)``.


EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver(3);

julia> d = [2, 3]; dstar = [d];

julia> is_HN_type(Q, d, dstar)
true

```
"""
function is_HN_type(
    Q::Quiver,
    d::AbstractVector{Int},
    dstar::Vector{<:AbstractVector{Int}},
    theta::AbstractVector{Int} = canonical_stability(Q, d),
    denom::Function = sum,
)::Bool
    if sum(dstar) != d
        return false
    end

    if !all(
        slope(dstar[i], theta, denom) > slope(dstar[i+1], theta, denom) for
        i = 1:length(dstar)-1
    )
        return false
    end

    if !all(has_semistables(Q, dstari, theta, denom) for dstari in dstar)
        return false
    end
    return true
end



"""
Returns the codimension of the given HN stratum.

EXAMPLES:
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
function codimension_HN_stratum(Q::Quiver, stratum::Vector{<:AbstractVector{Int}})
    if length(stratum) == 1
        return 0
    else
        return -sum(
            Euler_form(Q, stratum[i], stratum[j]) for i in 1:length(stratum)-1 for
            j = i+1:length(stratum)
        )
    end
end

"""
Checks wether the dimension vector ``d`` is amply stable
with respect to the slope function `theta`/`denominator`.

This means that the codimension of the unstable locus
in the parameter space is at least ``2``.
"""
function is_amply_stable(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    denom::Function = sum,
)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, denom))
    return all(stratum -> codimension_HN_stratum(Q, stratum) >= 2, HN)
end

#####################################################
# Canonical decomposition
#####################################################

"""
Computes the dimension of the ``\\mathrm{Ext}^1`` group between generic representations
of dimension vectors ``a`` and ``b``.

According to [Theorem 5.4, MR1162487]
(https://mathscinet.ams.org/mathscinet/relay-station?mr=1162487),
we have

```math
ext(a,b)=max\\{-\\langle c,b\\rangle~~|~~c~\\text{is a generic subdimension vector of }a\\}.
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


# TODO add examples
"""
Computes the canonical decomposition of the dimension vector ``d``
for the given quiver ``Q``.

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
        if d - e in generic_subdimensions &&
           generic_ext(Q, e, d - e) == 0 &&
           generic_ext(Q, d - e, e) == 0
            return vcat(canonical_decomposition(Q, e), canonical_decomposition(Q, d - e))
        end
    end
    return [d] # if nothing above worked then d is a Schur root.
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
function in_fundamental_domain(Q::Quiver, d::AbstractVector{Int}; interior::Bool = false)
    # https://arxiv.org/abs/2209.14791 uses a strict inequality,
    # while https://arxiv.org/abs/2310.15927 uses a non-strict.
    # here we set it to non-strict by default.

    simples = [unit_vector(nvertices(Q), i) for i in 1:nvertices(Q)]
    if interior
        return all(
            simple -> Euler_form(Q, d, simple) + Euler_form(Q, simple, d) < 0,
            simples,
        )
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
    Schur_subroots =
        filter(e -> is_Schur_root(Q, e), all_subdimension_vectors(d, nonzero = true))
    # @info d, Schur_subroots

    out = []
    for e in Schur_subroots, fstar in all_Schurian_decompositions(Q, d - e)
        push!(out, [e, fstar...])
    end
    return out
end

function admits_semistables(Q::Quiver, d::AbstractVector{Int})
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
Checks if the stability parameter ``\\theta`` belongs to the cone of parameters
admitting stable representations of dimension vector ``d``.

Assumes that the dimension vector ``d`` is Schurian (for now).
"""
function in_stable_cone(
    Q::Quiver,
    d::AbstractVector{Int},
    theta::AbstractVector{Int},
    strict::Bool = false,
)
    if !is_Schur_root(Q, d)
        throw(ArgumentError("d is not Schurian"))
    end
    if strict
        return all(e -> theta' * e < 0, all_generic_subdimension_vectors(Q, d))
    end
    return all(e -> theta' * e <= 0, all_generic_subdimension_vectors(Q, d))
end


#################
# Technical tools
#################

"""
	zero_vector(n::Int)

Create a zero vector of length `n`.

INPUT:
- `n::Int`: The length of the zero vector.

OUTPUT:
- A zero vector of length `n`.
"""
@memoize Dict function zero_vector(n::Int)
    return coerce_vector(zeros(Int, n))
end

"""
	thin_dimension_vector(Q::Quiver)

Compute the thin dimension vector for a given quiver `Q`.

INPUT:
- `Q::Quiver`: The input quiver.

OUTPUT:
- A vector of ones of length `n`.

"""
function thin_dimension_vector(Q::Quiver)
    return coerce_vector(ones(Int, nvertices(Q)))
end

"""
	all_subdimension_vectors(d::AbstractVector{Int})

Compute all subdimension vectors of a given dimension vector `d`.

INPUT:
- `d::AbstractVector{Int}`: The input dimension vector.
- `nonzero::Bool=false`: wether to exclude the zero vector.
- `strict::Bool=false`: wether to exclude the input vector `d`.

OUTPUT:
- An array of all subdimension vectors of `d`, with or without the zero vector and `d`.

Examples
```jldoctest
julia> all_subdimension_vectors([2, 3])
12-element Vector{AbstractVector{Int64}}:
 [0, 0]
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_subdimension_vectors([2, 3], nonzero=true)
11-element Vector{AbstractVector{Int64}}:
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]
 [2, 3]

julia> QuiverTools.all_subdimension_vectors([2, 3], nonzero=true, strict=true)
10-element Vector{AbstractVector{Int64}}:
 [1, 0]
 [2, 0]
 [0, 1]
 [1, 1]
 [2, 1]
 [0, 2]
 [1, 2]
 [2, 2]
 [0, 3]
 [1, 3]
 ```
"""
@memoize Dict function all_subdimension_vectors(
    d::AbstractVector{Int};
    nonzero::Bool = false,
    strict::Bool = false,
    )::Array{AbstractVector{Int}}

    d = coerce_vector(d)

    subdims = coerce_vector.(collect(Iterators.product(map(di -> 0:di, d)...)))
    subdims = filter(e -> true, subdims) #really now
    if nonzero
        subdims = filter(e -> any(ei != 0 for ei in e), subdims)
    end
    if strict
        subdims = filter(e -> e != d, subdims)
    end
    return subdims
end

function is_subdimension_vector(e::AbstractVector{Int}, d::AbstractVector{Int})
    return all(ei <= di for (ei, di) in zip(e, d))
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

function coerce_matrix(m::AbstractMatrix)
    return SMatrix{size(m)...}(m)
end


#######################################################
# Include all the submodules
#######################################################

include("constructors.jl")
include("moduli.jl")
include("teleman.jl")

######################
# end of QuiverTools
######################
end
