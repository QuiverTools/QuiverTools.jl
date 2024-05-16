module QuiverTools

import Base.show
using Memoize, LinearAlgebra
using Nemo

export Quiver
export nvertices, narrows, indegree, outdegree, is_acyclic, is_connected, is_sink, is_source
export Euler_form, canonical_stability, is_coprime, slope
export is_Schur_root, generic_ext, generic_hom, canonical_decomposition
export all_HN_types, has_semistables, has_stables, is_amply_stable
export all_Teleman_bounds, all_weights_endomorphisms_universal_bundle, all_weights_universal_bundle, all_weights_irreducible_component_canonical
export Hodge_diamond, Hodge_polynomial, Picard_rank
export mKronecker_quiver, loop_quiver, subspace_quiver, three_vertex_quiver


# TODO coframed_quiver()
# TODO full_subquiver()
# TODO bipartite quiver
# TODO wall and chamber decomposition

"""
A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij})``,
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows i → j.

Attributes:

- `adjacency` is the adjacency matrix of the quiver
- `name` is the name of the quiver, defaults to `""`.
"""
mutable struct Quiver
    adjacency::Matrix{Int}
    name::String

    function Quiver(adjacency::Matrix{Int}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else 
            new(adjacency, name)
        end
    end 
    function Quiver(adjacency::Matrix{Int})
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
        print(io, Q.name*", with adjacency matrix ")
    end
    println(io, Q.adjacency)
end


"""
Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.
"""
function underlying_graph(Q::Quiver)
    return Matrix{Int}(Q.adjacency + transpose(Q.adjacency) - diagm(diag(Q.adjacency)))
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
```julia-repl
julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0])

julia> is_connected(Q)
true

julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2])
false

# The 4-Kronecker quiver:
julia> Q = mKroneckerquiver(4)

julia> is_connected(Q)
true

# The 4-loop quiver:
julia> Q = LoopQuiver(4)

julia> is_connected(Q)
true

# The 4-subspace quiver:
julia> Q = SubspaceQuiver(4)

julia> is_connected(Q)
true

# The A10 quiver:
julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                        0 0 1 0 0 0 0 0 0 0;
                        0 0 0 1 0 0 0 0 0 0;
                        0 0 0 0 1 0 0 0 0 0;
                        0 0 0 0 0 1 0 0 0 0;
                        0 0 0 0 0 0 1 0 0 0;
                        0 0 0 0 0 0 0 1 0 0;
                        0 0 0 0 0 0 0 0 1 0;
                        0 0 0 0 0 0 0 0 0 1;
                        0 0 0 0 0 0 0 0 0 0] )

julia> is_connected(A10)
true

# The A10 quiver without one arrow:
julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                        0 0 1 0 0 0 0 0 0 0;
                        0 0 0 1 0 0 0 0 0 0;
                        0 0 0 0 1 0 0 0 0 0;
                        0 0 0 0 0 1 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 1 0 0;
                        0 0 0 0 0 0 0 0 1 0;
                        0 0 0 0 0 0 0 0 0 1;
                        0 0 0 0 0 0 0 0 0 0] )

julia> is_connected(A10)
false
```
"""
function is_connected(Q::Quiver)
    paths = underlying_graph(Q)
    for i in 2:nvertices(Q) - 1
        paths += paths*underlying_graph(Q)
    end
    for i in 1:nvertices(Q), j in 1:nvertices(Q)
            if i != j && paths[i, j] == 0 && paths[j, i] == 0
                return false
            end
    end
    return true
end

# the docstrings on these functions are from the file quiver.py

"""
Returns the number of incoming arrows to the vertex j.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)

julia> indegree(Q, 1)
0

julia> indegree(Q, 2)
4
```
"""
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j])

"""
Returns the number of outgoing arrows from the vertex i.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)

julia> outdegree(Q, 1)
4

julia> outdegree(Q, 2)
0
```
"""
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :])

"""
Checks if the vertex i is a source, i.e., a vertex with no incoming arrows.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)

julia> is_source(Q, 1)
true

julia> is_source(Q, 2)
false
```
"""
is_source(Q::Quiver, i::Int) = indegree(Q, i) == 0

"""
Checks if the vertex j is a sink, i.e., a vertex with no outgoing arrows.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)

julia> is_sink(Q, 1)
false

julia> is_sink(Q, 2)
true
```
"""
is_sink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

function arrows(Q::Quiver)
    return reduce(vcat, [[i,j] for k in 1:Q.adjacency[i,j]] for i in 1:nvertices(Q) for j in 1:nvertices(Q) if Q.adjacency[i,j] > 0)
end

"""
Returns the Euler matrix of the quiver.

The Euler matrix of a quiver Q is defined as 
```math
E = I - A,
```
where ``A`` is the adjacency matrix of Q and ``I`` is the identity matrix of the same size as ``A``.
"""
@memoize Dict Euler_matrix(Q::Quiver) = Matrix{Int}(LinearAlgebra.I, nvertices(Q), nvertices(Q)) - Q.adjacency

"""
Computes the Euler form of the quiver for vectors x and y.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where E is the Euler matrix of the quiver.
"""
Euler_form(Q::Quiver, x::Vector{Int}, y::Vector{Int}) = x'*Euler_matrix(Q)*y

"""
The canonical stability parameter for the couple ``(Q,d)`` is given by ``<d,- > - < - ,d>``
"""
canonical_stability(Q::Quiver, d::Vector{Int})::Vector{Int} = -(-transpose(Euler_matrix(Q)) + Euler_matrix(Q))*d

"""Checks wether the given dimension vector ``d`` is ``theta``-coprime for the stability parameter ``theta``."""
function is_coprime(d::Vector{Int}, theta::Vector{Int})
    return all(e -> theta'*e != 0, all_proper_subdimension_vectors(d))
end

"""Checks if the gcd of all the entries of d is ``1``."""
function is_coprime(d::Vector{Int})
    return gcd(d) == 1
end

# TODO should this return a function?
"""
Returns the slope of the dimension vector ``d`` with respect to the stability parameter ``theta``
and a choice of a denominator function.
"""
function slope(d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    return (theta'*d)//slope_denominator(d)
end
"""
Returns the subdimension vectors of ``d`` with a strictly larger slope than ``d``.
"""
@memoize Dict function all_destabilizing_subdimension_vectors(d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    return filter(e -> slope(e, theta, slope_denominator) > slope(d, theta, slope_denominator), all_nonzero_subdimension_vectors(d))
end

"""
Returns the list of all sequences ``(d^1,...,d^l)`` which sum to d such that ``\\mu(d^1) > ... > \\mu(d^l).``

Examples:
```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_slope_decreasing_sequences(Q, d, theta)
8-element Array{Array{Vector}}:
    [[[2, 3]],
    [[1, 1], [1, 2]],
    [[2, 2], [0, 1]],
    [[2, 1], [0, 2]],
    [[1, 0], [1, 3]],
    [[1, 0], [1, 2], [0, 1]],
    [[1, 0], [1, 1], [0, 2]],
    [[2, 0], [0, 3]]]
```
"""
function all_slope_decreasing_sequences(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, denominator::Function = sum)

    # List all subdimension vectors e of bigger slope than d.
    subdimensions = filter(e -> slope(e,theta,denominator) > slope(d,theta,denominator), all_nonzero_subdimension_vectors(d))

    # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    subdimensions = sort(subdimensions, by = e -> slope(e,theta,denominator))
    # The slope decreasing sequences which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

    # I will rewrite this as functional programming later
    allSlopeDecreasing = []
    for e in subdimensions
        for fstar in filter(fstar -> slope(e,theta,denominator) > slope(fstar[1],theta, denominator), all_slope_decreasing_sequences(Q, d-e, theta, denominator))
        push!(allSlopeDecreasing, [e, fstar...])
        end
    end
    # Add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
    return [[d], allSlopeDecreasing...]
end


"""Checks if there is a theta-semistable representation of dimension vector d.

Examples:
```julia-repl
julia> A2 = mKroneckerquiver(1); theta = [1,-1];

julia> has_semistables(A2, [1,1], theta)
true

julia> has_semistables(A2, [2,2], theta)
true

julia> has_semistables(A2, [1,2], theta)
false

julia> has_semistables(A2, [0,0], theta)
true

# The 3-Kronecker quiver:
julia> K3 = mKroneckerquiver(3); theta = [3,-2];

julia> has_semistables(K3, [2,3], theta)
true

julia> has_semistables(K3, [1,4], theta)
false
```
"""
@memoize Dict function has_semistables(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, slope_denominator)
        subdimensionsBiggerSlope = filter(e -> slope(e, theta, slope_denominator) > slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d, algorithm="schofield"), subdimensionsBiggerSlope)
    end
end

"""Checks if Q has a theta-stable representation of dimension vector d.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3]; theta = [3,-2];

julia> has_stables(Q, d, theta)
true

julia> Q = mKroneckerquiver(2); d = [2,2]; theta = [1,-1];

julia> has_stables(Q, d, theta)
false

julia> has_semistables(Q, d, theta)
true
```
"""
@memoize Dict function has_stables(Q::Quiver, d::Vector{Int}, theta::Vector{Int}; slope_denominator::Function = sum)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, slope_denominator)
        subdimensions_bigger_or_equal_slope = filter(e -> slope(e, theta, slope_denominator) >= slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensions_bigger_or_equal_slope)
    end
end

"""
Checks if d is a Schur root for Q.

By a lemma of Schofield (See Lemma 4.2 of https://arxiv.org/pdf/0802.2147.pdf),
this is equivalent to the existence of a stable representation of dimension vector d
for the canonical stability parameter.
    
Examples:

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3];

julia> is_Schur_root(Q, d)
true
```
"""
is_Schur_root(Q::Quiver, d::Vector{Int}) = has_stables(Q, d, canonical_stability(Q, d))

function is_real_root(Q, d)
    return Euler_form(Q, d, d) == 1
end

function is_imaginary_root(Q, d)
    return Euler_form(Q, d, d) <= 0
end

function is_isotropic_root(Q, d)
    return Euler_form(Q, d, d) == 0
end


"""Checks if e is a generic subdimension vector of d.

A dimension vector e is called a generic subdimension vector of d if a generic representation
of dimension vector d possesses a subrepresentation of dimension vector e.

By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf)
e is a generic subdimension vector of d if and only if
``<e',d-e> \\geq 0``
for all generic subdimension vectors e' of e.
"""
@memoize Dict function is_generic_subdimension_vector(Q::Quiver, e::Vector{Int}, d::Vector{Int}; algorithm::String = "schofield")
    if e == d
        return true
    elseif all(ei == 0 for ei in e)
        # return false
        return true
    else
    # considering subdimension vectors that violate the numerical condition
    # TODO this filtering is inefficent.
    # For fixed d-e, this is a LINEAR form, we KNOW which eprimes violate the condition.
    # We should just check those.
        Euler_matrix_temp = Euler_matrix(Q) * (d - e) #to speed up computation of <eprime,d-e>
        subdimensions = filter(eprime -> eprime' * Euler_matrix_temp < 0, all_subdimension_vectors(e))
        # none of the subdimension vectors violating the condition should be generic
        return !any(eprime -> is_generic_subdimension_vector(Q, eprime, e), subdimensions)
    end
end

"""
Returns the list of all generic subdimension vectors of d.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3];

julia> QuiverTools.all_generic_subdimension_vectors(Q, d)
7-element Vector{Vector{Int64}}:
 [0, 0]
 [0, 1]
 [0, 2]
 [1, 2]
 [0, 3]
 [1, 3]
 [2, 3]
```
"""
function all_generic_subdimension_vectors(Q::Quiver, d::Vector{Int}) 
    return filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end

"""
Returns a list of all the Harder Narasimhan types of representations of Q with dimension vector d, with respect to the slope function theta/slope_denominator.

Examples:
```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> all_HN_types(Q, d, theta)
8-element Vector{Vector{Vector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

julia> all_HN_types(Q, [3,0], [0,0])
1-element Vector{Vector{Vector{Int64}}}:
 [[3, 0]]

julia> Q = three_vertex_quiver(1,4,1); d = [4,1,4];

julia> theta = canonical_stability(Q, d);

julia> all_HN_types(Q, d, theta)
106-element Vector{Vector{Vector{Int64}}}:
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
 [[2, 1, 1], [2, 0, 3]]
 [[2, 1, 1], [1, 0, 1], [1, 0, 2]]
 [[2, 1, 1], [2, 0, 2], [0, 0, 1]]
 [[4, 1, 2], [0, 0, 2]]
 [[2, 0, 1], [2, 1, 3]]
 [[2, 0, 1], [0, 1, 0], [2, 0, 3]]
 [[2, 0, 1], [1, 0, 1], [1, 1, 2]]
 [[2, 0, 1], [1, 1, 1], [1, 0, 2]]
 [[2, 0, 1], [2, 0, 2], [0, 1, 1]]
 [[2, 0, 1], [2, 1, 2], [0, 0, 1]]
 [[2, 0, 1], [2, 1, 1], [0, 0, 2]]
 [[4, 0, 2], [0, 1, 1], [0, 0, 1]]
 [[4, 0, 2], [0, 1, 0], [0, 0, 2]]
 [[3, 1, 1], [1, 0, 3]]
 [[3, 1, 1], [1, 0, 2], [0, 0, 1]]
 [[3, 1, 1], [1, 0, 1], [0, 0, 2]]
 ⋮
 [[2, 0, 0], [1, 0, 1], [1, 1, 3]]
 [[2, 0, 0], [1, 0, 1], [1, 0, 2], [0, 1, 1]]
 [[2, 0, 0], [1, 0, 1], [1, 1, 2], [0, 0, 1]]
 [[2, 0, 0], [1, 1, 1], [1, 0, 3]]
 [[2, 0, 0], [1, 1, 1], [1, 0, 2], [0, 0, 1]]
 [[2, 0, 0], [2, 0, 2], [0, 1, 1], [0, 0, 1]]
 [[2, 0, 0], [2, 1, 2], [0, 0, 2]]
 [[2, 0, 0], [2, 1, 1], [0, 0, 3]]
 [[2, 0, 0], [2, 0, 1], [0, 1, 1], [0, 0, 2]]
 [[2, 0, 0], [2, 0, 1], [0, 1, 0], [0, 0, 3]]
 [[2, 0, 0], [1, 1, 0], [1, 0, 4]]
 [[2, 0, 0], [1, 1, 0], [1, 0, 3], [0, 0, 1]]
 [[2, 0, 0], [1, 1, 0], [1, 0, 2], [0, 0, 2]]
 [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 0, 3]]
 [[3, 0, 0], [1, 1, 4]]
 [[3, 0, 0], [1, 1, 3], [0, 0, 1]]
 [[3, 0, 0], [1, 0, 2], [0, 1, 1], [0, 0, 1]]
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
@memoize Dict function all_HN_types(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function=sum; ordered=true)

    if all(di == 0 for di in d)
        return [[d]]
    else
        # We consider just proper subdimension vectors which admit a semistable representation and for which μ(e) > μ(d)
        # Note that we also eliminate d by the following
        subdimensions = filter(e -> has_semistables(Q,  e, theta, slope_denominator), all_destabilizing_subdimension_vectors(d, theta, slope_denominator)) 
        
        # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        if ordered
            subdimensions = sort(subdimensions, by = e -> slope(e, theta, slope_denominator))
        end

        # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN type of f = d-e and μ(e) > μ(f^1) holds.

        alltypes = Vector{Vector{Int}}[[e, efstar...] for e in subdimensions for efstar in filter(fstar -> slope(e, theta, slope_denominator) > slope(fstar[1], theta, slope_denominator), all_HN_types(Q, d-e, theta, slope_denominator, ordered=ordered))]

        # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
        if has_semistables(Q, d, theta, slope_denominator)
            pushfirst!(alltypes, [d])
        end
        return alltypes
    end
end

"""
Returns the codimension of the given HN stratum.

Examples:
```julia-repl
julia> Q = mKronecker_quiver(3); d = [2,3]; theta = [3,-2];

julia> HN = all_HN_types(Q, d, theta)
8-element Vector{Vector{Vector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]

 julia> for hntype in HN
 println(QuiverTools.codimension_HN_stratum(Q, hntype))
 end
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
function codimension_HN_stratum(Q::Quiver, stratum::Vector{Vector{Int}})
    if length(stratum) == 1
        return 0
    else
        return -sum(Euler_form(Q, stratum[i], stratum[j]) for i in 1:length(stratum)-1 for j in i+1:length(stratum))
    end
end

"""
Checks wether the dimension vector d is amply stable with respect to the slope function theta/denominator.
    
This means that the codimension of the unstable locus in the parameter space is at least 2.
"""
function is_amply_stable(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    # We say that representations of a given dimension vector d are amply stable (for any notion of stability) if the codimension of the semistable locus is at least 2.
    # We verify this by computing the codimension of each HN stratum.
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, slope_denominator))
    return all(stratum -> codimension_HN_stratum(Q, stratum) >= 2, HN)
end


######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################


""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type."""
function Teleman_bound_onstratum(Q::Quiver, hntype::Vector{Vector{Int}}, theta::Vector{Int}, slope_denominator::Function = sum)::Int
    if length(hntype) == 1
        throw(ArgumentError("Weight not defined for HN type of length 1."))
    end
    slopes = map(h -> slope(h, theta, slope_denominator), hntype)
    slopes = lcm(denominator.(slopes)) .* slopes
    return sum((slopes[t] - slopes[s]) * Euler_form(Q, hntype[s], hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype))
end

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given Q, d, theta and denominator."""
function all_Teleman_bounds(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, slope_denominator))
    return Dict([hntype, Teleman_bound_onstratum(Q, hntype, theta, slope_denominator)] for hntype in HN)
end

"""Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type."""
function weights_universal_bundle_onstratum(theta::Vector{Int}, a::Vector{Int}, hntype, slope_denominator::Function = sum)::Vector{Int}
    slopes = map(h -> slope(h, theta, slope_denominator), hntype)
    slopes *= lcm(denominator.(slopes))

    constant_term = sum(slopes[i]* (a' * hntype[i]) for i in eachindex(hntype))

    return -constant_term .+ slopes
end
"""Computes the weights of the universal bundle ``U_i(a)`` for the linearization ``a``
on all the non-dense Harder-Narasimhan strata for each 1-PS corresponding to each HN type."""
function all_weights_universal_bundle(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, a::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, slope_denominator))
    return Dict([hntype, weights_universal_bundle_onstratum(theta, a, hntype, slope_denominator)] for hntype in HN)
end


"""Computes the weight of the irreducible component of ``\\omega_R|_Z``
on a Harder-Narasimhan stratum for the 1-PS corresponding to each HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weight of
the pullback of O(H) on the given stratum."""
function weight_irreducible_component_canonical_on_stratum(Q::Quiver,d::Vector{Int},hntype::Vector{Vector{Int}},theta::Vector{Int},slope_denominator::Function = sum)::Int
    kweights = map(di -> slope(di,theta, slope_denominator), hntype)
    kweights = kweights * lcm(denominator.(kweights))

    dd = sum( kweights[m] .* hntype[m] for m in 1:length(hntype))
    # The Fano paper shows that under appropriate conditions,
    # the canonical bundle is given by linearizing with minus
    # the canonical stability parameter.
    can = canonical_stability(Q, d)
    can /= gcd(can)
    return can' * dd
end

"""Computes the weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of O(H) on each stratum."""
function all_weights_irreducible_component_canonical(Q::Quiver,d::Vector{Int},theta::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q,d,theta))
    return Dict([hntype, weight_irreducible_component_canonical_on_stratum(Q, d, hntype, theta, slope_denominator)] for hntype in HN)
end

"""Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type."""
function weights_endomorphism_universal_bundle_on_stratum(hntype::Vector{Vector{Int}},theta::Vector{Int},slope_denominator::Function = sum)::Vector{Int}
    # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is slope of first term in the HN type - slope of the last term in the HN type
    kweights = map(di -> slope(di,theta, slope_denominator), hntype)
    kweights = kweights * lcm(denominator.(kweights))
    # return kweights[1] - kweights[end] # this is the largest one
    return [kweights[i] - kweights[j] for i in 1:length(hntype) for j in 1:length(hntype)]
end

"""Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type."""
function all_weights_endomorphisms_universal_bundle(Q::Quiver,d::Vector{Int},theta::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], all_HN_types(Q, d, theta, slope_denominator))
    return Dict([hntype, weights_endomorphism_universal_bundle_on_stratum(hntype, theta, slope_denominator)] for hntype in HN)
end


#####################################################
# Canonical decomposition
#####################################################

"""
Computes the dimension of the ext group between generic representations
of dimension vectors a and b.

According to Thm. 5.4 in Schofield's 'General representations of quivers', we have
``ext(a,b) = max\\{- \\langle c, b\\rangle \\| c gen. subdimension vector of a\\}``.
"""
function generic_ext(Q::Quiver, a::Vector{Int}, b::Vector{Int})
    return maximum(-Euler_form(Q, c, b) for c in all_generic_subdimension_vectors(Q, a))
end

"""
Computes the dimension of the hom group between generic representations
of dimension vectors a and b.
"""
function generic_hom(Q::Quiver, a::Vector{Int}, b::Vector{Int})
    return Euler_form(Q, a, b) + generic_ext(Q, a, b)
end


# TODO add references
# TODO implement Derksen-Weyman?
# TODO add examples
"""
Computes the canonical decomposition of the dimension vector d for the given quiver Q.
If ``\\beta_1, \\dots, \\beta_{\\ell}`` is a sequence
of Schur roots such that, for all ``i \\neq j``, one has
``ext(\\beta_i, \\beta_j) = ext(\\beta_j, \\beta_i) = 0``,
then the general representation of dimension ``\\sum_i \\beta_i`` is a direct sum
of irreducible representations of dimension vectors ``\\beta_i``.
Such a decomposition is called therefore the canonical decomposition.
"""
function canonical_decomposition(Q::Quiver, d::Vector{Int})
    # if is_Schur_root(Q, d)
    #     return [d]
    # end
    generic_subdimensions = filter(e -> e != d, all_generic_subdimension_vectors(Q, d))
    for e in generic_subdimensions
        if d - e in generic_subdimensions && generic_ext(Q, e, d-e) == 0 && generic_ext(Q, d-e, e) == 0
            return vcat(canonical_decomposition(Q, e), canonical_decomposition(Q, d-e))
        end
    end
    return [d] # if nothing above worked then d is a Schur root.
end

# """
# Checks wether the stability parameter theta is on a wall with respect to the wall-and-chamber decomposition for the dimension vector d.
# The wall and chamber decomposition is described in Section 2.2, MR4352662
# """
# function is_on_a_fake_wall(d::Vector{Int}, theta::Vector{Int}) 
#     return any(e -> e'*theta == 0, all_proper_subdimension_vectors(d))
# end

"""
Given an HN type the quiver Q, returns the upper triangular matrix whose ij entry is ext(d^1,d^j) where (d^1,...,d^l) is the HN type.
The sum of all the entries is the codimension of the HN stratum; the sum of all the rectangles starting on the "up-diagonal" (where the 1s go in a Jordan form) and going all the way is at least 1.
Teleman is satisfied for this stratum iif one of these rectangles sums to 2 or more.
"""
function extension_matrix(Q::Quiver, hntype::Vector{Vector{Int}})
    n = length(hntype)

    if n <= 1
        throw(ArgumentError("HN type must have length at least 2, this makes no sense for the dense stratum"))
    else
        M = zeros(Int, n, n)
        for i in 1:n-1, j in i+1:n
            M[i,j] = -Euler_form(Q, hntype[i], hntype[j])
        end
        return M
    end
end



function is_luna_type(Q::Quiver, tau::Vector{Tuple{Vector{Int},Int}}, theta::Vector{Int})
    n = nvertices(Q)
    zeroVector = Vector{Int}(zeros(Int, n))
    d = sum(sum(tupl[1] .* tupl[2]) for tupl in tau) 
    if d == zeroVector
        return tau == [(zeroVector, 1)]
    else
        dstar = [tupl[1] for tupl in tau]
        equalSlope = all(e -> slope(e, theta, denominator) == slope(d, theta, denominator), dstar)
        semistable = all(e -> has_stable_representation(Q, e, theta, algorithm="schofield"), dstar)
        return equalSlope && semistable
    end
end

function all_luna_types(Q::Quiver, d::Vector{Int}, theta::Vector{Int})
    throw(ArgumentError("not implemented"))
end

function semistable_equals_stable(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, algorithm::String = "schofield")
    throw(ArgumentError("not implemented"))
end

"""
Checks if the dimension vector d is in the fundamental domain of the quiver Q.

The fundamental domain is the cone of dimension vectors in ``\\mathbb{Z}^{Q_0}``
such that the symmetric Tits form is negative on all the simple roots, i.e.,
for all vertices i, 

``(i, d) := \\langle d, s_i\\rangle + \\langle s_i, d\\rangle  \\leq 0``,

where ``s_i`` is the dimension vector with all entries set to 0, and the i-th
set to 1.
"""
function in_fundamental_domain(Q::Quiver, d::Vector{Int}; interior::Bool=false)
    # https://arxiv.org/abs/2209.14791 uses a strict inequality, while https://arxiv.org/abs/2310.15927 uses a non-strict.
    # here we set it to non-strict by default.
    
    simples = [unit_vector(nvertices(Q), i) for i in 1:nvertices(Q)]
    if interior
        return all(simple -> Euler_form(Q, d, simple) + Euler_form(Q, simple, d) < 0, simples)
    end
    return all(simple -> Euler_form(Q, d, simple) + Euler_form(Q, simple, d) <= 0, simples)
end

######################################################################################
# Below lie methods to compute Hodge diamonds translated from the Hodge diamond cutter.
# In turn, these are based on M. Reineke's paper
# "The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli", 
# https://doi.org/10.1007/s00222-002-0273-4
######################################################################################


###################################################
# auxiliary functions for HOdge_polynomial() below

"""
Solve Ax=b for A upper triangular via back substitution
"""
function solve(A,b)
    n = length(b)
    x = Vector{Any}(zeros(n))
    
    x[n] = b[n] / A[n,n]
    
    for i in n-1:-1:1
        x[i] = (b[i] - sum(A[i,j]*x[j] for j in i+1:n)) / A[i,i]
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
function CardinalRd(Q::Quiver, d::Vector{Int}, q)
    return q^sum(d[i] * d[j] * Q.adjacency[i,j] for i in 1:nvertices(Q), j in 1:nvertices(Q))
end

"""
Cardinality of product of general linear groups \$\\mathrm{GL}_{d}(\\mathbb{F}_q)\$.
"""
@memoize Dict function CardinalGd(d::Vector{Int}, q)
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

function Td(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, q)
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
    
# auxiliary functions for Hodge_polynomial() above
###################################################
    
"""
Returns the Hodge polynomial of the moduli space of theta-semistable representations of Q with dimension vector d.
The algorithm is based on [MR1974891](https://doi.org/10.1007/s00222-002-0273-4),
and the current implementation is translated from the [Hodge diamond cutter]
(https://zenodo.org/doi/10.5281/zenodo.3893509).
"""
function Hodge_polynomial(Q::Quiver, d::Vector{Int}, theta::Vector{Int})
    
    # safety checks
    if theta' * d == 0 && !is_coprime(d)
        throw(ArgumentError("d is not coprime"))
    elseif theta' * d != 0 && gcd(theta' * d, sum(d)) != 1
        throw(ArgumentError("d is not coprime in the sense of Definition 6.3."))
    elseif !is_acyclic(Q)
        throw(ArgumentError("Q is not acyclic."))
    end

    R, q = polynomial_ring(Nemo.QQ, 'q')
    F = fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?
    
    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])
    
    result = numerator(solve(T, one_at_the_end)[1] * (1-v))

    S,(x, y) = polynomial_ring(Nemo.QQ, ["x", "y"])
    return result(x*y)
    # return [coeff(result,i) for i in 0:degree(result)] # this is actually all we need for the Hodge diamond because the matrix is diagonal for quiver moduli
end

"""
Returns the Hodge diamond of the moduli space of theta-semistable representations of Q with dimension vector d.
"""
function Hodge_diamond(Q::Quiver, d::Vector{Int}, theta::Vector{Int})
    g = Hodge_polynomial(Q, d, theta)
    return map(ind -> coeff(g, [ind[1] - 1, ind[2] - 1]).num, Iterators.product(1:degree(g, 1) + 1, 1:degree(g, 2) + 1))
end

"""
Picard rank of the moduli space of theta-semistable representations of Q with dimension vector d.
"""
function Picard_rank(Q::Quiver, d::Vector{Int}, theta::Vector{Int})
    # TODO If over the complex numbers this should be h^{1,1}, since the moduli space is rational.
    # TODO This should follow from the long exact sequence in cohomology given by the exponential short exact sequence.

    return coeff(Hodge_polynomial(Q, d, theta), 2, 2).num
end

function _Hodge_polynomial_fast(Q::Quiver, d::Vector{Int}, theta::Vector{Int})
    # unsafe, curate input!
    # this is about 2% faster than the above, and occupies about 2% less memory.
    
    R, q = polynomial_ring(QQ, 'q')
    F = fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?
    
    T = Td(Q, d, theta, v)

    one_at_the_end = unit_vector(size(T)[1], size(T)[1])
    
    result = numerator(solve(T, one_at_the_end)[1] * (1-v))
    return [coeff(result, i) for i in 0:degree(result)] # this is actually all we need for the Hodge diamond because the matrix is diagonal for quiver moduli
end



#################
# Constructors
#################

function mKronecker_quiver(m::Int)
    return Quiver([0 m; 0 0], string(m)*"-Kronecker quiver")
end

function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "An acyclic 3-vertex quiver")
end

function loop_quiver(m::Int)
    return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m)*"-loop quiver")
end

function subspace_quiver(m::Int)
    A = zeros(Int, m+1, m+1)
    for i in 1:m
        A[i, m+1] = 1
    end
    return Quiver(A, string(m)*"-subspace quiver")
end

function DynkinQuiver(Tn::String)
    #parse the string Tn
    T = Tn[1:end-1]
    n = parse(Int, Tn[end])
    return DynkinQuiver(T,n)
end

function DynkinQuiver(T::String, n::Int)

    if T == "A"
        if !(n >= 1)
            throw(ArgumentError("$n is out of bounds"))
        end
        if n == 1
#            return Quiver([[1]], "Dynkin quiver of type A1")
            return LoopQuiver(1)
        else
            M = zeros(Int, n, n)
            for i in 1:n-1
                M[i, i+1] = 1
            end
            return Quiver(M, "Dynkin quiver of type A$n")
        end
    elseif T == "D"
        if !(n >= 3)
            throw(ArgumentError("$n is out of bounds."))
        end
        M = zeros(Int, n, n)
        for i in 1:n-2
            M[i, i+1] = 1
        end
        M[n-2, n] = 1

        return Quiver(M, "Dynkin quiver of type D$n")
    elseif T == "E"
        if !(n in [6, 7, 8])
            throw(ArgumentError("$n is out of bounds."))
        end
        if n == 6
            return Quiver([ 0 1 0 0 0 0 0;
                            0 0 1 0 0 0 0;
                            0 0 0 1 1 0 0;
                            0 0 0 0 0 0 0;
                            0 0 0 0 0 1 0;
                            0 0 0 0 0 0 1;
                            0 0 0 0 0 0 0], "Dynkin quiver of type E6")
        elseif n == 7
            return Quiver([ 0 1 0 0 0 0 0 0;
                            0 0 1 0 0 0 0 0;
                            0 0 0 1 1 0 0 0;
                            0 0 0 0 0 0 0 0;
                            0 0 0 0 0 1 0 0;
                            0 0 0 0 0 0 1 0;
                            0 0 0 0 0 0 0 1;
                            0 0 0 0 0 0 0 0], "Dynkin quiver of type E7")
        elseif n == 8
            return Quiver([ 0 1 0 0 0 0 0 0 0;
                            0 0 1 0 0 0 0 0 0;
                            0 0 0 1 1 0 0 0 0;
                            0 0 0 0 0 0 0 0 0;
                            0 0 0 0 0 1 0 0 0;
                            0 0 0 0 0 0 1 0 0;
                            0 0 0 0 0 0 0 1 0;
                            0 0 0 0 0 0 0 0 1;
                            0 0 0 0 0 0 0 0 0], "Dynkin quiver of type E8")
        end
    else
        throw(ArgumentError("not implemented"))
    end
end

#TODO: constructors
ExtendedDynkinQuiver(T::String) = throw(ArgumentError("not implemented"))

CyclicQuiver(n::Int) = throw(ArgumentError("not implemented"))

BipartiteQuiver(m::Int, n::Int) = throw(ArgumentError("not implemented"))

""""
Returns a Quiver with the same vertices and an arrow ``j \\to i`` for every arrow  ``i \\to j`` in the original quiver.
"""
opposite_quiver(Q::Quiver) = Quiver(Matrix{Int}(transpose(Q.adjacency)), "Opposite of "*Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) = Quiver(Q.adjacency + Matrix{Int}(transpose(Q.adjacency)), "Double of "*Q.name)


#################
# Technical tools
#################

@memoize Dict function zero_vector(n::Int)
    return Vector{Int}(zeros(Int, n))
end

function thin_dimension_vector(Q::Quiver)
    return ones(Int, nvertices(Q))
end

@memoize Dict function all_subdimension_vectors(d::Vector{Int})
    if any(di < 0 for di in d)
        throw(ArgumentError("dimension vector must be nonnegative"))
    end
    return collect.(Iterators.product(map(di -> 0:di, d)...))
end
@memoize Dict function all_nonzero_subdimension_vectors(d::Vector{Int})::Vector{Vector{Int}}
    return filter(e->!all(ei == 0 for ei in e), all_subdimension_vectors(d))
end

@memoize Dict function all_proper_subdimension_vectors(d::Vector{Int})::Vector{Vector{Int}}
    return filter(e -> any(ei != 0 for ei in e) && e != d, all_subdimension_vectors(d))
end

@memoize Dict function unit_vector(n::Int, i::Int)
    v = zeros(Int, n)
    v[i] = 1
    return v
end

function unit_vector(Q::Quiver, i::Int)
    return unit_vector(nvertices(Q), i)
end


end