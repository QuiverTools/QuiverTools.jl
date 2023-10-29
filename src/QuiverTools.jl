# This code is a translation of the Sage code in quivers.py. It is not necessarily meant to be published, but rather to help with high-performance hypothesis testing and profiling.

module QuiverTools
export Quiver, GeneralizedKroneckerQuiver, LoopQuiver, SubspaceQuiver, ThreeVertexQuiver, all_slope_decreasing_sequences, has_semistable_representation, all_harder_narasimhan_types, all_weight_bounds, does_rigidity_inequality_hold, is_luna_type, all_luna_types, semistable_equals_stable, slope, all_subdimension_vectors, all_forbidden_subdimension_vectors, is_coprime_for_stability_parameter, is_indivisible, is_acyclic, is_connected, indegree, outdegree, is_source, is_sink, euler_matrix, euler_form, opposite_quiver, double_quiver, canonical_decomposition, canonical_stability_parameter, is_harder_narasimhan_type, codimension_of_harder_narasimhan_stratum, is_amply_stable, GeneralizedKroneckerQuiver, LoopQuiver, SubspaceQuiver, thin_dimension_vectors, all_generic_subdimension_vectors, generic_ext_vanishing, is_generic_subdimension_vector, number_of_arrows, number_of_vertices, underlying_graph, ZeroVector, DynkinQuiver

using LinearAlgebra


"""
A quiver is represented by its adjacency matrix (a_ij) in M_{n x n}(N) where Q_0 = {1,...,n} and a_{ij} is the number of arrows i --> j.

Variables:
`adjacency::Matrix{Int64}` is the adjacency matrix of the quiver
`name = None`
"""
mutable struct Quiver
    adjacency::Matrix{Int64}
    name::String

    
    function Quiver(adjacency::Matrix{Int64}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[1])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        elseif !(all(a >= 0 for a in adjacency))
            throw(DomainError(adjacency, "adjacency matrix must have non-negative entries"))
        else 
            new(adjacency, name)
        end
    end 
    function Quiver(adjacency::Matrix{Int64})
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        elseif !(all(a >= 0 for a in adjacency))
            throw(DomainError(adjacency, "adjacency matrix must have non-negative entries"))
        else
            new(adjacency, "")
        end
    end
end

"""
Returns the adjacency matrix of the quiver.

OUTPUT: A square matrix M whose entry M[i,j] is the number of arrows from the vertex i to the vertex j.
"""
function adjacency_matrix(Q::Quiver)
    return Q.adjacency
end

"""
Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.

OUTPUT: A square, symmetric matrix M whose entry M[i,j] = M[j,i] is the number of edges between the vertices i and j.
"""
function underlying_graph(Q::Quiver)
    return Matrix{Int64}(Q.adjacency + transpose(Q.adjacency) - diagm(diag(Q.adjacency)))
end
"""
Returns the number of vertices of the quiver.
"""
number_of_vertices(Q::Quiver) = size(Q.adjacency)[1]

"""
Returns the number of arrows of the quiver.
"""
number_of_arrows(Q::Quiver) = sum(Q.adjacency)

"""
Returns true if the quiver is acyclic, false otherwise.
"""
is_acyclic(Q::Quiver) = (Q.adjacency^number_of_vertices(Q) == zeros(Int64, number_of_vertices(Q), number_of_vertices(Q)))

"""
Returns true if the quiver is connected, false otherwise.

Examples:
    julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0])
    julia> is_connected(Q)
    true

    julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2])
    false

    The 4-Kronecker quiver:
    julia> Q = GeneralizedKroneckerQuiver(4)
    julia> is_connected(Q)
    true

    The 4-loop quiver:
    julia> Q = LoopQuiver(4)
    julia> is_connected(Q)
    true

    The 4-subspace quiver:
    julia> Q = SubspaceQuiver(4)
    julia> is_connected(Q)
    true

    The A10 quiver:

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

    The A10 quiver without one arrow:
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
    
"""
function is_connected(Q::Quiver)
    paths = underlying_graph(Q)
    for i in 2:number_of_vertices(Q) - 1
        paths += paths*underlying_graph(Q)
    end
    for i in 1:number_of_vertices(Q), j in 1:number_of_vertices(Q)
            if i != j && paths[i,j] == 0 && paths[j,i] == 0
                return false
            end
    end
    return true
end

# the docstrings on these functions are from the file quiver.py

"""
Returns the number of incoming arrows to the vertex j.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> indegree(Q, 1)
0
julia> indegree(Q, 2)
4
"""
indegree(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[:,j]) : throw(DomainError(j, "vertex index out of bounds"))

"""
Returns the number of outgoing arrows from the vertex i.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> outdegree(Q, 1)
4
julia> outdegree(Q, 2)
0
"""
outdegree(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[i,:]) : throw(DomainError(i, "vertex index out of bounds"))

"""
Returns true if the vertex i is a source, i.e. there are no incoming arrows into i, false otherwise.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> is_source(Q, 1)
true
julia> is_source(Q, 2)
false
"""
is_source(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? indegree(Q, i) == 0 : throw(DomainError(i, "vertex index out of bounds"))

"""
Returns true if the vertex j is a sink, i.e. there are no outgoing arrows from j, false otherwise.

Examples:

julia> Q = GeneralizedKroneckerQuiver(4)
julia> is_sink(Q, 1)
false
julia> is_sink(Q, 2)
true
"""
is_sink(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? outdegree(Q, j) == 0 : throw(DomainError(j, "vertex index out of bounds"))

"""
Returns the Euler matrix of the quiver.
"""
euler_matrix(Q::Quiver) = Matrix{Int64}(I, number_of_vertices(Q), number_of_vertices(Q)) - adjacency_matrix(Q)

"""
Returns the value of the Euler bilinear form of the quiver computed on the vectors x and y.
"""
euler_form(Q::Quiver, x::Vector{Int64}, y::Vector{Int64}) = (length(x) == number_of_vertices(Q) & length(y) == number_of_vertices(Q)) ? x'*euler_matrix(Q)*y : throw(DomainError("dimension vectors must have length equal to number of vertices"))

"""
The opposite quiver is given by the transpose of the adjacency matrix of the original quiver.

Returns a Quiver object with the same vertices and an arrow from j to i for every arrow from i to j in the original quiver.
"""
opposite_quiver(Q::Quiver) = Quiver(Matrix{Int64}(transpose(adjacency_matrix(Q))), "Opposite of "*Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) = Quiver(adjacency_matrix(Q) + Matrix{Int64}(transpose(adjacency_matrix(Q))), "Double of "*Q.name)


## Everything that comes before this has been properly translated from the Sage code and should work.

thin_dimension_vectors(Q::Quiver) = [1 for i in 1:number_of_vertices(Q)]

"""
The canonical stability parameter is given by <d,_> - <_,d>
"""
canonical_stability_parameter(Q::Quiver, d::Vector{Int64})::Vector{Int64} = -(-transpose(euler_matrix(Q)) + euler_matrix(Q))*d


"""
Returns the list of all sequences (d^1,...,d^l) which sum to d such that slope(d^1) > ... > slope(d^l)

Examples:

julia> Q = GeneralizedKroneckerQuiver(3)
julia> d = [2,3]
julia> theta = [3,-2]
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
    """
function all_slope_decreasing_sequences(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)

    # List all subdimension vectors e of bigger slope than d.
    subdimensions = filter(e -> (e != ZeroVector(number_of_vertices(Q))) && (slope(e,theta,denominator=denominator) > slope(d,theta,denominator=denominator)), all_subdimension_vectors(d))

    # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
    subdimensions = sort(subdimensions, by = e -> slope(e,theta,denominator=denominator))
    # The slope decreasing sequences which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

    # I will rewrite this as functional programming later
    allSlopeDecreasing = []
    for e in subdimensions
        for fstar in filter(fstar -> slope(e,theta,denominator=denominator) > slope(fstar[1],theta, denominator=denominator), all_slope_decreasing_sequences(Q, d-e, theta, denominator=denominator))
        push!(allSlopeDecreasing, [e, fstar...])
        end
    end
    # Add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
    return [[d], allSlopeDecreasing...]
end


"""Checks if there is a theta-semistable representation of dimension vector d.

Examples:


julia> A2 = GeneralizedKroneckerQuiver(1)
julia> theta = [1,-1]
julia> d = [1,1]
julia> has_semistable_representation(A2, d, theta)
true

julia> d = [2,2]
julia> has_semistable_representation(A2, d, theta)
true

julia> d = [1,2]
julia> has_semistable_representation(A2, d, theta)
false

julia> d = [0,0]
julia> has_semistable_representation(A2, d, theta)
true

The 3-Kronecker quiver:
julia> K3 = GeneralizedKroneckerQuiver(3)
julia> theta = [3,-2]
julia> d = [2,3]
julia> has_semistable_representation(K3, d, theta)
true

julia> d = [1,4]
julia> has_semistable_representation(K3, d, theta)
false
"""
function has_semistable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function =  sum, algorithm::String = "schofield") 
    if algorithm == "reineke"
        throw(ArgumentError("reineke algorithm not implemented"))

    elseif algorithm == "schofield"
        # collect the list of all subdimension vectors e of bigger slope than d
        subdimensionsBiggerSlope = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e, theta, denominator=denominator) > slope(d, theta, denominator=denominator), all_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return !any(e -> is_generic_subdimension_vector(Q, e, d), subdimensionsBiggerSlope)
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end

function has_stable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum, algorithm::String = "schofield")

    if algorithm == "al"
        throw(ArgumentError("al algorithm not implemented"))
    elseif algorithm == "schofield"
        if d == ZeroVector(number_of_vertices(Q))
            return false
        else
            subdimensionsSlopeNoLess = filter(e -> e != zeroVector && e != d && slope(e, theta, denominator=denominator) >= slope(d, theta, denominator=denominator), all_subdimension_vectors(d))
            return !any(e -> is_generic_subdimension_vector(Q, e, d), subdimensionsSlopeNoLess)
        end
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end


is_schur_root(Q::Quiver, d::Vector{Int64}) = has_stable_representation(Q, d, canonical_stability_parameter(Q, d))

"""Checks if e is a generic subdimension vector of d.
        # using notation from Section 5 of https://arxiv.org/pdf/0802.2147.pdf
    A dimension vector e is called a generic subdimension vector of d if a generic representation of dimension vector d possesses a subrepresentation of dimension vector e.
    By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf) e is a generic subdimension vector of d if and only if <e',d-e> is non-negative for all generic subdimension vectors e' of e."""
function is_generic_subdimension_vector(Q::Quiver, e::Vector{Int64}, d::Vector{Int64})
    if e == d
        return true
    else
    # considering subdimension vectors that violate the numerical condition
    subdimensions = filter(eprime -> euler_form(Q,eprime, d-e) < 0, all_subdimension_vectors(e))
    # none of the subdimension vectors violating the condition should be generic
    return !any(eprime -> is_generic_subdimension_vector(Q,eprime,e), subdimensions)
end
end
all_generic_subdimension_vectors(Q::Quiver, d::Vector{Int64}) = filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))


generic_ext_vanishing(Q::Quiver, a::Vector{Int64}, b::Vector{Int64}) = is_generic_subdimension_vector(Q, a, a+b)

function canonical_decomposition(Q::Quiver, d::Vector{Int64}, algorithm::String = "derksen-weyman")
    if algorithm == "derksen-weyman"
        throw(ArgumentError("derksen-weyman algorithm not implemented"))
    elseif algorithm == "schofield-1"
        throw(ArgumentError("schofield-1 algorithm not implemented"))
    elseif algorithm == "schofield-2"
        throw(ArgumentError("schofield-2 algorithm not implemented"))
    else
        throw(ArgumentError("algorithm not recognized"))
    end
end

"""
Checks wether dstar is a Harder--Narasimhan type of Q, with dimension vector d, with respect to the slope function theta/denominator
"""
function is_harder_narasimhan_type(Q::Quiver, dstar::Vector{Vector{Int64}}, theta::Vector{Int64}; denominator::Function = sum)
    @warn "This method is not relevant: it does not guarantee that the tuple of dimension vectors is a Harder--Narasimhan type for some representation."
    if length(dstar) == 1
        return has_semistable_representation(Q, dstar[1], theta)
    else
        for i in 1:length(dstar)-1
            if slope(dstar[i],theta, denominator=denominator) <= slope(dstar[i+1],theta, denominator=denominator)
                return false
            end
        end
        return all(e -> has_semistable_representation(Q, e, theta, denominator=denominator), dstar)
    end    
end

"""
Returns a list of all the Harder Narasimhan types of representations of Q with dimension vector d, with respect to the slope function theta/denominator.
"""
function all_harder_narasimhan_types(Q::Quiver,d::Vector{Int64},theta::Vector{Int64}; denominator::Function = sum)

    if d == ZeroVector(number_of_vertices(Q))
        return [[ZeroVector(number_of_vertices(Q))]]
    else

        # We consider just those subdimension vectors which are not zero, whose slope is bigger than the slope of d and which admit a semi-stable representation
        # Note that we also eliminate d by the following
        subdimensions = filter(e -> (e != ZeroVector(number_of_vertices(Q))) && (slope(e, theta, denominator=denominator) > slope(d,theta,denominator=denominator)) && has_semistable_representation(Q, e, theta, denominator=denominator, algorithm="schofield"), all_subdimension_vectors(d))
        # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        subdimensions = sort(subdimensions, by = e -> slope(e,theta,denominator=denominator))
        # The HN types which are not of the form (d) are given by (e,f^1,...,f^s) where e is a proper subdimension vector such that mu_theta(e) > mu_theta(d) and (f^1,...,f^s) is a HN type of f = d-e such that mu_theta(e) > mu_theta(f^1) holds.

        allHNtypes = []
        for e in subdimensions
            append!(allHNtypes, map(effestar -> [e,effestar...], filter(fstar -> slope(e,theta, denominator=denominator) > slope(fstar[1],theta, denominator=denominator) ,all_harder_narasimhan_types(Q, d-e, theta, denominator=denominator))))
#                push!(allHNtypes, [e, fstar...])

        end
        
        # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
        if has_semistable_representation(Q, d, theta, denominator=denominator, algorithm="schofield")
            pushfirst!(allHNtypes, [d])
        end
        return allHNtypes
    end
end

"""
Returns the codimension of the given HN stratum.
"""
function codimension_of_harder_narasimhan_stratum(Q::Quiver, stratum::Vector{Vector{Int64}})
    return -sum(euler_form(Q, stratum[i],stratum[j]) for i in 1:length(stratum)-1 for j in i+1:length(stratum))
end

"""
Checks wether the stability parameter theta is on a wall with respect to the wall-and-chamber decomposition for the dimension vector d.
The wall and chamber decomposition is described in Section 2.2, MR4352662
"""
is_on_a_wall(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}) = any(e -> sum(e .* theta) == 0,all_proper_subdimension_vectors(d))


"""
Checks wether the dimension vector d is amply stable with respect to the slope function theta/denominator. This means that the codimension of the unstable locus in the parameter space is at least 2.
"""
function is_amply_stable(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)
    # We say that representations of a given dimension vector d are amply stable (for any notion of stability) if the codimension of the semistable locus is at least 2.
    # We verify this by computing the codimension of each HN stratum.
    HN = filter(hntype -> hntype != [d], all_harder_narasimhan_types(Q, d, theta, denominator=denominator))
    return all(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) >= 2, HN)
end


function all_weight_bounds(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)

    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], all_harder_narasimhan_types(Q, d, theta, denominator=denominator))

    return map(hntype -> -sum((slope(hntype[s],theta,denominator=denominator) - slope(hntype[t],theta,denominator=denominator))*euler_form(Q, hntype[s],hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype) ), HN)

end
    
function does_rigidity_inequality_hold(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)
    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], all_harder_narasimhan_types(Q, d, theta, denominator=denominator))

    # The following code is here, commented out, for readability. It is equivalent to the one-liner below.

    # # We compute the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type
    # weights = map(hntype -> -sum((slope(hntype[s],theta,denominator=denominator) - slope(hntype[t],theta,denominator=denominator))*euler_form(Q, hntype[s],hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype) ), HN)
    
    # # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is slope of first term in the HN type - slope of the last term in the HN type
    # tensorWeights = map(hntype -> slope(first(hntype),theta,denominator=denominator) - slope(last(hntype),theta,denominator=denominator), HN)

    # return all(weights[i] >= tensorWeights[i] for i in 1:length(HN))


    return all( hntype -> sum((slope(hntype[t],theta,denominator=denominator) - slope(hntype[s],theta,denominator=denominator))*euler_form(Q, hntype[s],hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype) ) > slope(first(hntype), theta, denominator=denominator) - slope(last(hntype), theta, denominator=denominator), HN)

end



function is_luna_type(Q::Quiver, tau::Vector{Tuple{Vector{Int64},Int64}}, theta::Vector{Int64})
    n = number_of_vertices(Q)
    zeroVector = Vector{Int64}(zeros(Int64, n))
    d = sum(sum(tupl[1] .* tupl[2]) for tupl in tau) 
    if d == zeroVector
        return tau == [(zeroVector, 1)]
    else
        dstar = [tupl[1] for tupl in tau]
        equalSlope = all(e -> slope(e, theta, denominator=sum) == slope(d, theta, denominator=sum), dstar)
        semistable = all(e -> has_stable_representation(Q, e, theta, algorithm="schofield"), dstar)
        return equalSlope && semistable
    end
end

function all_luna_types(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    throw(ArgumentError("not implemented"))
end

function semistable_equals_stable(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield")
    throw(ArgumentError("not implemented"))
end

slope(d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum) = (length(d) == length(theta) && denominator(d)>0) ? (theta'*d)//denominator(d) : throw(DomainError("dimension vector and stability parameter must have same length"))

function in_fundamental_domain(Q::Quiver, d::Vector{Int64}; strict=false)
    # https://arxiv.org/abs/2209.14791 uses a strict inequality, while https://arxiv.org/abs/2310.15927 uses a non-strict?
    # here we set it to non-strict by default because.

    # there has to be a way to do this better
    simples = [ZeroVector(number_of_vertices(Q)) for i in ZeroVector(number_of_vertices(Q))]
    for i in 1:number_of_vertices(Q)
        simples[i][i] = 1
    end
    if strict
        return all(euler_form(Q, d, simple) + euler_form(Q, simple, d) < 0 for simple in simples)
    else
        return all(euler_form(Q, d, simple) + euler_form(Q, simple, d) <= 0 for simple in simples)
    end
end

function all_forbidden_subdimension_vectors(d::Vector{Int64}, theta::Vector{Int64})
    zeroVector = Vector{Int64}(zeros(Int64, length(d)))
    properSubdimensions = filter(e -> e != d && e != zeroVector, all_subdimension_vectors(d))
    return filter(e -> slope(e, theta) > slope(d, theta), properSubdimensions)
end

function is_coprime_for_stability_parameter(d::Vector{Int64}, theta::Vector{Int64})
    zeroVector = Vector{Int64}(zeros(Int64, length(d)))
    properSubdimensions = filter(e -> e != d && e != zeroVector, all_subdimension_vectors(d))
    return all(e -> sum(theta .*e) != 0, properSubdimensions)
end

function is_indivisible(d::Vector{Int64})
    return gcd(d) == 1
end

function GeneralizedKroneckerQuiver(m::Int64)
    return Quiver([0 m; 0 0], string(m)*"-Kronecker quiver")
end

KroneckerQuiver() = GeneralizedKroneckerQuiver(2)

ThreeVertexQuiver(m12::Int64, m13::Int64, m23::Int64) = Quiver([0 m12 m13; 0 0 m23; 0 0 0], "An acyclic 3-vertex quiver")

function LoopQuiver(m::Int64)
    @warn "Assess behaviour of 1x1 matrices in Julia. Add tests!"
    return Quiver(Matrix{Int64}(reshape([m],1,1)), string(m)*"-loop quiver")
end

function SubspaceQuiver(m::Int64)
    A = zeros(Int64, m+1, m+1)
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

function DynkinQuiver(T::String,n::Int64)

    if T == "A"
        if !(n >= 1)
            throw(ArgumentError("$n is out of bounds"))
        end
        if n == 1
#            return Quiver([[1]], "Dynkin quiver of type A1")
            return LoopQuiver(1)
        else
            M = zeros(Int64, n, n)
            for i in 1:n-1
                M[i, i+1] = 1
            end
            return Quiver(M, "Dynkin quiver of type A$n")
        end
    elseif T == "D"
        if !(n >= 3)
            throw(ArgumentError("$n is out of bounds."))
        end
        M = zeros(Int64, n, n)
        for i in 1:n-2
            M[i, i+1] = 1
        end
        M[n-2, n] = 1

        return Quiver(M, "Dynkin quiver of type D$n")
    elseif T == "E"
        if !(n in [6,7,8])
            throw(ArgumentError("$n is out of bounds."))
        end
        if n == 6
            return Quiver([0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 1 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1; 0 0 0 0 0 0 0], "Dynkin quiver of type E6")
        elseif n == 7
            return Quiver([0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 1 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 0 0], "Dynkin quiver of type E7")
        elseif n == 8
            return Quiver([0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0; 0 0 0 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 0 0 0], "Dynkin quiver of type E8")
        end
    else
        throw(ArgumentError("not implemented"))
    end
end

ExtendedDynkinQuiver(T::String) = throw(ArgumentError("not implemented"))
CyclicQuiver(n::Int64) = throw(ArgumentError("not implemented"))
BipartiteQuiver(m::Int64, n::Int64) = throw(ArgumentError("not implemented"))

ZeroVector(n::Int64) = Vector{Int64}(zeros(Int64, n))

# subdimension vectors

all_subdimension_vectors(d::Vector{Int64}) = collect(collect.(Iterators.product((0:di for di in d)...)))
all_nonzero_subdimension_vectors(d::Vector{Int64}) = filter(e -> e != ZeroVector(length(d)), all_subdimension_vectors(d))
all_proper_subdimension_vectors(d::Vector{Int64}) = filter(e -> e != ZeroVector(length(d)) && e != d, all_subdimension_vectors(d))



"""
This is a method for studying various properties of the moduli problem (Q,d,theta). It saves the results on a text file, or prints them without saving.
This method exists to avoid hardcoding multiple scripts that do the same thing, but with different parameters, or slighty different things each time.
"""
function CaseStudy(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum, coprimality=false,returning=false,rigidity=false, strong_AS=false)

    # Here there should be several boolean variables as input that determine what to return. For now, I will just return everything.
    @info "Running CaseStudy for $(adjacency_matrix(Q)) with dimension vector $d and stability parameter $theta."

    
    if is_coprime_for_stability_parameter(d, theta) && coprimality
         @info "The dimension vector $d and stability parameter $theta are coprime."
    elseif coprimality 
        @info "The dimension vector $d and stability parameter $theta are not coprime."
    end

    if has_semistable_representation(Q, d, theta, denominator=denominator)
        @info "The dimension vector $d admits a semistable representation with respect to the stability parameter $theta."
        @info "The moduli space has dimension $(1 - euler_form(Q,d,d))."
    else
        @info "The dimension vector $d does not admit a semistable representation with respect to the stability parameter $theta."
    end

    candidates_strong_AS = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e,theta) >= slope(d-e,theta),all_subdimension_vectors(d))
    strong_Ample_Stability = all(e -> euler_form(Q,e,d-e) <= -2, candidates_strong_AS)


    if strong_Ample_Stability
        @info "Strong ample stability holds. Ample stability and rigidity also hold."
    elseif !strong_Ample_Stability
        failed_element_strong_AS = findfirst(e -> euler_form(Q,e,d-e) > -2, candidates_strong_AS)
        if strong_AS
            @info "Strong ample stability failed for e = $(candidates_strong_AS[failed_element_strong_AS]): slope(e,theta) = $(slope(candidates_strong_AS[failed_element_strong_AS],theta)), slope(d-e,theta) = $(slope(d-candidates_strong_AS[failed_element_strong_AS],theta)), euler_form(Q,e,d-e) = $(euler_form(Q,candidates_strong_AS[failed_element_strong_AS],d-candidates_strong_AS[failed_element_strong_AS]))."
        end
        HN = all_harder_narasimhan_types(Q, d, theta, denominator=denominator)
        @info "There are $(length(HN)) Harder--Narasimhan types: $HN."
        ample_stability = all(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) >= 2, filter(hntype -> hntype != [d], HN))
        @info "The lowest codimension is $(minimum(map(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum), filter(hntype -> hntype != [d], HN)))), and the stratum that gives it is $(argmin(stratum -> codimension_of_harder_narasimhan_stratum(Q,stratum), filter(hntype -> hntype != [d], HN)))."

        rigidity_inequality = all(stratum -> sum((slope(stratum[t],theta,denominator=denominator) - slope(stratum[s],theta,denominator=denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) > slope(first(stratum), theta, denominator=denominator) - slope(last(stratum), theta, denominator=denominator), filter(hntype -> hntype != [d], HN))
        
        if !ample_stability
           failed_element_ample_stability = findfirst(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) < 2, filter(hntype -> hntype != [d], HN))
              @info "Ample stability failed for element: $failed_element_ample_stability."
        else 
            @info "Ample stability holds."
        end

        if !rigidity_inequality && rigidity
            failed_element_rigidity = findfirst(stratum -> sum((slope(stratum[t],theta,denominator=denominator) - slope(stratum[s],theta,denominator=denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) <= slope(first(stratum), theta, denominator=denominator) - slope(last(stratum), theta, denominator=denominator), filter(hntype -> hntype != [d], HN))
            @info "Rigidity failed for HN type $(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]): slope(first(stratum), theta, denominator=denominator) = $(slope(first(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]), theta, denominator=denominator)), slope(last(stratum), theta, denominator=denominator) = $(slope(last(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]), theta, denominator=denominator)), sum((slope(stratum[t],theta,denominator=denominator) - slope(stratum[s],theta,denominator=denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) = $(sum((slope(filter(hntype -> hntype != [d], HN)[failed_element_rigidity][t],theta,denominator=denominator) - slope(filter(hntype -> hntype != [d], HN)[failed_element_rigidity][s],theta,denominator=denominator))*euler_form(Q, filter(hntype -> hntype != [d], HN)[failed_element_rigidity][s],filter(hntype -> hntype != [d], HN)[failed_element_rigidity][t]) for s in 1:length(filter(hntype -> hntype != [d], HN)[failed_element_rigidity])-1 for t in s+1:length(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]) ))."
        elseif rigidity
            @info "Rigidity holds."
        end

    end

    if returning
        if !strong_Ample_Stability
            return strong_Ample_Stability, ample_stability, rigidity_inequality, HN
        else
            return strong_Ample_Stability, ample_stability, rigidity_inequality
        end
    end
end


end