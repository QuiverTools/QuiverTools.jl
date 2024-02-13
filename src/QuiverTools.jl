module QuiverTools
export Quiver, GeneralizedKroneckerQuiver, LoopQuiver, SubspaceQuiver, ThreeVertexQuiver, all_slope_decreasing_sequences, has_semistable_representation, AllHarderNarasimhanTypes, all_weight_bounds, does_teleman_inequality_hold, is_luna_type, all_luna_types, semistable_equals_stable, slope, all_subdimension_vectors, all_forbidden_subdimension_vectors, is_coprime_for_stability_parameter, is_indivisible, IsAcyclic, IsConnected, InDegree, OutDegree, IsSource, IsSink, euler_matrix, euler_form, opposite_quiver, double_quiver, canonical_decomposition, canonical_stability_parameter, is_harder_narasimhan_type, codimension_of_harder_narasimhan_stratum, is_amply_stable, GeneralizedKroneckerQuiver, LoopQuiver, SubspaceQuiver, thin_dimension_vectors, all_generic_subdimension_vectors, generic_ext_vanishing, is_generic_subdimension_vector, number_of_arrows, number_of_vertices, underlying_graph, ZeroVector, DynkinQuiver,CaseStudy, extension_matrix, HodgeDiamond

using LinearAlgebra, Nemo, Memoize


"""
A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij}),
``
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows i → j.

Attributes:

- `adjacency` is the adjacency matrix of the quiver
- `name` is the name of the quiver, defaults to `""`.
"""

mutable struct Quiver
    adjacency::Matrix{Int64}
    name::String

    
    function Quiver(adjacency::Matrix{Int64}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[1])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else 
            new(adjacency, name)
        end
    end 
    function Quiver(adjacency::Matrix{Int64})
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
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
    #TODO this is not necessary really. We should just use Q.adjacency
    return Q.adjacency
end

"""
Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.
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
Checks wether the quiver is acyclic, i.e. has no oriented cycles.
"""
IsAcyclic(Q::Quiver) = all(entry == 0 for entry in Q.adjacency^number_of_vertices(Q))

"""
Checks wether the underlying graph of the quiver is connected.

Examples:
```julia-repl
    julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0])
    julia> IsConnected(Q)
    true

    julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2])
    false

    # The 4-Kronecker quiver:
    julia> Q = GeneralizedKroneckerQuiver(4)
    julia> IsConnected(Q)
    true

    # The 4-loop quiver:
    julia> Q = LoopQuiver(4)
    julia> IsConnected(Q)
    true

    # The 4-subspace quiver:
    julia> Q = SubspaceQuiver(4)
    julia> IsConnected(Q)
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
    julia> IsConnected(A10)
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
    julia> IsConnected(A10)
    false
```
"""
function IsConnected(Q::Quiver)
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
```julia-repl
julia> Q = GeneralizedKroneckerQuiver(4)
julia> InDegree(Q, 1)
0
julia> InDegree(Q, 2)
4
```
"""
InDegree(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[:,j]) : throw(DomainError(j, "vertex index out of bounds"))
# TODO unnecessary checks
"""
Returns the number of outgoing arrows from the vertex i.

Examples:
```julia-repl
julia> Q = GeneralizedKroneckerQuiver(4)
julia> OutDegree(Q, 1)
4
julia> OutDegree(Q, 2)
0
```
"""
OutDegree(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? sum(adjacency_matrix(Q)[i,:]) : throw(DomainError(i, "vertex index out of bounds"))
# TODO unnecessary checks
"""
Checks if the vertex i is a source, i.e., a vertex with no incoming arrows.

Examples:
```julia-repl
julia> Q = GeneralizedKroneckerQuiver(4)
julia> IsSource(Q, 1)
true
julia> IsSource(Q, 2)
false
```
"""
IsSource(Q::Quiver, i::Int64) = (1 <= i & i<= number_of_vertices(Q)) ? InDegree(Q, i) == 0 : throw(DomainError(i, "vertex index out of bounds"))

"""
Checks if the vertex j is a sink, i.e., a vertex with no outgoing arrows.

Examples:
```julia-repl
julia> Q = GeneralizedKroneckerQuiver(4)
julia> IsSink(Q, 1)
false
julia> IsSink(Q, 2)
true
```
"""
IsSink(Q::Quiver, j::Int64) = (1 <= j & j<= number_of_vertices(Q)) ? OutDegree(Q, j) == 0 : throw(DomainError(j, "vertex index out of bounds"))

"""
Returns the Euler matrix of the quiver.


The Euler matrix of a quiver Q is defined as 
```math
E = I - A,
```
where ``A`` is the adjacency matrix of Q and ``I`` is the identity matrix of the same size as ``A``.
"""
@memoize Dict euler_matrix(Q::Quiver) = Matrix{Int64}(I, number_of_vertices(Q), number_of_vertices(Q)) - adjacency_matrix(Q)

"""
Computes the Euler form of the quiver for vectors x and y.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where E is the Euler matrix of the quiver.
"""
@memoize Dict euler_form(Q::Quiver, x::Vector{Int64}, y::Vector{Int64}) = (length(x) == number_of_vertices(Q) & length(y) == number_of_vertices(Q)) ? x'*euler_matrix(Q)*y : throw(DomainError("dimension vectors must have length equal to number of vertices"))

euler_form_first(Q::Quiver, y::Vector{Int64}) = x -> x'*euler_matrix(Q)*y

euler_form_second(Q::Quiver, x::Vector{Int64}) = y -> x'*euler_matrix(Q)*y


"""
Returns a Quiver with the same vertices and an arrow ``j \\to i`` for every arrow  ``i \\to j`` in the original quiver.
"""
opposite_quiver(Q::Quiver) = Quiver(Matrix{Int64}(transpose(adjacency_matrix(Q))), "Opposite of "*Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) = Quiver(adjacency_matrix(Q) + Matrix{Int64}(transpose(adjacency_matrix(Q))), "Double of "*Q.name)


## Everything that comes before this has been properly translated from the Sage code and should work.

thin_dimension_vectors(Q::Quiver) = ones(Int64,number_of_vertices(Q))

"""
The canonical stability parameter for the couple ``(Q,d)`` is given by ``<d,- > - < - ,d>``
"""
canonical_stability_parameter(Q::Quiver, d::Vector{Int64})::Vector{Int64} = -(-transpose(euler_matrix(Q)) + euler_matrix(Q))*d


"""
Returns the list of all sequences ``(d^1,...,d^l)`` which sum to d such that ``\\mu(d^1) > ... > \\mu(d^l).``

Examples:
```julia-repl
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
```
"""
function all_slope_decreasing_sequences(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum)

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

# the following code serves the computation of the directed graph of generic subdimension vectors. It is only called with the "gianni" algorithm for now.


function _indices_for_generic_subdimension_graph(d::Vector{Int64})::Dict{Vector{Int64},Int64}
    all_subdim_vectors = all_nonzero_subdimension_vectors(d,sorted=true)
    return Dict{Vector{Int64},Int64}(all_subdim_vectors[i] => i for i in eachindex(all_subdim_vectors))
end

"""
Returns the adjacency matrix of the directed graph of generic roots for dimension vectors up to d for the quiver Q.
"""
function _generic_subdimension_graph(Q::Quiver, d::Vector{Int64})::Matrix{Bool}
    # sort by sum to have adjacency matrix upper triangular.
    all_subdim_vectors = all_nonzero_subdimension_vectors(d,sorted=true)

    adjacency = zeros(Bool,(length(all_subdim_vectors), length(all_subdim_vectors)))
    # euler matrix here to not call it thousands of times in the loop
    # euler_matrix = copy(QuiverTools.euler_matrix(Q)) # maybe useful for threading
    euler_matrix = QuiverTools.euler_matrix(Q)

    for i in eachindex(all_subdim_vectors)
        #the list of all generic subdimension vectors of the one corresponding to i.
        generic_subdim_indices = [l for l in 1:i if adjacency[l,i] || l == i]

        for j in i:length(all_subdim_vectors)
            # this eliminates all cases where i is not componentwise smaller than j
            if all(all_subdim_vectors[i][k] <= all_subdim_vectors[j][k] for k in eachindex(all_subdim_vectors[i]))
                # fast euler form
                euler_matrix_temp = euler_matrix * (all_subdim_vectors[j] - all_subdim_vectors[i])
                adjacency[i,j] = all(all_subdim_vectors[l]' * euler_matrix_temp >= 0 for l in generic_subdim_indices)
            end
        end
    end
    return adjacency
end

"""Checks if there is a theta-semistable representation of dimension vector d.

Examples:
```julia-repl
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
```
"""
@memoize Dict function has_semistable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum; algorithm="schofield")
    if all(di == 0 for di in d)
        return true
    else
        if algorithm == "schofield"
            # collect the list of all subdimension vectors e of bigger slope than d
            subdimensionsBiggerSlope = filter(e -> slope(e, theta, denominator) > slope(d, theta, denominator), all_proper_subdimension_vectors(d))
            # to have semistable representations, none of the vectors above must be generic subdimension vectors.
            return all(e -> !is_generic_subdimension_vector(Q, e, d, algorithm="schofield"), subdimensionsBiggerSlope)

        elseif algorithm == "gianni"
            generic_roots = _generic_subdimension_graph(Q,d)
            indices_for = _indices_for_generic_subdimension_graph(d)

            # list of subdimension vectors with bigger slope than d
            subdimensionsBiggerSlope = filter(e -> slope(e, theta, denominator) > slope(d, theta, denominator), all_proper_subdimension_vectors(d))

            # none of the subdimension vectors with bigger slope than d should be a generic subdimension.
            return all(subdim -> !generic_roots[indices_for[subdim], indices_for[d]],subdimensionsBiggerSlope)
        else
            throw(ArgumentError("algorithm not recognized"))
        end
    end
end

"""Checks if Q has a theta-stable representation of dimension vector d."""
function has_stable_representation(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum, algorithm::String = "schofield")
    if all(di == 0 for di in d)
        return true
    else
        if algorithm == "schofield"
            # collect the list of all subdimension vectors e of bigger slope than d
            subdimensions_bigger_or_equal_slope = filter(e -> slope(e, theta, denominator) >= slope(d, theta, denominator), all_proper_subdimension_vectors(d))
            # to have semistable representations, none of the vectors above must be generic subdimension vectors.
            return all(e -> !is_generic_subdimension_vector(Q, e, d, algorithm="schofield"), subdimensions_bigger_or_equal_slope)

        elseif algorithm == "gianni"
            generic_roots = _generic_subdimension_graph(Q,d)
            indices_for = _indices_for_generic_subdimension_graph(d)
            
            # list of subdimension vectors with bigger slope than d
            subdimensions_bigger_or_equal_slope = filter(e -> slope(e, theta, denominator) >= slope(d, theta, denominator), all_proper_subdimension_vectors(d))
            # none of the subdimension vectors with bigger slope than d should be a generic subdimension.
            return all(subdim -> !generic_roots[indices_for[subdim], indices_for[d]],subdimensions_bigger_or_equal_slope)
        else
            throw(ArgumentError("algorithm not recognized"))
        end
    end
end

"""
Checks if d is a Schur root for Q.
By a lemma of Schofield (See Lemma 4.2 of https://arxiv.org/pdf/0802.2147.pdf),
this is equivalent to the existence of a stable representation of dimension vector d.
"""
is_schur_root(Q::Quiver, d::Vector{Int64}) = has_stable_representation(Q, d, canonical_stability_parameter(Q, d))

"""Checks if e is a generic subdimension vector of d.

A dimension vector e is called a generic subdimension vector of d if a generic representation
of dimension vector d possesses a subrepresentation of dimension vector e.

By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf)
e is a generic subdimension vector of d if and only if
``<e',d-e> \\geq 0``
for all generic subdimension vectors e' of e.
"""
@memoize Dict function is_generic_subdimension_vector(Q::Quiver, e::Vector{Int64}, d::Vector{Int64}; algorithm::String = "schofield")
    if e == d
        return true
    elseif all(ei==0 for ei in e)
        return false
    else
        if algorithm == "schofield"
        # considering subdimension vectors that violate the numerical condition
        # TODO this filtering is inefficent. For fixed d-e, this is a LINEAR form, we KNOW which eprimes violate the condition. We should just check those.
            euler_matrix_temp = euler_matrix(Q) * (d-e) #to speed up computation of <eprime,d-e>
            subdimensions = filter(eprime -> eprime'*euler_matrix_temp < 0, all_nonzero_subdimension_vectors(e))
            # none of the subdimension vectors violating the condition should be generic
            return !any(eprime -> is_generic_subdimension_vector(Q,eprime,e,algorithm="schofield"), subdimensions)

        elseif algorithm == "gianni"
            indices = _indices_for_generic_subdimension_graph(d)
            return _generic_subdimension_graph(Q,d)[indices[e], indices[d]]
        else
            throw(ArgumentError("algorithm not recognized"))
        end            
    end
end

"""
Returns the list of all generic subdimension vectors of d.
"""
function all_generic_subdimension_vectors(Q::Quiver, d::Vector{Int64}) 
    @warn "rewrite this function to use the directed graph of generic subdimension vectors"
    return filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end

function indices_for_generic_root_tree(d::Vector{Int64})::Dict{Vector{Int64},Int64}
    all_subdim_vectors = all_nonzero_subdimension_vectors(d,sorted=true)
    return Dict{Vector{Int64},Int64}(all_subdim_vectors[i] => i for i in eachindex(all_subdim_vectors))
end

"""
Returns the list of all dimension vectors of d that admit semistable representations.
"""
function all_semistable_dimension_vectors(Q::Quiver,d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)
    all_subdimensions = all_nonzero_subdimension_vectors(d, sorted=true)
    # the following line is extremely slow, we compute all cases at once with the generic root tree instead
    # stable = map(e -> has_semistable_representation(Q,e,theta,denominator), all_subdimensions)
    generic_roots = _generic_subdimension_graph(Q,d)
    index_for = indices_for_generic_root_tree(d)
    slopes = map(e -> slope(e,theta,denominator),all_subdimensions)

    # retains all the dimension vectors which only have generic subdimensions with smaller slope
    return filter(e -> all(slopes[j] <= slopes[index_for[e]] for j in 1:index_for[e] if generic_roots[index_for[e],j]), all_subdimensions)
end


generic_ext_vanishing(Q::Quiver, a::Vector{Int64}, b::Vector{Int64}) = is_generic_subdimension_vector(Q, a, a+b)

generic_hom_vanishing(Q::Quiver, a::Vector{Int64}, b::Vector{Int64}) = throw(ArgumentError("not implemented"))
"""
Two dimension vectors `a` and `b` are said to be left orthogonal if `hom(a,b) = ext(a,b) = 0`.
"""
function is_left_orthogonal(Q::Quiver, a::Vector{Int64}, b::Vector{Int64})
    if generic_ext_vanishing(Q, a, b)
        return EulerForm(Q, a, b) == 0
    else
        return false
    end
end

function is_real_root(Q::Quiver, a::Vector{Int64})
    return EulerForm(Q, a, a) == 1
end

function rearrange_dw_decomposition!(Q,decomposition,i,j)
    # apply Lemma 11.9.10 of Derksen--Weyman to rearrange the roots from i to j as k, k+1
    S = []
    for k in i:j-1
        # if k has a ``path to j'' satisfying the hypothesis of Lemma 11.9.10, we keep it
        # m is the j-k times j-k matrix with entry m[i,j] = 1 if !generic_hom_vanishing(Q,decomposition[j][2],decomposition[i][2]) and 0 otherwise
        m = zeros(Int64,j-k,j-k)
        for l in k:j-1 for s in l+1:j
            if !generic_hom_vanishing(Q,decomposition[s][2],decomposition[l][2])
                m[l-k+1,s-k+1] = 1
            end
        end
        paths = zeros(Int64,j-k,j-k) # paths[i,j] is the number of paths from k + i - 1 to k + j - 1 of length at most j-k
        for l in 1:j-k
            paths = paths + m^l
        end
        if paths[1,j-k] > 0
            push!(S,k)
        end
    end
    rearrangement = [l for l in i+1:j-1 if l ∉ S]
    final = [S...,i,j,rearrangement...]
    permute!(decomposition,final)
    return i + length(S) + 1 #this is the index of the element i now.
end
# the function below is written in Julia, and below we translate it to Sagemath


function canonical_decomposition(Q::Quiver, d::Vector{Int64}, algorithm::String = "derksen-weyman")
    if algorithm == "derksen-weyman"
        if any(Q.adjacency[i,j] != 0 for (i,j) in eachindex(Q.adjacency) if i <= j)
            throw(ArgumentError("Algorithm is meant for strictly lower triangular adjacency matrix"))
        end
        decomposition = [[d[i],UnitVector(length(d),i)] for i in 1:number_of_vertices(Q)]
        while true
            
            # check P1 to P4
            # None of this is actually used in the algorithm
            # might be useful for debugging though

            # if any(root[1] < 0 for root in decomposition)
            #     throw(ErrorException("P1 violated"))
            # elseif any(!is_schur_root(Q,root[2]) for root in decomposition)
            #     throw(ErrorException("P2 violated"))
            # elseif !all(is_left_orthogonal(Q,decomposition[i][2],decomposition[j][2]) for i in eachindex(decomposition) for j in i+1:length(decomposition))
            #     throw(ErrorException("P3 violated"))
            # elseif all(root[1] != 1 for root in decomposition if EulerForm(Q,root[2],root[2]) < 0)
            #     throw(ErrorException("P4 violated"))
            # end
            
            decomposition = filter(root -> root[1] > 0,decomposition)
            violating_pairs = [(i,j) for i in 1:length(decomposition)-1 for j in i+1:length(decomposition) if EulerForm(Q,decomposition[j][2],decomposition[i][2]) < 0]            
            if isempty(violating_pairs)
                break
            end
            violating_pairs = sort(violating_pairs, by = (i,j) -> j - i) # only need smallest actually, not to sort them
            
            i,j = violating_pairs[1]

            i = rearrange_dw_decomposition!(Q,decomposition,i,j)
            p,xi = decomposition[i]
            q,eta = decomposition[i+1]
        
            xireal = is_real_root(Q,xi)
            etareal = is_real_root(Q,eta)
            zeta = p*xi + q*eta
            if xireal & etareal
                # both roots real
                # TODO figure out the vague instructions in Derksen--Weyman
                discriminant = EulerForm(Q,zeta,zeta)
                if discriminant > 0
                    # TODO what
                elseif discriminant == 0
                    zetaprime = zeta ÷ gcd(zeta)
                    # replace xi, eta by zetaprime
                    deleteat!(decomposition, i+1)
                    decomposition[i] = [1,zetaprime]
                else
                    # replace xi, eta by zeta 
                    deleteat!(decomposition, i+1)
                    decomposition[i] = [1,zeta]
                end
            elseif xireal & !etareal
                # xi real, eta imaginary
                # as per P4, here the coefficent of eta is 1
                if p + q*EulerForm(Q,eta,xi) >= 0
                    # replace (xi,eta) with (eta - <eta,xi>xi,xi)
                    deleteat!(decomposition, i+1)
                    decomposition[i] = [1,eta - EulerForm(Q,eta,xi)*xi]
                else
                    # replace (xi,eta) with zeta
                    deleteat!(decomposition, i+1)
                    decomposition[i] = [1,zeta]
                end
            elseif !xireal & etareal
                # xi imaginary, eta real
                if q + p*EulerForm(Q,eta,xi) >= 0
                    # replace (xi,eta) with (eta,xi - <eta,xi>eta)
                    decomposition[i] = [1,eta]
                    decomposition[i+1] = [1,xi - EulerForm(Q,eta,xi)*eta]
                else
                    # replace (xi,eta) with zeta
                    deleteat!(decomposition, i+1)
                    decomposition[i] = [1,zeta]
                end
            elseif !xireal & !etareal
                # both roots imaginary
                # replace (xi,eta) with zeta
                deleteat!(decomposition, i+1)
                decomposition[i] = [1,zeta]
            end


        end
        return decomposition
    
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
    if length(dstar) == 1
        return has_semistable_representation(Q, dstar[1], theta,denominator)
    else
        for i in 1:length(dstar)-1
            if slope(dstar[i],theta, denominator) <= slope(dstar[i+1],theta, denominator)
                return false
            end
        end
        return all(e -> has_semistable_representation(Q, e, theta, denominator), dstar)
    end    
end

"""
Returns a list of all the Harder Narasimhan types of representations of Q with dimension vector d, with respect to the slope function theta/denominator.
"""
@memoize Dict function AllHarderNarasimhanTypes(Q::Quiver,d::Vector{Int64},theta::Vector{Int64}, denominator=sum; ordered=true)

    if all(di == 0 for di in d)
        return [[d]]
    else

        # We consider just proper subdimension vectors which admit a semistable representation and for which μ(e) > μ(d)
        # Note that we also eliminate d by the following
        subdimensions = filter(e -> has_semistable_representation(Q, e, theta, denominator), all_forbidden_subdimension_vectors(d,theta,denominator)) 
        
        # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        if ordered
            subdimensions = sort(subdimensions, by = e -> slope(e,theta,denominator))
        end

        # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN type of f = d-e and μ(e) > μ(f^1) holds.

        allHNtypes =  [[e,efstar...] for e in subdimensions for efstar in filter(fstar -> slope(e,theta, denominator) > slope(fstar[1],theta, 
        denominator) ,AllHarderNarasimhanTypes(Q, d-e, theta, denominator,ordered=ordered)) ]

        # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
        if has_semistable_representation(Q, d, theta, denominator)
            pushfirst!(allHNtypes, [d])
        end
        return allHNtypes
    end
end

"""
Returns the codimension of the given HN stratum.
"""
function codimension_of_harder_narasimhan_stratum(Q::Quiver, stratum::Vector{Vector{Int64}})
    if length(stratum) == 1
        return 0
    else
        return -sum(euler_form(Q, stratum[i],stratum[j]) for i in 1:length(stratum)-1 for j in i+1:length(stratum))
    end
end

"""
Checks wether the stability parameter theta is on a wall with respect to the wall-and-chamber decomposition for the dimension vector d.
The wall and chamber decomposition is described in Section 2.2, MR4352662
"""
function is_on_a_fake_wall(d::Vector{Int64}, theta::Vector{Int64}) 
    return any(e -> e'*theta == 0,all_proper_subdimension_vectors(d))
end

"""
Checks wether the dimension vector d is amply stable with respect to the slope function theta/denominator.
This means that the codimension of the unstable locus in the parameter space is at least 2.
"""
function is_amply_stable(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)
    # We say that representations of a given dimension vector d are amply stable (for any notion of stability) if the codimension of the semistable locus is at least 2.
    HN = filter(hntype -> hntype != [d], AllHarderNarasimhanTypes(Q, d, theta, denominator))
    return all(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) >= 2, HN)
end


""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda`` corresponding to the given HN type."""
function weight_bound(hntype::Vector{Vector{Int}}, theta::Vector{Int}, denominator::Function = sum)
    if length(hntype) == 1
        throw(ArgumentError("Weight not defined for HN type of length 1."))
    end
    return -sum((slope(hntype[s],theta,denominator) - slope(hntype[t],theta,denominator))*EulerForm(Q, hntype[s],hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype) )
end

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each HN type for the given Q, d, theta and denominator."""
function all_weight_bounds(Q::Quiver, d::Vector{Int}, theta::Vector{Int},denominator::Function = sum)

    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], AllHarderNarasimhanTypes(Q, d, theta, denominator))

    return map(hntype -> weight_bound(hntype,theta, denominator), HN)
end

function does_teleman_inequality_hold(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}; denominator::Function = sum)
    #This is only relevant on the unstable locus
    
    # We compute the weights of the 1-PS lambda on det(N_{S/R}|_Z) for each HN type
    weights = all_weight_bounds(Q, d, theta, denominator)
    
    HN = filter(hntype -> hntype != [d], AllHarderNarasimhanTypes(Q, d, theta, denominator))

    # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is slope of first term in the HN type - slope of the last term in the HN type
    tensorWeights = map(hntype -> slope(first(hntype),theta,denominator) - slope(last(hntype),theta,denominator), HN)

    @warn "Inequality has to be strict pheraps."
    return all(weights[i] >= tensorWeights[i] for i in eachindex(HN))

end

function _fano_paper_picard_rank(Q::Quiver,d::Vector{Int64})
    # MR4352662 gives the following easy computation. What to do for others?
    #TODO This should really be a method for a QuiverModuliSpace object. Translate the rest of the code?
    @info "Do this as in the hodge diamond cutter."
    
    theta=canonical_stability_parameter(Q,d)
    if is_coprime_for_stability_parameter(d,theta) && is_amply_stable(Q,d,theta)
        return number_of_vertices(Q) - 1
    else
       Throw(NotImplementedError("Not implemented for this stability parameter."))
    end
end


function does_Mukai_inequality_hold(Q::Quiver, d::Vector{Int64}; theta::Vector{Int64} = canonical_stability_parameter(Q,d), safe=false)
    # the first condition should really verify that theta is in the canonical chamber, but that is not implemented yet.
    # TODO: implement the canonical chamber check
    if safe
        run = theta == canonical_stability_parameter(Q,d) && is_coprime_for_stability_parameter(d, theta) && is_amply_stable(Q, d, theta)
    else
        run = true
    end
    if run
        PicardRank = number_of_vertices(Q) - 1
        Index = gcd(theta)
        return 1 - euler_form(Q,d,d) >= PicardRank*(Index - 1)
    else
        Throw(NotImplementedError("Not implemented for this stability parameter."))
    end
end


"""
Given an HN type the quiver Q, returns the upper triangular matrix whose ij entry is ext(d^1,d^j) where (d^1,...,d^l) is the HN type.
The sum of all the entries is the codimension of the HN stratum; the sum of all the rectangles starting on the "up-diagonal" (where the 1s go in a Jordan form) and going all the way is at least 1.
Teleman is satisfied for this stratum iif one of these rectangles sums to 2 or more.
"""
function extension_matrix(Q::Quiver, hntype::Vector{Vector{Int64}})
    n = length(hntype)

    if n <= 1
        throw(ArgumentError("HN type must have length at least 2, this makes no sense for the dense stratum"))
    else
        M = zeros(Int64, n, n)
        for i in 1:n-1, j in i+1:n
            M[i,j] = -euler_form(Q, hntype[i], hntype[j])
        end
        return M
    end
end



function is_luna_type(Q::Quiver, tau::Vector{Tuple{Vector{Int64},Int64}}, theta::Vector{Int64})
    n = number_of_vertices(Q)
    zeroVector = Vector{Int64}(zeros(Int64, n))
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

function all_luna_types(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    throw(ArgumentError("not implemented"))
end

function semistable_equals_stable(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64}, algorithm::String = "schofield")
    throw(ArgumentError("not implemented"))
end

slope(d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum) = (length(d) == length(theta) && denominator(d)>0) ? (theta'*d)//denominator(d) : throw(DomainError("different length entries or zero denominator"))

function in_fundamental_domain(Q::Quiver, d::Vector{Int64}; interior=false)
    # https://arxiv.org/abs/2209.14791 uses a strict inequality, while https://arxiv.org/abs/2310.15927 uses a non-strict?
    # here we set it to non-strict by default because.

    # there has to be a way to do this better
    simples = [zeros(Int64,number_of_vertices(Q)) for i in 1:number_of_vertices(Q)]

    for i in 1:number_of_vertices(Q)
        simples[i][i] = 1
    end
    if interior
        return all(simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) < 0, simples)
    elseif !interior
        return all(simple -> euler_form(Q, d, simple) + euler_form(Q, simple, d) <= 0, simples)
    end
    throw(ArgumentError("interior must be true or false"))
end

@memoize Dict function all_forbidden_subdimension_vectors(d::Vector{Int64}, theta::Vector{Int64}, denominator::Function = sum)
    return filter(e -> slope(e, theta,denominator) > slope(d, theta,denominator), all_proper_subdimension_vectors(d))
end

function is_coprime_for_stability_parameter(d::Vector{Int64}, theta::Vector{Int64})
    return all(e -> theta'*e != 0, all_proper_subdimension_vectors(d))
end

function is_indivisible(d::Vector{Int64})
    return gcd(d) == 1
end

# Below lie methods to compute Hodge diamonds translated from the Hodge diamond cutter.

"""
Solve Ax=b for A upper triangular via back substitution
"""
function solve(A,b)
    #TODO type
    n = length(b)
    x = Vector{Any}(zeros(n))

    x[n] = b[n] / A[n,n]

    for i in n-1:-1:1
        x[i] = (b[i] - sum(A[i,j]*x[j] for j in i+1:n)) / A[i,i]
    end
    return x
end

"""
Cardinality of general linear group \$\\mathrm{GL}_n(\\mathbb{F}_v)
"""
@memoize function CardinalGl(n::Int64,q)
    #TODO type
    if n == 0
        return 1
    else
        return prod(q^n - q^i for i in 0:n-1)
    end
end

"""
Cardinality of representation space \$\\mathrm{R}(Q,d)\$, over \$\\mathbb{F}_q\$.
"""
function CardinalRd(Q::Quiver, d::Vector{Int64},q)
    #TODO type
    return q^sum(d[i] * d[j] * adjacency_matrix(Q)[i,j] for i in 1:number_of_vertices(Q), j in 1:number_of_vertices(Q))
end

"""
Cardinality of product of general linear groups \$\\mathrm{GL}_{d}(\\mathbb{F}_q)\$.
"""
@memoize function CardinalGd(d::Vector{Int64},q)
    #TODO type
    return prod(CardinalGl(di,q) for di in d)
end

"""Entry of the transfer matrix, as per Corollary 6.9"""
function TransferMatrixEntry(Q, e, f, q)
    #TODO type
    fe = f - e

    if all(fei >= 0 for fei in fe)
        return q^euler_form(Q,-fe,e) * CardinalRd(Q, fe,q) / CardinalGd(fe,q)
    else
        return 0
    end
end

function TdChangeThisName(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64},q)

    # indexing set for the transfer matrix
    I = filter(e -> slope(e, theta) > slope(d, theta), all_proper_subdimension_vectors(d))
    append!(I, [d])
    pushfirst!(I, ZeroVector(number_of_vertices(Q)))

    # TODO this should just be a placeholder for the matrix
    T = Matrix{Any}(zeros(length(I),length(I)))

    for (i, Ii) in enumerate(I)
        for j in i:length(I)  # upper triangular
            T[i, j] = TransferMatrixEntry(Q, Ii, I[j],q)
        end
    end

    return T
end
"""Returns the Hodge polynomial of the moduli space of theta-semistable representations of Q with dimension vector d. 
"""
function hodge_polynomial(Q::Quiver,d::Vector{Int64},theta::Vector{Int64})

    # safety checks
    if !is_coprime_for_stability_parameter(d,theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif !IsAcyclic(Q)
        throw(ArgumentError("Q is not acyclic"))
    end

    R,q = polynomial_ring(Nemo.QQ,'q')
    F = fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?

    T = TdChangeThisName(Q,d,theta,v)

    # there HAS to be a better way to do this
    one_at_the_end = [0 for i in 1:size(T)[1]]
    one_at_the_end[end] = 1

    result = numerator(solve(T, one_at_the_end)[1] * (1-v))
    S,(x, y) = PolynomialRing(Nemo.QQ, ["x", "y"])
    return result(x*y)
    # return [coeff(result,i) for i in 0:degree(result)] # this is actually all we need for the Hodge diamond because the matrix is diagonal for quiver moduli
end

function hodge_diamond_from_polynomial(g)
    #TODO can be bypassed by knowing that this matrix is diagonal. Do this before high performance computations!
    result = zeros(degree(g,1)+1,degree(g,2)+1)
    for i in 1:degree(g,1)+1
        for j in 1:degree(g,2)+1
            # .num is to extract the Int64 
            result[i,j] = coeff(g,[i-1,j-1]).num
        end
    end
    return result
end

"""
Returns the Hodge diamond of the moduli space of theta-semistable representations of Q with dimension vector d.
"""
function HodgeDiamond(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    return hodge_diamond_from_polynomial(hodge_polynomial(Q,d,theta))
end

"""
Picard rank of the moduli space of theta-semistable representations of Q with dimension vector d.
"""
function PicardRank(Q::Quiver, d::Vector{Int64}, theta::Vector{Int64})
    # TODO If over the complex numbers this should be h^{1,1}, since the moduli space is rational.
    # TODO This should follow from the long exact sequence in cohomology given by the exponential short exact sequence.

    # safety checks
    if !is_coprime_for_stability_parameter(d,theta)
        throw(ArgumentError("d and theta are not coprime"))
    elseif !IsAcyclic(Q)
        throw(ArgumentError("Q is not acyclic"))
    end

    # find Hodge polynomial as in th diamond cutter
    R,q = polynomial_ring(Nemo.QQ,'q')
    F = fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?
    
    T = TdChangeThisName(Q,d,theta,v)
    
    # there HAS to be a better way to do this
    one_at_the_end = [0 for i in 1:size(T)[1]]
    one_at_the_end[end] = 1

    # extract coefficent of Hodge polynomial of degree 1 in v (= h^{1,1})
    return coeff(numerator(solve(T, one_at_the_end)[1] * (1-v)),1).num
end

function _picard_rank_fast(Q::Quiver,d::Vector{Int64},theta::Vector{Int64})
    # internal, comes with no safety checks.
    R,q = polynomial_ring(Nemo.QQ,'q')
    F = fraction_field(R)
    v = F(q) # worsens performance by ~8%. Necessary?

    T = TdChangeThisName(Q,d,theta,v)

    # there HAS to be a better way to do this
    one_at_the_end = [0 for i in 1:size(T)[1]]
    one_at_the_end[end] = 1

    return coeff(numerator(solve(T, one_at_the_end)[1] * (1-v)),1).num
end


function GeneralizedKroneckerQuiver(m::Int64)
    return Quiver([0 m; 0 0], string(m)*"-Kronecker quiver")
end

function KroneckerQuiver()
    return GeneralizedKroneckerQuiver(2)
end

function ThreeVertexQuiver(m12::Int64, m13::Int64, m23::Int64)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "An acyclic 3-vertex quiver")
end

function LoopQuiver(m::Int64)
    # @warn "Assess behaviour of 1x1 matrices in Julia. Add tests!"
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

@memoize function ZeroVector(n::Int64)
    return Vector{Int64}(zeros(Int64, n))
end

# subdimension vectors @memoize Dict 
@memoize Dict function all_subdimension_vectors(d::Vector{Int64})
    return collect.(Iterators.product(map(di -> 0:di, d)...))
end
@memoize Dict function all_nonzero_subdimension_vectors(d::Vector{Int64};sorted::Bool=false)::Vector{Vector{Int64}}
    if sorted
        return sort(all_nonzero_subdimension_vectors(d), by=e -> sum(e))
    end
    return filter(e->!all(ei == 0 for ei in e),collect.(Iterators.product(map(di -> 0:di, d)...)))
end
@memoize Dict function all_proper_subdimension_vectors(d::Vector{Int64})::Vector{Vector{Int64}}
    return filter(e -> any(ei != 0 for ei in e) && e != d, all_subdimension_vectors(d))
end
@memoize Dict function UnitVector(n::Int64, i::Int64)
    v = zeros(Int64, n)
    v[i] = 1
    return v
end


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

    if has_semistable_representation(Q, d, theta, denominator)
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
        HN = AllHarderNarasimhanTypes(Q, d, theta, denominator)
        @info "There are $(length(HN)) Harder--Narasimhan types: $HN."
        ample_stability = all(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) >= 2, filter(hntype -> hntype != [d], HN))
        @info "The lowest codimension is $(minimum(map(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum), filter(hntype -> hntype != [d], HN)))), and the stratum that gives it is $(argmin(stratum -> codimension_of_harder_narasimhan_stratum(Q,stratum), filter(hntype -> hntype != [d], HN)))."

        teleman_inequality = all(stratum -> sum((slope(stratum[t],theta,denominator) - slope(stratum[s],theta,denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) > slope(first(stratum), theta, denominator) - slope(last(stratum), theta, denominator), filter(hntype -> hntype != [d], HN))
        
        if !ample_stability
           failed_element_ample_stability = findfirst(stratum -> codimension_of_harder_narasimhan_stratum(Q, stratum) < 2, filter(hntype -> hntype != [d], HN))
              @info "Ample stability failed for element: $failed_element_ample_stability."
        else 
            @info "Ample stability holds."
        end

        if !teleman_inequality && rigidity
            failed_element_rigidity = findfirst(stratum -> sum((slope(stratum[t],theta,denominator) - slope(stratum[s],theta,denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) <= slope(first(stratum), theta, denominator) - slope(last(stratum), theta, denominator), filter(hntype -> hntype != [d], HN))
            @info "Rigidity failed for HN type $(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]): slope(first(stratum), theta, denominator) = $(slope(first(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]), theta, denominator)), slope(last(stratum), theta, denominator) = $(slope(last(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]), theta, denominator)), sum((slope(stratum[t],theta,denominator) - slope(stratum[s],theta,denominator))*euler_form(Q, stratum[s],stratum[t]) for s in 1:length(stratum)-1 for t in s+1:length(stratum) ) = $(sum((slope(filter(hntype -> hntype != [d], HN)[failed_element_rigidity][t],theta,denominator) - slope(filter(hntype -> hntype != [d], HN)[failed_element_rigidity][s],theta,denominator))*euler_form(Q, filter(hntype -> hntype != [d], HN)[failed_element_rigidity][s],filter(hntype -> hntype != [d], HN)[failed_element_rigidity][t]) for s in 1:length(filter(hntype -> hntype != [d], HN)[failed_element_rigidity])-1 for t in s+1:length(filter(hntype -> hntype != [d], HN)[failed_element_rigidity]) ))."
        elseif rigidity
            @info "Rigidity holds."
        end

    end

    if returning
        if !strong_Ample_Stability
            return strong_Ample_Stability, ample_stability, teleman_inequality, HN
        else
            return strong_Ample_Stability, ample_stability, teleman_inequality
        end
    end
end


end