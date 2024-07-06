### Here go the constructors for standard quivers.

#################
# Constructors
#################

export mKronecker_quiver, loop_quiver, subspace_quiver, three_vertex_quiver


"""
    mKronecker_quiver(m::Int)

Constructs a Kronecker quiver with `m` vertices.

# Arguments
- `m::Int`: (Default = 2) The number of vertices in the Kronecker quiver.

# Returns
A Kronecker quiver with `m` vertices.

# Example

```jldoctest
julia> mKronecker_quiver(3)
"""
function mKronecker_quiver(m::Int=2)
    return Quiver([0 m; 0 0], string(m) * "-Kronecker quiver")
end

"""
    three_vertex_quiver(m12::Int, m13::Int, m23::Int)

Constructs a three-vertex quiver with the given edge weights.

# Arguments
- `m12::Int`: The weight of the edge from vertex 1 to vertex 2.
- `m13::Int`: The weight of the edge from vertex 1 to vertex 3.
- `m23::Int`: The weight of the edge from vertex 2 to vertex 3.

# Returns
A three-vertex quiver object.

"""
function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "Acyclic 3-vertex quiver")
end

"""
    loop_quiver(m::Int)

Constructs a loop quiver with `m` vertices.

# Arguments
- `m::Int`: The number of vertices in the loop quiver.

# Returns
A loop quiver with `m` vertices.
"""
function loop_quiver(m::Int)
    return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m) * "-loop quiver")
end

"""
    subspace_quiver(m::Int)

Constructs a subspace quiver with `m` vertices.

# Arguments
- `m::Int`: The number of subspace-vertices.

# Returns
A subspace quiver with `m` subspaces.

# Examples
"""
function subspace_quiver(m::Int)
    A = zeros(Int, m + 1, m + 1)
    for i = 1:m
        A[i, m+1] = 1
    end
    return Quiver(A, string(m) * "-subspace quiver")
end

function Dynkin_quiver(Tn::String)
    #parse the string Tn
    T = Tn[1:end-1]
    n = parse(Int, Tn[end])
    return Dynkin_quiver(T, n)
end

function Dynkin_quiver(T::String, n::Int)

    if T == "A"
        if !(n >= 1)
            throw(ArgumentError("$n is out of bounds"))
        end
        if n == 1
            #            return Quiver([[1]], "Dynkin quiver of type A1")
            return loop_quiver(1)
        else
            M = zeros(Int, n, n)
            for i = 1:n-1
                M[i, i+1] = 1
            end
            return Quiver(M, "Dynkin quiver of type A$n")
        end
    elseif T == "D"
        if !(n >= 3)
            throw(ArgumentError("$n is out of bounds."))
        end
        M = zeros(Int, n, n)
        for i = 1:n-2
            M[i, i+1] = 1
        end
        M[n-2, n] = 1

        return Quiver(M, "Dynkin quiver of type D$n")
    elseif T == "E"
        if !(n in [6, 7, 8])
            throw(ArgumentError("$n is out of bounds."))
        end
        if n == 6
            return Quiver(
                [
                    0 1 0 0 0 0 0
                    0 0 1 0 0 0 0
                    0 0 0 1 1 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 1 0
                    0 0 0 0 0 0 1
                    0 0 0 0 0 0 0
                ],
                "Dynkin quiver of type E6",
            )
        elseif n == 7
            return Quiver(
                [
                    0 1 0 0 0 0 0 0
                    0 0 1 0 0 0 0 0
                    0 0 0 1 1 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 1 0 0
                    0 0 0 0 0 0 1 0
                    0 0 0 0 0 0 0 1
                    0 0 0 0 0 0 0 0
                ],
                "Dynkin quiver of type E7",
            )
        elseif n == 8
            return Quiver(
                [
                    0 1 0 0 0 0 0 0 0
                    0 0 1 0 0 0 0 0 0
                    0 0 0 1 1 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 1 0 0 0
                    0 0 0 0 0 0 1 0 0
                    0 0 0 0 0 0 0 1 0
                    0 0 0 0 0 0 0 0 1
                    0 0 0 0 0 0 0 0 0
                ],
                "Dynkin quiver of type E8",
            )
        end
    else
        throw(ArgumentError("not implemented"))
    end
end

#TODO: constructors
ExtendedDynkinQuiver(T::String) = throw(ArgumentError("not implemented"))

function CyclicQuiver(n::Int)
    if n < 1
        throw(ArgumentError("n must be greater than 0"))
    end
    A = zeros(Int, n, n)
    for i = 1:n-1
        A[i, i+1] = 1
    end
    A[n, 1] = 1
    return Quiver(A, "Cyclic quiver on $n vertices")
end

function BipartiteQuiver(m::Int, n::Int)
    if m < 1 || n < 1
        throw(ArgumentError("m and n must be greater than 0"))
    end
    A = zeros(Int, m+n, m+n)
    for i = 1:m
        for j = m+1:m+n
            A[i, j] = 1
        end
    end
    return Quiver(A, "Bipartite quiver on $m and $n vertices")
end

""""
Returns a Quiver with the same vertices and an arrow
``j \\to i`` for every arrow  ``i \\to j`` in the original quiver.

# Arguments
- `Q::Quiver`: The quiver to be reversed.

# Returns
A quiver with the same vertices and reversed arrows.
"""
opposite_quiver(Q::Quiver) =
    Quiver(Matrix{Int}(transpose(Q.adjacency)), "Opposite of " * Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of
the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) =
    Quiver(Q.adjacency + Matrix{Int}(transpose(Q.adjacency)), "Double of " * Q.name)
