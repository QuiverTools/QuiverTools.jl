#################
# Constructors
#################

export mKronecker_quiver,
    loop_quiver,
    subspace_quiver,
    three_vertex_quiver,
    cyclic_quiver,
    bipartite_quiver,
    opposite_quiver,
    double_quiver,
    Dynkin_quiver


"""
    mKronecker_quiver(m::Int)

Constructs a Kronecker quiver with `m` vertices.

INPUT:
- `m`: (Default = 2) The number of arrows in the Kronecker quiver.

OUTPUT:
A Kronecker quiver with `m` vertices.

EXAMPLES:
```jldoctest
julia> mKronecker_quiver(3)
3-Kronecker quiver, with adjacency matrix [0 3; 0 0]
```
"""
function mKronecker_quiver(m::Int = 2)
    return Quiver([0 m; 0 0], string(m) * "-Kronecker quiver")
end

"""
    three_vertex_quiver(m12::Int, m13::Int, m23::Int)

Constructs a three-vertex quiver with the given edge weights.

INPUT:
- `m12`: The number of arrows from vertex 1 to vertex 2.
- `m13`: The number of arrows from vertex 1 to vertex 3.
- `m23`: The number of arrows from vertex 2 to vertex 3.

OUTPUT:
A three-vertex quiver with the specified arrows.

EXAMPLES:
```jldoctest
julia> three_vertex_quiver(1, 2, 3)
Acyclic 3-vertex quiver, with adjacency matrix [0 1 2; 0 0 3; 0 0 0]
```
"""
function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "Acyclic 3-vertex quiver")
end

"""
    loop_quiver(m::Int)

Constructs a loop quiver with `m` vertices.

INPUT:
- `m`: The number of vertices in the loop quiver.

OUTPUT:
A loop quiver with `m` vertices.

EXAMPLES:
```jldoctest
julia> loop_quiver(4)
4-loop quiver, with adjacency matrix [4;;]
```
"""
function loop_quiver(m::Int)
    return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m) * "-loop quiver")
end

"""
    subspace_quiver(m::Int)

Constructs a subspace quiver with `m` vertices.

INPUT:
- `m`: The number of subspace-vertices.

OUTPUT:
A subspace quiver with `m` subspaces.


EXAMPLES:
```jldoctest
julia> subspace_quiver(3)
3-subspace quiver, with adjacency matrix [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0]
```
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

"""
    Dynkin_quiver(T, n)

Constructs the Dynkin quiver, with arbitrary orientation of the arrows.

EXAMPLES:
```jldoctest
julia> Dynkin_quiver("D4")
Dynkin quiver of type D4, with adjacency matrix [0 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 0 0]
```
"""
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
"""
    cyclic_quiver(n)

Returns a cyclic quiver on n vertices.

EXAMPLES:
```jldoctest
julia> cyclic_quiver(4)
cyclic quiver on 4 vertices, with adjacency matrix [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0]
```
"""
function cyclic_quiver(n::Int)
    if n < 1
        throw(ArgumentError("n must be greater than 0"))
    end
    A = zeros(Int, n, n)
    for i = 1:n-1
        A[i, i+1] = 1
    end
    A[n, 1] = 1
    return Quiver(A, "cyclic quiver on $n vertices")
end

"""
    bipartite_quiver(m, n)

Constructs the bipartite quiver on `m` and `n` vertices.

EXAMPLES:
```jldoctest
julia> bipartite_quiver(2, 3).adjacency
5×5 Matrix{Int64}:
 0  0  1  1  1
 0  0  1  1  1
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
```
"""
function bipartite_quiver(m::Int, n::Int)
    if m < 1 || n < 1
        throw(ArgumentError("m and n must be greater than 0"))
    end
    A = zeros(Int, m + n, m + n)
    for i = 1:m
        for j = m+1:m+n
            A[i, j] = 1
        end
    end
    return Quiver(A, "bipartite quiver on $m and $n vertices")
end

""""
Returns a Quiver with the same vertices and an arrow
``j \\to i`` for every arrow  ``i \\to j`` in the original quiver.

INPUT:
- `Q::Quiver`: The quiver to be reversed.

OUTPUT:
A quiver with the same vertices and reversed arrows.

EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver()
2-Kronecker quiver, with adjacency matrix [0 2; 0 0]

julia> opposite_quiver(Q)
opposite of 2-Kronecker quiver, with adjacency matrix [0 0; 2 0]
```
"""
opposite_quiver(Q::Quiver) =
    Quiver(Matrix{Int}(transpose(Q.adjacency)), "opposite of " * Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of
the adjacency matrix of the original quiver and its transpose.


EXAMPLES:
```jldoctest
julia> Q = mKronecker_quiver();

julia> double_quiver(Q)
double of 2-Kronecker quiver, with adjacency matrix [0 2; 2 0]
```
"""
double_quiver(Q::Quiver) =
    Quiver(Q.adjacency + Matrix{Int}(transpose(Q.adjacency)), "double of " * Q.name)
