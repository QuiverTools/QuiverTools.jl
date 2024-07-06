### Here go the constructors for standard quivers.

#################
# Constructors
#################

function mKronecker_quiver(m::Int)
    return Quiver([0 m; 0 0], string(m) * "-Kronecker quiver")
end

function three_vertex_quiver(m12::Int, m13::Int, m23::Int)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "Acyclic 3-vertex quiver")
end

function loop_quiver(m::Int)
    return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m) * "-loop quiver")
end

function subspace_quiver(m::Int)
    A = zeros(Int, m + 1, m + 1)
    for i in 1:m
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
            return Quiver([0 1 0 0 0 0 0;
                    0 0 1 0 0 0 0;
                    0 0 0 1 1 0 0;
                    0 0 0 0 0 0 0;
                    0 0 0 0 0 1 0;
                    0 0 0 0 0 0 1;
                    0 0 0 0 0 0 0], "Dynkin quiver of type E6")
        elseif n == 7
            return Quiver([0 1 0 0 0 0 0 0;
                    0 0 1 0 0 0 0 0;
                    0 0 0 1 1 0 0 0;
                    0 0 0 0 0 0 0 0;
                    0 0 0 0 0 1 0 0;
                    0 0 0 0 0 0 1 0;
                    0 0 0 0 0 0 0 1;
                    0 0 0 0 0 0 0 0], "Dynkin quiver of type E7")
        elseif n == 8
            return Quiver([0 1 0 0 0 0 0 0 0;
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
opposite_quiver(Q::Quiver) = Quiver(Matrix{Int}(transpose(Q.adjacency)), "Opposite of " * Q.name)

"""
The adjacency matrix of the double of a quiver is the sum of the adjacency matrix of the original quiver and its transpose.
"""
double_quiver(Q::Quiver) = Quiver(Q.adjacency + Matrix{Int}(transpose(Q.adjacency)), "Double of " * Q.name)
