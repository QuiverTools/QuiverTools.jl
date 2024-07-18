########################################################################################
# Definitions of types and primitive constructors for quivers and moduli spaces
########################################################################################

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
        pairs = map(p -> [parse(Int, p[1]), parse(Int, p[end]), length(p) - 1], pairs)

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





abstract type QuiverModuli end

# TODO consider this:
# https://stackoverflow.com/questions/71738970/in-julia-declare-abstractvectorabstractvector
# this is also necessary to be able to type function outputs correctly.
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
