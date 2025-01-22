#####################################################
# Methods that deal with quivers and their properties
#####################################################

function deglex_key(Q::Quiver, e::AbstractVector{Int})
  b = maximum(e) + 1
  n = n_vertices(Q)

  return Int(sum(e[i] * b^(n - i) for i in 1:length(e)) + sum(e) * b^n)
end

# TODO why is this a symmetric matrix? this should be `underlying_undirectd_graph`?
# a "graph", as opposed to a "directed graph (i.e., a quiver)", is undirected - and thus the adjacency matrix is symmetric
"""
    underlying_graph(Q::Quiver)

Return the (necessarily symmetric) adjacency matrix of the underlying graph of `Q`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> underlying_graph(Q) == [0 4; 4 0]
true
```
"""
function underlying_graph(Q::Quiver)
  return Matrix{Int}(Q.adjacency + transpose(Q.adjacency) - diagonal(Q.adjacency))
end

"""
    n_vertices(Q::Quiver)

Return the number of vertices of `Q`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> n_vertices(Q) == 2
true
```
"""
n_vertices(Q::Quiver) = size(Q.adjacency)[1]

"""
    n_arrows(Q::Quiver)

Return the number of arrows of `Q`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> n_arrows(Q) == 4
true
```
"""
n_arrows(Q::Quiver) = sum(Q.adjacency)

"""
    is_acyclic(Q::Quiver)

Check whether `Q` is acyclic, i.e., has no oriented cycles.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_acyclic(Q)
true
```
"""
is_acyclic(Q::Quiver) = all(entry == 0 for entry in Q.adjacency^n_vertices(Q))

"""
    is_connected(Q::Quiver)

Check whether `Q` is connected.

# Examples

```jldoctest
julia> is_connected(Quiver([0 1 0; 0 0 1; 1 0 0]))
true

julia> is_connected(Quiver([0 1 0; 1 0 0; 0 0 2]))
false

julia> is_connected(kronecker_quiver(4))
true

julia> is_connected(loop_quiver(4))
true

julia> is_connected(subspace_quiver(4))
true
```
"""
function is_connected(Q::Quiver)
  paths = underlying_graph(Q)
  for i in 2:(n_vertices(Q) - 1)
    paths += paths * underlying_graph(Q)
  end
  for i in 1:n_vertices(Q), j in 1:n_vertices(Q)
    if i != j && paths[i, j] == 0 && paths[j, i] == 0
      return false
    end
  end
  return true
end

"""
    indegree(Q::Quiver, j::Int)

Return the number of incoming arrows to the vertex `j`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> indegree(Q, 1)
0

julia> indegree(Q, 2)
4
```
"""
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j])

"""
    outdegree(Q::Quiver, i::Int)

Return the number of outgoing arrows from the vertex `i`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> outdegree(Q, 1)
4

julia> outdegree(Q, 2)
0
```
"""
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :])

"""
    is_source(Q::Quiver, i::Int)

Check whether the vertex `i` is a source, i.e., a vertex with no incoming arrows.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_source(Q, 1)
true

julia> is_source(Q, 2)
false
```
"""
is_source(Q::Quiver, i::Int) = indegree(Q, i) == 0

"""
    is_sink(Q::Quiver, j::Int)

Checks if the vertex `j` is a sink, i.e., a vertex with no outgoing arrows.

# Examples

```jldoctest
julia> Q = kronecker_quiver(4);

julia> is_sink(Q, 1)
false

julia> is_sink(Q, 2)
true
```
"""
is_sink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

"""
    arrows(Q::Quiver)

Return a list of all arrows of `Q`.

# Examples

```jldoctest
julia> arrows(kronecker_quiver(3))
3-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 2)
 (1, 2)

julia> arrows(loop_quiver(3))
3-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 1)
 (1, 1)
```
"""
function arrows(Q::Quiver)
  n = n_vertices(Q)
  return reduce(
    vcat,
    [(i, j) for k in 1:Q.adjacency[i, j]]
    for i in 1:n for j in 1:n if Q.adjacency[i, j] > 0
  )
end

"""
    first_hochschild_cohomology(Q::Quiver)

Compute the first Hochschild cohomology group.

The Hochschild cohomology groups of an acyclic quiver `Q` are described by Happel in
[Proposition 1.6, MR1035222](https://mathscinet.ams.org/mathscinet/relay-station?mr=1035222)
to be

```math
\\mathrm{H}^{0}(Q) = k,
\\quad \\mathrm{H}^{1}(Q) = 1 - n + \\sum_{i,j} #\\{\\text{paths from } i \\text{ to } j\\},
\\quad \\mathrm{H}^{n}(Q) = 0 \\text{ for } n > 1.
```

# Examples

On the 3-Kronecker quiver:

```jldoctest
julia> first_hochschild_cohomology(kronecker_quiver(3))
8

julia> first_hochschild_cohomology(subspace_quiver(7))
0

julia> first_hochschild_cohomology(three_vertex_quiver(2, 1, 3))
18
```
"""
function first_hochschild_cohomology(Q::Quiver)
  !is_acyclic(Q) && throw(ArgumentError("The quiver must be acyclic."))
  return 1 - n_vertices(Q) + sum(
    (Q.adjacency^n)[i, j]
    for n in 1:(n_vertices(Q) - 1), (i, j) in arrows(Q)
  )
end
