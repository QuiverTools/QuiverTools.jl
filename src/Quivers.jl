#####################################################
# Methods that deal with quivers and their properties
#####################################################

function deglex_key(Q::Quiver, e::AbstractVector{Int})::Int
  b = maximum(e) + 1
  n = n_vertices(Q)

  return (sum(e[i] * b^(n - i) for i in 1:length(e)) + sum(e) * b^n)
end

"""
    underlying_graph(Q::Quiver)

Returns the (necessarily symmetric) adjacency matrix
of the underlying graph of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the adjacency matrix of the underlying graph of the quiver

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

Returns the number of vertices of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the number of vertices of the quiver

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

Returns the number of arrows of the quiver.

# Input

- `Q::Quiver` a quiver

# Output

- the number of arrows of the quiver

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

Checks whether the quiver is acyclic, i.e. has no oriented cycles.

# Input

- `Q::Quiver` a quiver

# Output

- `true` if the quiver is acyclic, `false` otherwise

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

Checks whether the quiver is connected.

# Input

- `Q::Quiver` a quiver

# Output

- `true` if the underlying graph of the quiver is connected, `false` otherwise

# Examples

```jldoctest
julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0]);

julia> is_connected(Q)
true

julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2]);

julia> is_connected(Q)
false

julia> # The 4-Kronecker quiver:

julia> Q = kronecker_quiver(4);

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

Returns the number of incoming arrows to the vertex `j`.

# Input

- `Q::Quiver` a quiver
- `j::Int`: an integer representing the vertex

# Output

- the number of incoming arrows to the vertex `j`

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

Returns the number of outgoing arrows from the vertex `i`.

# Input

- `Q::Quiver` a quiver
- `i::Int`: an integer representing the vertex

# Output

- the number of outgoing arrows from the vertex `i`

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

Checks if the vertex `i` is a source, i.e., a vertex with no incoming arrows.

# Input

- `Q::Quiver` a quiver
- `i::Int`: an integer representing the vertex

# Output

- `true` if the vertex `i` is a source, `false` otherwise

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

# Input

- `Q::Quiver` a quiver
- `j::Int`: an integer representing the vertex

# Output

- `true` if the vertex `j` is a sink, `false` otherwise

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

Returns a list of all arrows of the quiver `Q`.

# Input

- `Q::Quiver` a quiver

# Output

- a list of all arrows of the quiver `Q`.

# Examples

```jldoctest
julia> Q = kronecker_quiver(3);

julia> arrows(Q)
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 2]
 [1, 2]
```
"""
function arrows(Q::Quiver)
  n = n_vertices(Q)
  return reduce(
    vcat,
    [[i, j] for k in 1:Q.adjacency[i, j]] for i in 1:n for
    j in 1:n if Q.adjacency[i, j] > 0
  )
end
