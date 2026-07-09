#####################################################
# Methods that deal with quivers and their properties
#####################################################

function deglex_key(Q::Quiver, e::AbstractVector{Int})
  b = maximum(e) + 1
  n = n_vertices(Q)

  return Int(sum(e[i] * b^(n - i) for i in 1:length(e)) + sum(e) * b^n)
end

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
  return Q.adjacency + transpose(Q.adjacency) - diagonal(Q.adjacency)
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
n_vertices(Q::Quiver) = size(Q.adjacency, 1)

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
n_arrows(Q::Quiver) = sum(Q.adjacency; init=0)

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
  paths = sum(underlying_graph(Q)^k for k in 0:(n_vertices(Q) - 1))
  return all(p -> p > 0, paths)
end

"""
    strongly_connected_components(Q::Quiver)

Compute the strongly connected components of `Q`.

Two vertices belong to the same strongly connected component if and only if
they are connected by paths in both directions. The reachability relation is
computed as the reflexive-transitive closure of the adjacency relation, using the
Floyd--Warshall algorithm in its original, Boolean, form
[[Warshall](https://doi.org/10.1145/321105.321107)]; its ``O(n^3)`` running time
is not an issue for the quivers we consider.

# Input

- `Q::Quiver`: a quiver.

# Output

- a list of the strongly connected components, each given as the list of its vertices.

# Examples

```jldoctest
julia> strongly_connected_components(cyclic_quiver(3))
1-element Vector{Vector{Int64}}:
 [1, 2, 3]

julia> strongly_connected_components(kronecker_quiver(3))
2-element Vector{Vector{Int64}}:
 [1]
 [2]

julia> strongly_connected_components(Quiver("1-2,2-1,2-3"))
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3]
```
"""
function strongly_connected_components(Q::Quiver)
  n = n_vertices(Q)
  # reflexive-transitive closure by Floyd--Warshall [doi:10.1145/321105.321107]
  reachable = [i == j || Q.adjacency[i, j] > 0 for i in 1:n, j in 1:n]
  for k in 1:n, i in 1:n, j in 1:n
    reachable[i, j] |= reachable[i, k] && reachable[k, j]
  end
  return unique([findall(j -> reachable[i, j] && reachable[j, i], 1:n) for i in 1:n])
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
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j]; init=0)

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
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :]; init=0)

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
[[Proposition 1.6, MR1035222](https://mathscinet.ams.org/mathscinet/relay-station?mr=1035222)]
to be

```math
\\begin{aligned}
\\mathrm{HH}^{0}(Q) &= k,\\\\
\\mathrm{HH}^{1}(Q) &= 1 - n + \\sum_{\\alpha \\in Q_1} \\# \\{ \\text{paths from } s(\\alpha) \\text{ to } t(\\alpha) \\},\\\\
\\mathrm{HH}^{n}(Q) &= 0 \\text{ for } n > 1.
\\end{aligned}
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

julia> first_hochschild_cohomology(Quiver([0;;]))
0
```
"""
function first_hochschild_cohomology(Q::Quiver)
  !is_acyclic(Q) && throw(ArgumentError("The quiver must be acyclic."))
  n = n_vertices(Q)
  n == 1 && return 0
  !is_connected(Q) && throw(ArgumentError("The quiver must be connected."))

  path_matrix = Q.adjacency
  path_matrix += sum((Q.adjacency)^k for k in 2:(n - 1); init=zeros(Int, n, n))

  return 1 - n + sum(Q.adjacency .* path_matrix; init=0)
end
