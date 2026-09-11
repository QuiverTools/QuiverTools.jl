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

# Compass ports used to spread self-loops around a vertex: Graphviz stacks
# multiple loops on the same side otherwise, which reads as a single blob.
const _LOOP_PORTS = ("n", "s", "e", "w", "ne", "sw", "se", "nw")

function _escape_dot_string(value::AbstractString)
  return replace(
    value,
    '\\' => "\\\\",
    '"' => "\\\"",
    '\n' => "\\n",
    '\r' => "\\r",
  )
end

"""
    to_dot(Q::Quiver)

Return a Graphviz DOT description of `Q` as a `String`.

Every arrow becomes its own `i -> j` line, so parallel arrows and loops (a
quiver's defining features) are drawn as themselves; self-loops are distributed
around their vertex with compass ports. Isolated vertices are emitted explicitly.
Quiver names are preserved as UTF-8 labels; DOT syntax characters and physical
line breaks are escaped.

Feed the result to Graphviz, e.g. `using GraphViz; GraphViz.Graph(to_dot(Q))`,
which also makes `Q` render as SVG in notebooks and VS Code once GraphViz is
loaded. With the DOT string in hand you can equally pipe it to the `dot`
command line tool.

# Examples

```jldoctest
julia> print(to_dot(kronecker_quiver(2)))
digraph {
  label="2-Kronecker quiver";
  node [shape=circle];
  1;
  2;
  1 -> 2;
  1 -> 2;
}

julia> print(to_dot(loop_quiver(2)))
digraph {
  label="2-loop quiver";
  node [shape=circle];
  1;
  1:n -> 1:n;
  1:s -> 1:s;
}
```
"""
function to_dot(Q::Quiver)
  n = n_vertices(Q)
  io = IOBuffer()
  println(io, "digraph {")
  isempty(Q.name) || println(io, "  label=\"", _escape_dot_string(Q.name), "\";")
  println(io, "  node [shape=circle];")
  for i in 1:n
    println(io, "  ", i, ";")
  end
  for i in 1:n, j in 1:n
    if i == j
      for t in 1:Q.adjacency[i, i]
        p = _LOOP_PORTS[mod1(t, length(_LOOP_PORTS))]
        println(io, "  ", i, ":", p, " -> ", i, ":", p, ";")
      end
    else
      for _ in 1:Q.adjacency[i, j]
        println(io, "  ", i, " -> ", j, ";")
      end
    end
  end
  println(io, "}")
  return String(take!(io))
end

"""
    draw(Q::Quiver)

Render `Q` with Graphviz and open the drawing in the system's default viewer
(browser or image application). Writes a temporary SVG file and returns its path.

Requires GraphViz: run `using GraphViz` to enable this method. In notebooks and
VS Code you do not need `draw`, evaluating `Q` displays the drawing inline.

# Example

```julia
using QuiverTools, GraphViz

draw(loop_quiver(2))   # opens the drawing; returns the temporary file path
```
"""
function draw end

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
