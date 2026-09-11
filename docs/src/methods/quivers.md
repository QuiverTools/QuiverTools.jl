# Quivers

Quivers are represented via their adjacency matrix.
Vertices are numbered from ``1`` to ``n``.

```@docs
Quiver
```

For all built-in quivers and the different ways to construct them,
see [Constructors](@ref).

## Properties

QuiverTools offers various methods to study some graph-theoretical properties of quivers.

```@docs
n_vertices
n_arrows
arrows
indegree
outdegree
is_acyclic
is_connected
is_sink
is_source
underlying_graph
```

## Visualization

`to_dot` renders a quiver as a [Graphviz](https://graphviz.org) DOT string,
drawing parallel arrows and loops as themselves. Running `using GraphViz`
alongside QuiverTools additionally makes a `Quiver` display as an SVG drawing in
notebooks (Pluto, IJulia) and VS Code; from the REPL, `draw` opens that drawing
in the system viewer.

```@docs
to_dot
draw
```
