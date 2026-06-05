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
