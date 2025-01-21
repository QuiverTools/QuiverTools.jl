# Quivers

```@meta
CurrentModule = QuiverTools
```

Quivers are represented via their adjacency matrix.
Vertices are numbered from ``1`` to ``n``.

## Constructors

QuiverTools implements constructors for most known quivers.

```@docs
kronecker_quiver
loop_quiver
subspace_quiver
three_vertex_quiver
cyclic_quiver
bipartite_quiver
opposite_quiver
double_quiver
dynkin_quiver
```

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
