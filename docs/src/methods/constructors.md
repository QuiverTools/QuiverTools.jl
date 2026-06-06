# Constructors

There are three ways to build a quiver in QuiverTools: from an adjacency matrix,
from a compact string description, or with one of the many built-in named constructors.

## From an adjacency matrix

A quiver is determined by its ``n \times n`` adjacency matrix, where the entry ``a_{ij}``
is the number of arrows ``i \to j``. Pass the matrix to [`Quiver`](@ref); an optional
second argument gives the quiver a name.

```julia-repl
julia> using QuiverTools

julia> Quiver([0 3; 0 0])
Quiver with adjacency matrix [0 3; 0 0]

julia> Quiver([0 3; 0 0], "my favourite quiver")
my favourite quiver
```

## From a string

For quick experiments it is often handy to type a quiver directly. [`Quiver`](@ref)
accepts a string that is a comma-separated list of chains `i-j-k-...`:

- a **run of `r` hyphens** between two vertices encodes `r` parallel arrows;
- a **chain** `i-j-k` is read left to right, producing arrows `i → j` and `j → k`;
- vertex labels may be **arbitrary tokens**; they are numbered `1, …, n` in order of
  first appearance.

```julia-repl
julia> Quiver("1--2-3")
Quiver with adjacency matrix [0 2 0; 0 0 1; 0 0 0]

julia> Quiver("a---b")
Quiver with adjacency matrix [0 3; 0 0]

julia> Quiver("1--2,1---3,2----3")
Quiver with adjacency matrix [0 2 3; 0 0 4; 0 0 0]
```

The second example shows that `Quiver("a---b")` and `kronecker_quiver(3)` describe the
same quiver: three arrows between two vertices.

## Named quivers

QuiverTools ships constructors for most quivers that show up in practice.

```@docs
kronecker_quiver
loop_quiver
jordan_quiver
subspace_quiver
generalized_subspace_quiver
thickened_subspace_quiver
three_vertex_quiver
cyclic_quiver
bipartite_quiver
dynkin_quiver
extended_dynkin_quiver
```

## Operations on quivers

New quivers can also be built from existing ones.

```@docs
opposite_quiver
double_quiver
disjoint_union
```
