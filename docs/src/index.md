# QuiverTools

## Introduction

QuiverTools is a software suite for treatment of quiver representations,
their roots, their moduli spaces and computations of several invariants of these.

QuiverTools is available as a Julia package and as a Sage library.

## Constructing quivers

You can build a quiver from an adjacency matrix, from a compact string, or with one
of the many built-in named constructors:

```julia-repl
julia> using QuiverTools

julia> Quiver([0 3; 0 0])                 # from an adjacency matrix
Quiver with adjacency matrix [0 3; 0 0]

julia> Quiver("1--2-3")                   # from a string: a hyphen run is the arrow count
Quiver with adjacency matrix [0 2 0; 0 0 1; 0 0 0]

julia> kronecker_quiver(3)                # one of many built-in named quivers
3-Kronecker quiver
```

See [Constructors](@ref) for the full catalogue.

## Installation

You can install `QuiverTools` as any other Julia package, by running

```julia-repl
julia> using Pkg; Pkg.add("QuiverTools")

```

## Acknowledgements

QuiverTools is developed by [P. Belmans](https://pbelmans.ncag.info/),
H. Franzen and [G. Petrella](https://www.giannipetrella.eu).

The Julia version is developed and maintained by
[G. Petrella](https://giannipetrella.eu).

G. Petrella was supported by the Luxembourg National Research Fund (FNR–17953441).
