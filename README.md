# QuiverTools

[![tests](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Runtests.yml)
[![docs](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Documenter.yml)
[![style](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Format.yml)

QuiverTools is an open source Julia package for working
with moduli spaces of quiver representations.

## Installation

To install QuiverTools,
run the following commands in the Julia REPL:

```julia
julia> using Pkg; Pkg.add("QuiverTools");

```

## Quick start

You can build a quiver from an adjacency matrix, from a compact string, or with one
of the many built-in named constructors:

```julia
using QuiverTools

Quiver([0 3; 0 0])        # from an adjacency matrix
Quiver("1--2-3")          # from a string: a hyphen run is the number of arrows
kronecker_quiver(3)       # one of many built-in named quivers
```

See the [documentation](https://julia.quiver.tools)
for the full catalogue of quiver constructors.

## Documentation

The documentation for QuiverTools is available
[here](https://QuiverTools.github.io/QuiverTools.jl/).

## Acknowledgements

QuiverTools is developed by
[Pieter Belmans](https://pbelmans.ncag.info/),
Hans Franzen and
[Gianni Petrella](https://www.giannipetrella.eu/).

Gianni Petrella is supported by the Luxembourg National Research Fund (FNR-7953441).
