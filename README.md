# QuiverTools

[![tests](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Runtests.yml)
[![docs](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Documenter.yml)
[![style](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/QuiverTools/QuiverTools.jl/actions/workflows/Format.yml)

QuiverTools is an open source Julia package for working
with moduli spaces of quiver representations.

This branch contains an unstable version of QuiverTools.
It works exactly the same way as the stable,
but some internal methods rely on unstable dependencies
and outperform the stable release.

## Installation

To install this branch of QuiverTools,
one must checkout this branch of QuiverTools
and install the dependencies manually
in the desired environment.
Some familiarity in handling Julia environments is useful.

First, activate a new environment in a folder of choice by running

```bash
julia --project=.
```

Then, add the dependencies manually by running

```julia

julia> using Pkg

julia> Pkg.add(url="https://github.com/pseudoeffective/BumplessPipeDreams.jl");

julia> Pkg.add(url="https://github.com/pseudoeffective/SchubertPolynomials.jl");

julia> Pkg.add(url="https://github.com/QuiverTools/QuiverTools.jl", rev="schubert-polynomials");
```

Once this is done, QuiverTools will be ready to use in the selected environment.

Note: this conflicts with the stable version of QuiverTools.
This means that both can not be used in the same environment.

## Documentation

The documentation for QuiverTools is available
[here](https://QuiverTools.github.io/QuiverTools.jl/dev/).

## Acknowledgements

QuiverTools is developed by [Pieter Belmans](https://pbelmans.ncag.info/),
Hans Franzen and
[Gianni Petrella](https://www.giannipetrella.eu/).

Gianni Petrella is supported by the Luxembourg National Research Fund (FNR-7953441).
