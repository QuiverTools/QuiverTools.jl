using Documenter

using Pkg
Pkg.activate(@__DIR__)
# As long as it is not registered, this is nice, in general it locally always
# renders docs of the current version checked out in this repo.
Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))

using QuiverTools

```@meta
CurrentModule = QuiverTools
```
DocMeta.setdocmeta!(QuiverTools, :DocTestSetup, :(using QuiverTools))

makedocs(;
  sitename="QuiverTools",
  authors="Gianni Petrella",
  doctest=false,
  modules=[QuiverTools],
  format=Documenter.HTML(),

  # TODO this should detect whether we are working locally or not
  warnonly=true,

  # Use the following two parameters for local pdf build.
  # format=Documenter.LaTeX(),
  # remotes=nothing,
  pages=[
    "QuiverTools" => "index.md",
    "Tutorial" => "tutorial.md",
    "All methods" => [
      "Quivers" => "methods/quivers.md",
      "Quiver moduli" => "methods/quiver-moduli.md",
      "Representation theory" => "methods/representation-theory.md",
      "Teleman quantization" => "methods/teleman-quantization.md",
      "Chow rings" => "methods/chow-rings.md",
      "Walls and chamber decompositions" => "methods/walls-and-chambers.md",
    ],
    "Benchmarks" => "benchmarks.md"],
)

deploydocs(;
  branch="docs",
  repo="github.com/QuiverTools/QuiverTools.jl.git",
  cname="julia.quiver.tools",
  versions=["v#.#.#"],
)
