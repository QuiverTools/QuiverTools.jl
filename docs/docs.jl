using Documenter, QuiverTools
# To ensure that the version of the documentation being compiled
# is the latest committed one,
# the local installation of QuiverTools should be the one developed at path="../."

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
      "Constructors" => "methods/constructors.md",
      "Quiver moduli" => "methods/quiver-moduli.md",
      "Representation theory" => "methods/representation-theory.md",
      "Teleman quantization" => "methods/teleman-quantization.md",
      "Chow rings" => "methods/chow-rings.md",
      "Walls and chamber decompositions" => "methods/walls-and-chambers.md",
    ],
    "Internal methods" => "methods.md",
    "Benchmarks" => "benchmarks.md"],
)

deploydocs(;
  branch="docs",
  repo="github.com/QuiverTools/QuiverTools.jl.git",
  cname="julia.quiver.tools",
  versions=["stable" => "v^", "dev" => "dev", "v#.#.#"],
)

## READ BEFORE DEPLOYING
#
# Since each time I try do deploy the documentation it fails,
# here are the exact steps to deploy the docs for a new release.
#
# 0. Make a commit with the version bump in Project.toml
# 1. Push, and register the bump commit with @JuliaRegistrator register, with the release notes.
# 2. Once the registry PR merges, TagBot creates the version tag and release. This can take time!
# 3. From the Actions tab, run the Documenter workflow by hand **from the v2.0.0 tag**
#    (“Run workflow” → “Use workflow from” → Tags → v2.0.0).
