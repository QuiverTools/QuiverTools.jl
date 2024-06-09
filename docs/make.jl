using Documenter
using QuiverTools

# memo for myself: when changing dependencies of QuiverTools,
# update the GLOBAL Julia environment to have the documenter work.

```@meta
CurrentModule = QuiverTools
```
DocMeta.setdocmeta!(QuiverTools, :DocTestSetup, :(using QuiverTools))

makedocs(
    sitename = "QuiverTools",
    authors = "Gianni Petrella",
    format = Documenter.HTML(),
    # format = Documenter.LaTeX(), # builds pdf, does not like the github Documenter action for now. Use only in local build.
    modules = [QuiverTools],
    pages = [   "QuiverTools" => "index.md", 
    "Tutorial" => "tutorial.md",
    "All methods" => "methods.md",
    "Benchmarks" => "benchmarks.md"]
    )

# deploydocs(
#     branch = "docs",
#     repo = "github.com/quiver-tools/QuiverTools.jl.git")