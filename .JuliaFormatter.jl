using JuliaFormatter

files = [
  "src/Bundles.jl",
  "src/Chow.jl",
  "src/Constructors.jl",
  "src/Hodge.jl",
  "src/Misc.jl",
  "src/Moduli.jl",
  "src/Quivers.jl",
  "src/QuiverTools.jl",
  "src/RepresentationTheory.jl",
  "src/Stability.jl",
  "src/Teleman.jl",
  "src/Types.jl",
  "docs/docs.jl",
  "benchmark/benchmarks.jl",
  "test/runtests.jl",
  ".JuliaFormatter.jl",
]

for file in files
  if !format(file)
    exit(1)
  end
end
exit(0)
