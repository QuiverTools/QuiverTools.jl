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

clean = map(format, files)
all(clean) && exit(0)

# print non-formatted files
for x in files[.!(clean .== true)]
  println("Not formatted: ", x)
end
exit(1)
