using JuliaFormatter

files = [
  "src/QuiverTools.jl",
  "src/Types.jl",
  "src/Teleman.jl",
  "src/Moduli.jl",
  "src/Constructors.jl",
  "docs/make.jl",
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
