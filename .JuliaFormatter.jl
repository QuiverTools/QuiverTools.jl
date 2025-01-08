using JuliaFormatter

files = [
    "src/QuiverTools.jl",
    "src/QuiverTools-types.jl",
    "src/teleman.jl",
    "src/moduli.jl",
    "src/constructors.jl",
    "docs/make.jl",
    "benchmark/benchmarks.jl",
    "test/runtests.jl",
    ]

for file in files
  if !format(file) exit(1) end
end
exit(0)
