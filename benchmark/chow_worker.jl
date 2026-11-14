# Run by run_chow.py in a fresh process for one case and one algorithm.
using QuiverTools, TOML

"""
    save_chow_benchmark(destination::AbstractString, result::AbstractDict)

Save completed samples atomically, so a timeout leaves readable results.
"""
function save_chow_benchmark(destination::AbstractString, result::AbstractDict)
  open(destination * ".tmp", "w") do io
    TOML.print(io, result)
  end
  mv(destination * ".tmp", destination; force=true)
  return nothing
end

"""
    benchmark_chow(name::AbstractString, algorithm::Symbol, destination::AbstractString)

Measure one first call and three warm calls. Clear Memoization results and run
GC before each sample; compiled code and Singular's parent-ring caches persist.
The `_chern` cases also construct a fresh moduli object and compute all its
Chern numbers. Verify quotient dimensions outside the timed section.

Write phase markers for the driver's per-phase timeout, and save each sample
before starting its verification. No other Julia jobs should run concurrently.
"""
function benchmark_chow(
  name::AbstractString, algorithm::Symbol, destination::AbstractString
)
  cases = Dict(
    "P3" => (kronecker_quiver(4), [1, 1]),
    "Gr24" => (kronecker_quiver(4), [1, 2]),
    "Gr36" => (kronecker_quiver(6), [1, 3]),
    "K3_23" => (kronecker_quiver(3), [2, 3]),
    "K3_34" => (kronecker_quiver(3), [3, 4]),
    "K3_35" => (kronecker_quiver(3), [3, 5]),
    "subspace7" => (subspace_quiver(7), [ones(Int, 7); 2]),
    "K3_23_chern" => (kronecker_quiver(3), [2, 3]),
    "subspace7_chern" => (subspace_quiver(7), [ones(Int, 7); 2]),
    "chain_552" => (Quiver([0 2 0; 0 0 1; 0 0 0]), [5, 5, 2]),
    "fork_552" => (Quiver([0 2 1; 0 0 0; 0 0 0]), [5, 5, 2]),
  )
  Q, d = cases[name]
  include_chern_numbers = endswith(name, "_chern")

  # Identical generic warmup of both algorithms in every worker.
  for backend in (:weyl, :experimental_parabolic)
    chow_ring(kronecker_quiver(3), [1, 1]; algorithm=backend, unsafe=true)
  end
  result = Dict{String,Any}(
    "case" => name,
    "algorithm" => string(algorithm),
    "julia" => string(VERSION),
    "singular" => string(Base.pkgversion(QuiverTools.Singular)),
    "cpu" => first(Sys.cpu_info()).model,
    "threads" => Threads.nthreads(),
    "adjacency" => collect.(eachrow(Q.adjacency)),
    "d" => d,
    "theta" => collect(canonical_stability(Q, d)),
    "samples" => Dict{String,Any}[],
  )
  save_chow_benchmark(destination, result)
  for trial in 0:3
    QuiverTools.empty_all_caches!()
    GC.gc()
    write(destination * ".phase", "trial=$trial")
    timed = @timed if include_chern_numbers
      M = QuiverModuliSpace(Q, d)
      chow_ring(M; algorithm=algorithm, unsafe=true)
      chern_numbers(M; unsafe=true)
      M
    else
      chow_ring(Q, d; algorithm=algorithm, unsafe=true)[1][1]
    end
    write(destination * ".phase", "verify=$trial")
    sample = Dict{String,Any}(
      "trial" => trial, "seconds" => timed.time,
      "julia_bytes" => timed.bytes, "gc_seconds" => timed.gctime,
    )
    push!(result["samples"], sample)
    save_chow_benchmark(destination, result)
    ring = include_chern_numbers ? chow_ring(timed.value) : timed.value
    sample["chow_total_dimension"] = length(
      QuiverTools.gens(QuiverTools.Singular.kbase(QuiverTools.quotient_ideal(ring)))
    )
    save_chow_benchmark(destination, result)
  end
  write(destination * ".phase", "complete")
  return nothing
end

length(ARGS) == 3 || throw(ArgumentError("expected case, algorithm, and output path"))
benchmark_chow(ARGS[1], Symbol(ARGS[2]), ARGS[3])
