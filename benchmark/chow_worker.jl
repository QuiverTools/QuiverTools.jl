# Run via run_chow.py. Each process measures exactly one case and algorithm.
using QuiverTools, TOML
const QT = QuiverTools
name, algorithm, destination = ARGS
algorithm = Symbol(algorithm)
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
# Identical generic warmup in every process, including both algorithms.
for a in (:weyl, :experimental_parabolic)
  chow_ring(kronecker_quiver(3), [1, 1]; algorithm=a, unsafe=true)
end
result = Dict{String,Any}("case" => name, "algorithm" => string(algorithm),
  "julia" => string(VERSION), "singular" => string(Base.pkgversion(QT.Singular)),
  "cpu" => first(Sys.cpu_info()).model, "threads" => Threads.nthreads(),
  "adjacency" => [collect(row) for row in eachrow(Q.adjacency)], "d" => d,
  "theta" => collect(canonical_stability(Q, d)), "samples" => Any[])
function save()
  open(destination * ".tmp", "w") do io
    TOML.print(io, result)
  end
  mv(destination * ".tmp", destination; force=true)
end
save()
for trial in 0:3
  QT.empty_all_caches!()
  GC.gc()
  write(destination * ".phase", "trial=$trial")
  timed = @timed if endswith(name, "_chern")
    M = QuiverModuliSpace(Q, d)
    chow_ring(M; algorithm=algorithm, unsafe=true)
    chern_numbers(M; unsafe=true)
    M
  else
    chow_ring(Q, d; algorithm=algorithm, unsafe=true)
  end
  write(destination * ".phase", "verify=$trial")
  push!(
    result["samples"],
    Dict("trial" => trial, "seconds" => timed.time,
      "julia_bytes" => timed.bytes, "gc_seconds" => timed.gctime),
  )
  save()
  ring = endswith(name, "_chern") ? chow_ring(timed.value) : timed.value[1][1]
  result["samples"][end]["chow_total_dimension"] = length(
    QT.gens(QT.Singular.kbase(QT.quotient_ideal(ring)))
  )
  save()
end
write(destination * ".phase", "complete")
