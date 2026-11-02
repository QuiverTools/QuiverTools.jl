# Run with julia --project benchmark/quasimap_spectrum.jl 3 4 5.
# Larger values of n can be supplied explicitly.
# The spectrum is encoded exactly by the factored characteristic polynomial
# of anticanonical multiplication at q=1, including algebraic multiplicities.
using QuiverTools, TOML, SHA, Dates
const QT = QuiverTools
const N = QT.Singular.Nemo
const T, t = N.polynomial_ring(N.QQ, "t")
const output = get(ENV, "QT_SPECTRUM_RESULTS",
  joinpath(@__DIR__, "raw", "quasimap_spectrum_results.toml"))
const requested_samples = parse(Int, get(ENV, "QT_SPECTRUM_SAMPLES", "3"))
const repeat_threshold = parse(Float64, get(ENV, "QT_SPECTRUM_REPEAT_THRESHOLD", "20"))
const cases = isempty(ARGS) ? collect(3:5) : parse.(Int, ARGS)
const results = Dict{String,Any}(
  "julia_version" => string(VERSION),
  "julia_threads" => Threads.nthreads(),
  "package_source" => pathof(QT),
  "quantum_source_sha256" =>
    bytes2hex(sha256(read(joinpath(dirname(pathof(QT)), "Quantum.jl")))),
  "started_utc" => string(now(UTC)),
  "operator" => "c1(T_M) multiplication at q=1",
  "dimension_vector" => [2, 3],
  "requested_samples" => requested_samples,
  "repeat_threshold_seconds" => repeat_threshold,
  "cases" => Dict{String,Any}[],
)
function save_results()
  mkpath(dirname(output))
  open(output * ".tmp", "w") do io
    TOML.print(io, results; sorted=true)
  end
  mv(output * ".tmp", output; force=true)
end
function progress(message)
  println(message)
  flush(stdout)
end
function spectrum(n; announce=true)
  QT.empty_all_caches!()
  GC.gc()
  started = time_ns()
  presentation = @timed begin
    M = QuiverModuliSpace(kronecker_quiver(n), [2, 3])
    A = quasimap_quantum_cohomology(M)
    M, A
  end
  M, A = presentation.value
  announce && progress("PHASE n=$n presentation=$(presentation.time)")
  multiplication = @timed begin
    x = quantum_tautological_classes(A)
    c1 = -sum(canonical_stability(M.Q, M.d)[i] * x[i][1] for i in eachindex(M.d))
    quantum_multiplication_matrix(A, c1, [1])
  end
  matrix, basis, _ = multiplication.value
  announce && progress("PHASE n=$n matrix=$(multiplication.time) rank=$(length(basis))")
  characteristic = @timed N.charpoly(T, matrix)
  announce && progress("PHASE n=$n charpoly=$(characteristic.time)")
  factorization = @timed N.factor(characteristic.value)
  total = (time_ns() - started) / 1.0e9
  announce && progress("PHASE n=$n factorization=$(factorization.time)")
  @assert gw_comparison(A) == :equal
  @assert N.degree(characteristic.value) == length(basis)
  @assert prod(f^e for (f, e) in factorization.value; init=one(T)) == characteristic.value
  row = Dict{String,Any}(
    "n" => n, "dimension" => 6n - 12, "rank" => length(basis),
    "presentation_seconds" => presentation.time,
    "matrix_seconds" => multiplication.time,
    "charpoly_seconds" => characteristic.time,
    "factorization_seconds" => factorization.time,
    "total_seconds" => total,
    "julia_allocated_bytes" => sum(r.bytes for r in
                    (presentation, multiplication, characteristic, factorization)),
    "factor_degrees" => [N.degree(f) for (f, _) in factorization.value],
    "factor_multiplicities" => [e for (_, e) in factorization.value],
    "zero_eigenvalue_multiplicity" =>
      sum(e for (f, e) in factorization.value if iszero(f(0)); init=0),
  )
  return row, characteristic.value, factorization.value, M
end
save_results()
progress("WARMUP_START")
spectrum(3; announce=false)
progress("WARMUP_DONE")
for n in cases
  row = Dict{String,Any}("n" => n, "status" => "running", "samples" => Dict{String,Any}[])
  push!(results["cases"], row)
  save_results()
  progress("CASE_START n=$n")
  first, polynomial, factors, M = spectrum(n)
  push!(row["samples"], first)
  row["rank"] = first["rank"]
  row["dimension"] = first["dimension"]
  row["characteristic_polynomial"] = string(polynomial)
  row["factorization"] = string(factors)
  save_results()
  samples = first["total_seconds"] <= repeat_threshold ? requested_samples : 1
  for trial in 2:samples
    sample, check_polynomial, _, _ = spectrum(n)
    @assert check_polynomial == polynomial
    push!(row["samples"], sample)
    save_results()
  end
  @assert sum(betti_numbers(M)) == row["rank"]
  row["status"] = "complete"
  save_results()
  progress("CASE_DONE n=$n rank=$(row["rank"]) samples=$(length(row["samples"]))")
end
results["finished_utc"] = string(now(UTC))
save_results()
progress("BENCHMARK_DONE")
