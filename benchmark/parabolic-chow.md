# Experimental parabolic Chow benchmark

**LLM-suggested experiment; careful independent checking is required.**
These measurements do not establish correctness or a uniform performance improvement.
The proof and geometric hypotheses are in
[the experimental documentation](../docs/src/experimental-chow.md).

Measured 2026-11-14 on Apple M4, arm64, Julia 1.13.0,
Singular.jl 0.28.15, one Julia thread.
The reference is the unchanged `:weyl` backend from
`f1a12e766264204cca9d212c68f5193a5949da90`.
The experimental backend uses `adaptive=true`.
No tests or other Julia jobs were run concurrently with these serial measurements.

Each case/algorithm runs in a fresh process, with the same generic warmup of both
algorithms on the 3-arrow Kronecker $(1,1)$ model.
The first timed call includes remaining case-specific compilation.
The warm result is the median of three further constructions in the same process.
Before each sample, Memoization function-result caches are cleared and a full Julia GC is run.
Compiled code and Singular's polynomial parent-ring caches remain available.
Each call recomputes the relations and Gröbner bases;
the experimental algorithm also repeats its Betti computation.
The timings include any GC occurring during the computation.
Startup, generic warmup, cache clearing, the preceding GC, and subsequent dimension
checks are excluded.
`unsafe=true` is used for both algorithms, excluding the public geometric guards;
the experimental Betti certificate is included.
The two `_chern` cases include a fresh moduli object, ring construction, and all
Chern numbers, using the upstream `chern_numbers` implementation.
All other rows measure ring construction alone.

| Model / computation | Weyl first (s) | Experimental first (s) | Weyl warm median (ms) | Experimental warm median (ms) |
|---|---:|---:|---:|---:|
| $\mathbb P^3$ | 0.000178416 | 0.000364667 | 0.170916 | 0.294417 |
| $\mathrm{Gr}(2,4)$ | 0.118713 | 0.352805 | 0.202708 | 0.533083 |
| $\mathrm{Gr}(3,6)$ | 0.13229 | 0.382814 | 0.414625 | 1.73358 |
| 3-Kronecker, $(2,3)$ | 0.00220717 | 0.38357 | 1.95542 | 6.19208 |
| Seven-subspace fourfold | 1.00519 | 1.8305 | 8.88167 | 53.883 |
| Chain, $(5,5,2)$ | >120 (timeout) | 25.4782 | — | 23554.5 |
| Fork, $(5,5,2)$ | >120 (timeout) | 31.0579 | — | 29085.7 |
| 3-Kronecker, $(3,4)$ | 2.95846 | 1.86351 | 2871.04 | 1423.19 |
| 3-Kronecker, $(3,5)$ | >120 (timeout) | 4.80506 | — | 4471.92 |
| 3-Kronecker $(2,3)$, all Chern numbers | 0.033011 | 0.474405 | 3.56692 | 7.59263 |
| Seven-subspace fourfold, all Chern numbers | 1.29513 | 2.40271 | 18.7226 | 65.106 |

For **Chain, $(5,5,2)$**, the reference was still inside the first timed construction at the 120-second deadline. The first-call speedup is therefore **greater than 4.7×**; the exact speedup is unknown. The experimental warm constructions allocated about 21.09 GB cumulatively in Julia per call.

For **Fork, $(5,5,2)$**, the reference was still inside the first timed construction at the 120-second deadline. The first-call speedup is therefore **greater than 3.86×**; the exact speedup is unknown. The experimental warm constructions allocated about 26.76 GB cumulatively in Julia per call.

For **3-Kronecker, $(3,5)$**, the reference was still inside the first timed construction at the 120-second deadline. The first-call speedup is therefore **greater than 24.97×**; the exact speedup is unknown. The experimental warm constructions allocated about 3.39 GB cumulatively in Julia per call.

These are cumulative allocation counts, not peak resident memory.
Native Singular allocations and peak memory were not measured.
The small reference cases regress; their reference runtimes are already very small,
and the certificate and truncation add overhead.
The three-sample medians are exploratory measurements, not a statistical performance study.
The large-case timeout bounds compare first calls only; a censored first call is not
compared with an experimental warm median.
Basis-rank reductions must not be read as runtime factors.
The experiment combines several changes rather than isolating their individual effects.

The projective space, Grassmannian, favourite Kronecker sixfold, and seven-subspace
fourfold include examples from the existing Chow/Bundles doctests.
The additional Kronecker and rank $(5,5,2)$ examples extend that coverage.
The exact quiver matrices, dimension vectors, canonical stability parameters,
per-sample seconds, GC time, allocation counts, quotient dimensions, and timeout
statuses are in [the machine-readable results](parabolic-chow-results.json).
That file also records the resolved dependency versions and source tree hashes.
For both rank $(5,5,2)$ examples, tests separately check the graded dimensions
`[1,2,4,4,4,4,2,1]` and certified stopping at Schur degree one.
The old implementation did not finish on the two rank $(5,5,2)$ examples or
the 3-Kronecker $(3,5)$ example, so no whole-ideal comparison with it is claimed there.
`benchmark/compare_chow.jl` checks exact equality for all six distinct models
on which the reference completed, outside the timing runs.
The smaller tests compare complete ideals and known Chern numbers.

## Reproduce

From this branch's repository root, instantiate the package environment with Julia
1.13.0 and Singular.jl 0.28.15 (the checked-in results record these versions).
For example, use a disposable checkout and run
`julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.add(PackageSpec(name="Singular", version="0.28.15"))'`.
The `Pkg.add` is only to reproduce the measured version; do not commit its Project changes.
Python 3.11 or newer is needed for the TOML reader in the summarizer.

```sh
python3 benchmark/run_chow.py --julia /path/to/julia --timeout 120
python3 benchmark/summarize_chow.py
/path/to/julia --startup-file=no --project=. benchmark/compare_chow.jl
```

Generated logs and TOML samples are ignored in `benchmark/raw/`.
`--cases` selects models, and `--algorithms` selects backends.
Each construction has a 120-second deadline, independently of the startup/warmup
and verification phases.
A timed-out process is killed and reaped before the next benchmark starts.
The driver preserves completed samples if a later phase times out.
The checked-in JSON contains the final samples, including separately rerun censored
references with explicit construction/verification phase markers.
