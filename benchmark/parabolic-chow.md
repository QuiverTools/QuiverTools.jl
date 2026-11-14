# Experimental parabolic Chow benchmark

**LLM-suggested experiment; the proof and code need careful checking.**
The reduction and hypotheses are in [the documentation](../docs/src/experimental-chow.md).

Measured 2026-11-14 on Apple M4, Julia 1.13.0, Singular.jl 0.28.15, one Julia thread.
The Chow implementation is commit `771a2ed131c3a446bfdc88d32f46e4b0a418d433`; the worker adds the cutoff models listed in the raw data.
Both algorithms ran serially in fresh processes, with the same generic warmup on 3-Kronecker $(1,1)$.
The first timed call includes remaining case-specific compilation.
Each warm time is the median of three further constructions, with Memoization results cleared and full GC before every sample.
Compiled code and Singular's polynomial parent-ring caches persist; relations, Gröbner bases and experimental Betti numbers are recomputed.
Startup, warmup, cache clearing, preceding GC and subsequent quotient-dimension checks are outside the timer.
GC during construction is included. Both algorithms use `unsafe=true`; parabolic timings include the Betti check.
The two Chern rows include a fresh moduli object, ring construction and all Chern numbers.

Each phase has a three-minute maximum and a 20 GiB process-tree memory limit.
The additional cutoff runs have a 30-second phase limit; shared rows reuse the longer main-table runs.
For 5-Kronecker $(2,5)$, an additional probe allowed 180 seconds for the first call and 30 seconds for the next warm call. Its warm bound refers to one construction, not a three-sample median.
On macOS the watchdog measures both resident memory and physical footprint and kills the process group if either exceeds the limit.

## Main timings

Times are in seconds. A bound in the first-call column is a timeout, not an exact time.

| Model / computation | Weyl first | Parabolic first | Weyl warm | Parabolic warm | Weyl / parabolic peak RSS (GiB) |
|:--|--:|--:|--:|--:|--:|
| $\mathbb P^3$ | 0.00022375 | 0.000369458 | 0.000163583 | 0.0002535 | 0.697 / 0.702 |
| $\mathrm{Gr}(2,4)$ | 0.123908 | 0.182595 | 0.000201791 | 0.000420459 | 0.704 / 0.705 |
| $\mathrm{Gr}(3,6)$ | 0.126041 | 0.186369 | 0.000411333 | 0.00121013 | 0.706 / 0.709 |
| 3-Kronecker, $(2,3)$ | 0.00285246 | 0.199511 | 0.00197238 | 0.00373467 | 0.703 / 0.705 |
| Seven-subspace fourfold | 1.00601 | 1.50545 | 0.00861092 | 0.0504744 | 0.702 / 0.751 |
| Chain, $(5,5,2)$ | >180 | 12.569 | — | 10.3405 | 1.254 / 1.769 |
| Fork, $(5,5,2)$ | >180 | 14.4348 | — | 13.6319 | 1.175 / 2.049 |
| 3-Kronecker, $(3,4)$ | 2.98925 | 1.17072 | 2.89469 | 0.964042 | 1.237 / 0.919 |
| 3-Kronecker, $(3,5)$ | 126.043 | 2.22345 | 127.676 | 1.97867 | 1.817 / 1.105 |
| 3-Kronecker $(2,3)$, all Chern numbers | 0.0157211 | 0.202874 | 0.00344713 | 0.00543358 | 0.697 / 0.700 |
| Seven-subspace fourfold, all Chern numbers | 1.05122 | 1.63374 | 0.0186272 | 0.0602035 | 0.726 / 0.763 |

## Effect of the cleanup

The earlier measurements used the same machine, dependency versions and three-warm-sample convention.
This comparison measures the combined cleanup; it does not isolate the effect of individual edits.
The benchmark worker was also moved into a function during cleanup, which can affect overhead on the smallest rows.

| Model | Previous parabolic warm (s) | Current parabolic warm (s) | Improvement | Current cumulative Julia allocation (GB) |
|:--|--:|--:|--:|--:|
| Chain, $(5,5,2)$ | 23.5545 | 10.3405 | 2.278× | 8.629 |
| Fork, $(5,5,2)$ | 29.0857 | 13.6319 | 2.134× | 10.838 |
| 3-Kronecker, $(3,4)$ | 1.42319 | 0.964042 | 1.476× | 0.864 |
| 3-Kronecker, $(3,5)$ | 4.47192 | 1.97867 | 2.260× | 1.819 |
| 3-Kronecker, $(2,3)$ | 0.00619208 | 0.00373467 | 1.658× | 0.004 |
| Seven-subspace fourfold | 0.053883 | 0.0504744 | 1.068× | 0.045 |

For Chain, $(5,5,2)$, Weyl exceeded 180 seconds on its first call; the first-call speedup is therefore greater than 14.32×.

For Fork, $(5,5,2)$, Weyl exceeded 180 seconds on its first call; the first-call speedup is therefore greater than 12.46×.

Across all runs, maximum resident memory was 3.250 GiB and maximum physical footprint was 3.122 GiB.
Julia allocation counts are cumulative and exclude native Singular allocations; they are not peak memory.

## Cutoff assessment

Let $W=\prod_i d_i!$, let $F$ count the minimal forbidden vectors, and let $B$ count the minimal weighted monomials above the dimension in the ambient Chern ring, before imposing the normalization relation.
The proposed score is $S=FW^2/B$, with parabolic chosen when $S\geq32$.
The $W^2$ term reflects the reference basis enumeration and explicit alternating sums; $B$ approximates the truncation overhead in the parabolic method.
This is an empirical proxy for the current implementations, not a complexity theorem.

The original 29-model cutoff sample was rerun. There are 20 complete warm comparisons and 20 outside a 10% tie window.
Only complete, three-sample warm medians enter the threshold-fit table.
First-call timeouts provide first-call bounds; the targeted warm timeout is shown separately. Rows with both algorithms censored remain unresolved.

The largest score among clear warm Weyl wins is 4.65267; the smallest among clear warm parabolic wins is 25.0435.

The first-call column below is assessed separately, using completed first calls outside the same 10% tie window.

| Threshold | Clear warm mistakes | Worst warm slowdown | Clear first-call mistakes |
|--:|:--|--:|:--|
| 16 | none | 1.000× | `K5_14`, `K5_25`, `triangle_211_223` |
| 20 | none | 1.000× | `K5_14`, `triangle_211_223` |
| 24 | none | 1.000× | `K5_14`, `triangle_211_223` |
| 32 | `K5_14` | 1.426× | `triangle_211_223` |

At score 25.0435, 5-Kronecker $(1,4)$ takes 1.473 ms warm with Weyl and 1.033 ms with parabolic.
Its first timed calls go the other way: 0.0875 s with Weyl and 0.1910 s with parabolic.
A threshold for warmed-up calculations need not give the best first-call choice under this particular generic warmup.

| Model | Score | Weyl warm (s) | Parabolic warm (s) | Comparison |
|:--|--:|--:|--:|:--|
| 4-Kronecker, $(1,3)$ | 2.76923 | 0.000335084 | 0.000606625 | Weyl |
| 5-Kronecker, $(1,4)$ | 25.0435 | 0.00147342 | 0.00103342 | parabolic |
| 6-Kronecker, $(1,5)$ | 400 | 0.03449 | 0.00228137 | parabolic |
| 7-Kronecker, $(1,6)$ | 9094.74 | 2.99085 | 0.00576529 | parabolic |
| 8-Kronecker, $(1,4)$ | 1.29148 | 0.00316379 | 0.00950325 | Weyl |
| 10-Kronecker, $(1,5)$ | 4.65267 | 0.0868146 | 0.195966 | Weyl |
| 3-Kronecker, $(2,3)$ | 4.5 | 0.00197238 | 0.00373467 | Weyl |
| 4-Kronecker, $(2,3)$ | 0.870091 | 0.00370379 | 0.01898 | Weyl |
| 5-Kronecker, $(2,3)$ | 0.269916 | 0.00762258 | 0.113349 | Weyl |
| 7-Kronecker, $(2,3)$ | 0.0520796 | 0.0476772 | 10.5861 | Weyl |
| 3-Kronecker, $(2,5)$ | 8228.57 | 2.03946 | 0.00225879 | parabolic |
| 4-Kronecker, $(2,5)$ | 190.728 | 3.77976 | 0.396728 | parabolic |
| 3-Kronecker, $(3,4)$ | 78.5455 | 2.89469 | 0.964042 | parabolic |
| 4-Kronecker, $(3,4)$ | 5.10069 | 21.8548 | — | first-call Weyl / parabolic < 0.734 |
| 5-Kronecker, $(2,5)$ | 19.3353 | 10.4463 | >30 (one warm call) | one warm call: Weyl / parabolic < 0.348 |
| $Q(2,1,1), (2,2,3)$ | 38.8315 | 0.0366369 | 0.00663238 | parabolic |
| $Q(0,1,2), (2,3,3)$ | 576 | 0.451287 | 0.00637246 | parabolic |
| $Q(2,1,1), (2,3,3)$ | 136.421 | 2.46452 | 0.0387191 | parabolic |
| $Q(0,1,2), (1,3,4)$ | 1885.09 | 1.26402 | 0.00752854 | parabolic |
| $Q(3,2,1), (1,3,4)$ | 1792 | 20.3574 | 0.00806129 | parabolic |
| $Q(1,2,1), (2,3,4)$ | 2748.14 | — | 0.338731 | first-call Weyl / parabolic > 20.4 |
| $Q(0,1,2), (3,3,4)$ | 116640 | — | 0.0359487 | first-call Weyl / parabolic > 20.2 |
| 5-subspace, final rank 2 | 0.952381 | 0.00136838 | 0.00540588 | Weyl |
| 9-subspace, final rank 2 | 0.0394939 | 0.089164 | 0.990515 | Weyl |
| 4-Kronecker, $(3,5)$ | 47.4349 | — | — | unresolved |
| 6-Kronecker, $(2,5)$ | 3.71793 | — | — | unresolved |
| 5-Kronecker, $(3,4)$ | 0.785981 | — | — | unresolved |
| 8-Kronecker, $(2,3)$ | 0.0279097 | 0.144602 | — | first-call Weyl / parabolic < 0.00513 |
| 3-Kronecker, $(4,5)$ | 2749.45 | — | — | unresolved |

Including three additional fresh-process pairs with alternating algorithm order, parabolic wins 4 of 4 clear warm comparisons on 5-Kronecker $(1,4)$.
The Weyl/parabolic warm ratios range from 1.426 to 1.544.
These are repeat measurements of one model, not additional cutoff models.

## Preselected checks near the cutoff

After seeing the original sample, 6 further three-vertex cases were selected by score and rank-vector diversity, before measuring their timings.
All passed separate checks for nonemptiness, absence of properly semistable representations and ample stability.
The candidate warm threshold was fixed at 24 before running this group; 32 was retained for comparison.
These cases test the proposed revision; they are not included in the original-sample fit.

At threshold 24, clear warm mistakes are: `holdout_230_233`, `holdout_013_124`, `holdout_221_123`, `holdout_121_223`. The worst measured warm selection penalty is 5.198×.

At threshold 32, clear warm mistakes are: `holdout_230_233`, `holdout_013_124`, `holdout_032_124`, `holdout_221_123`, `holdout_121_223`. The worst measured warm selection penalty is 5.198×.


| Model | Score | Weyl warm (s) | Parabolic warm (s) | Comparison |
|:--|--:|--:|--:|:--|
| $Q(2,3,0), (2,3,3)$ | 21.8919 | 0.517019 | 0.0994721 | parabolic |
| $Q(0,1,3), (1,2,4)$ | 21.8907 | 0.143206 | 0.0414176 | parabolic |
| $Q(0,3,2), (1,2,4)$ | 27.3634 | 0.269748 | 0.08526 | parabolic |
| $Q(2,2,1), (1,2,3)$ | 20 | 0.00512258 | 0.0018035 | parabolic |
| $Q(1,2,1), (2,2,3)$ | 18.432 | 0.0804539 | 0.0325347 | parabolic |
| $Q(2,2,0), (2,2,3)$ | 32.3596 | 0.0179595 | 0.00456592 | parabolic |

## What cutoff to use

The 5-Kronecker $(2,5)$ probe has score 19.3353. Its parabolic first call finished in 97.906 s, then a cache-cleared warm construction exceeded 30 s. Weyl's three-sample warm median is 10.446 s.
The preselected $Q(1,2,1), (2,2,3)$ has the lower score 18.432 and favours parabolic by 2.473× in warm medians.
These observations do not support a sharp universal cutoff in this score.
A practical provisional guide is Weyl below 16 and parabolic above 32, with scores 16–32 treated as a gray zone for direct measurement. Some rows remain censored, so these bands are provisional; first-call choices also depend strongly on compilation.
No automatic selection has been added to the library.

For three-vertex rows, $Q(a,b,c)$ has adjacency matrix with upper-triangular entries $(a,b,c)$.
The original 29-model rerun reuses the fitting sample; the six preselected new cases provide a separate, limited check.

## Reproduce

The [raw results](parabolic-chow-results.json) contain matrices, dimension vectors, stability parameters, individual samples, resource measurements, cutoff features and source hashes.
Instantiate the repository environment with the recorded Julia and Singular versions and set a writable `JULIA_DEPOT_PATH`.
Python 3.11 or newer is required. Run the benchmark separately from tests and other Julia jobs.

```sh
python3 benchmark/run_chow.py --julia /path/to/julia --timeout 180 --memory-limit-gib 20
python3 benchmark/run_chow.py --julia /path/to/julia --timeout 30 --memory-limit-gib 20 \
  --cases K4_13 K5_14 K6_15 K7_16 K8_14 K10_15 \
    K3_23 K4_23 K5_23 K7_23 K3_25 K4_25 \
    K3_34 K4_34 K5_25 triangle_211_223 triangle_012_233 triangle_211_233 \
    triangle_012_134 triangle_321_134 triangle_121_234 triangle_012_334 subspace5 subspace9 \
    K4_35 K6_25 K5_34 K8_23 K3_45
python3 benchmark/run_chow.py --julia /path/to/julia --timeout 30 --cases holdout_230_233 holdout_013_124 holdout_032_124 holdout_221_123 holdout_121_223 holdout_220_223
python3 benchmark/run_chow.py --julia /path/to/julia --timeout 180 --warm-timeout 30 --cases K5_25 --algorithms experimental_parabolic
python3 benchmark/summarize_chow.py
/path/to/julia --project=. -e 'using Pkg; Pkg.test()'
```

The second command also reruns the shared small models with a 30-second limit; this does not affect their completed samples.
Generated TOML samples and logs are kept in `benchmark/raw/`.
Complete quotient dimensions agree between algorithms for every paired construction that finished.
The main test suite also checks full ideal equality on its reference models and known Chern numbers.
For the censored large Weyl cases, no reference ideal equality is claimed.
