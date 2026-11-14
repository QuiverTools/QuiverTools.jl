"""Collect run_chow.py output into the checked-in experimental report."""
from datetime import date
import json
import math
from pathlib import Path
import statistics
import subprocess
import tomllib

root = Path(__file__).resolve().parents[1]
raw = root / "benchmark/raw"
records = json.loads((raw / "status.json").read_text())
for record in records:
    stem = f"{record['case']}_{record['algorithm']}.toml"
    record["measurement"] = tomllib.loads((raw / stem).read_text())
    record["phase"] = Path(str(raw / stem) + ".phase").read_text()
result = dict(date=str(date.today()),
              reference_commit=subprocess.check_output(
                  ["git", "rev-parse", "main"], cwd=root, text=True).strip(),
              resolved_dependencies={name: [{key: entry[key] for key in
                  ("uuid", "version", "git-tree-sha1") if key in entry}
                  for entry in entries] for name, entries in
                  tomllib.loads((root / "Manifest.toml").read_text())["deps"].items()},
              records=records)
(root / "benchmark/parabolic-chow-results.json").write_text(
    json.dumps(result, indent=2) + "\n")

lookup = {(r["case"], r["algorithm"]): r for r in records}
labels = {
    "P3": "$\\mathbb P^3$",
    "Gr24": "$\\mathrm{Gr}(2,4)$",
    "Gr36": "$\\mathrm{Gr}(3,6)$",
    "K3_23": "3-Kronecker, $(2,3)$",
    "subspace7": "Seven-subspace fourfold",
    "chain_552": "Chain, $(5,5,2)$",
    "fork_552": "Fork, $(5,5,2)$",
    "K3_34": "3-Kronecker, $(3,4)$",
    "K3_35": "3-Kronecker, $(3,5)$",
    "K3_23_chern": "3-Kronecker $(2,3)$, all Chern numbers",
    "subspace7_chern": "Seven-subspace fourfold, all Chern numbers",
}
def first(record):
    samples = record["measurement"]["samples"]
    if samples:
        return f"{samples[0]['seconds']:.6g}"
    assert record["status"] == "timeout" and record["phase"] == "trial=0"
    return f">{record['timeout_seconds']:g} (timeout)"

def warm(record):
    samples = record["measurement"]["samples"][1:]
    return f"{1000 * statistics.median(s['seconds'] for s in samples):.6g}" if samples else "—"

rows = []
for case, label in labels.items():
    if (case, "weyl") not in lookup or (case, "experimental_parabolic") not in lookup:
        continue
    old, new = lookup[case, "weyl"], lookup[case, "experimental_parabolic"]
    old_dims = [s["chow_total_dimension"] for s in old["measurement"]["samples"]]
    new_dims = [s["chow_total_dimension"] for s in new["measurement"]["samples"]]
    assert len(set(old_dims + new_dims)) == 1
    rows.append(f"| {label} | {first(old)} | {first(new)} | {warm(old)} | {warm(new)} |")

env = records[0]["measurement"]
text = f"""# Experimental parabolic Chow benchmark

**LLM-suggested experiment; careful independent checking is required.**
These measurements do not establish correctness or a uniform performance improvement.
The proof and geometric hypotheses are in
[the experimental documentation](../docs/src/experimental-chow.md).

Measured {result['date']} on {env['cpu']}, arm64, Julia {env['julia']},
Singular.jl {env['singular']}, one Julia thread.
The reference is the unchanged `:weyl` backend from
`{result['reference_commit']}`.
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
""" + "\n".join(rows) + "\n\n"
for case in ("chain_552", "fork_552", "K3_35"):
    old, new = lookup[case, "weyl"], lookup[case, "experimental_parabolic"]
    if old["status"] == "timeout" and old["phase"] == "trial=0":
        seconds = new["measurement"]["samples"][0]["seconds"]
        allocated = statistics.median(s["julia_bytes"] for s in new["measurement"]["samples"][1:]) / 1e9
        text += (f"For **{labels[case]}**, the reference was still inside the first timed "
                 f"construction at the {old['timeout_seconds']:g}-second deadline. "
                 f"The first-call speedup is therefore **greater than {math.floor(100 * old['timeout_seconds'] / seconds) / 100:g}×**; "
                 f"the exact speedup is unknown. The experimental warm constructions allocated "
                 f"about {allocated:.2f} GB cumulatively in Julia per call.\n\n")
text += """These are cumulative allocation counts, not peak resident memory.
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
"""
(root / "benchmark/parabolic-chow.md").write_text(text)
