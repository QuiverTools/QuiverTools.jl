"""Summarize the measured timings and assess the proposed Chow algorithm cutoff."""

import argparse
from datetime import date
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import subprocess
import tomllib

BASE_LABELS = {
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


def label(case, measurement):
    """Give a short geometric or quiver label without assuming a classification."""
    if case in BASE_LABELS:
        return BASE_LABELS[case]
    ranks = ",".join(map(str, measurement["d"]))
    if case.startswith("K"):
        return f"{measurement['adjacency'][0][1]}-Kronecker, $({ranks})$"
    if case.startswith("subspace"):
        return (
            f"{len(measurement['d']) - 1}-subspace, final rank {measurement['d'][-1]}"
        )
    adjacency = measurement["adjacency"]
    return f"$Q({adjacency[0][1]},{adjacency[0][2]},{adjacency[1][2]}), ({ranks})$"


def weighted_border_count(ranks, dimension):
    """Count minimal monomials above the dimension, before normalization."""
    weights = [weight for rank in ranks for weight in range(1, rank + 1)]
    result = 0
    for excess in range(1, max(weights) + 1):
        counts = [1] + [0] * (dimension + excess)
        for weight in weights:
            if weight >= excess:
                for degree in range(weight, len(counts)):
                    counts[degree] += counts[degree - weight]
        result += counts[dimension + excess]
    return result


def cutoff_features(measurement):
    """Compute W, F, B and the proposed score using the presentation's conventions."""
    ranks, adjacency, theta = (measurement[key] for key in ("d", "adjacency", "theta"))
    n = len(ranks)
    dimension = (
        1
        + sum(adjacency[i][j] * ranks[i] * ranks[j] for i in range(n) for j in range(n))
        - sum(rank**2 for rank in ranks)
    )
    sources = [all(adjacency[j][i] == 0 for j in range(n)) for i in range(n)]
    sinks = [all(adjacency[i][j] == 0 for j in range(n)) for i in range(n)]
    forbidden = [
        e
        for e in itertools.product(*(range(rank + 1) for rank in ranks))
        if sum(theta[i] * e[i] for i in range(n)) > 0
    ]

    def precedes(f, e):
        return all(
            (not sources[i] or f[i] <= e[i])
            and (not sinks[i] or f[i] >= e[i])
            and (sources[i] or sinks[i] or f[i] == e[i])
            for i in range(n)
        )

    forbidden = [
        e for e in forbidden if not any(f != e and precedes(f, e) for f in forbidden)
    ]
    order = math.prod(math.factorial(rank) for rank in ranks)
    border = weighted_border_count(ranks, dimension)
    return dict(
        dimension=dimension,
        W=order,
        F=len(forbidden),
        B=border,
        score=order**2 * len(forbidden) / border,
    )


def warm(record):
    """Use a median only when all three verified warm samples are present."""
    samples = record["measurement"].get("samples", [])
    if len(samples) != 4 or any(
        "chow_total_dimension" not in sample for sample in samples
    ):
        return None
    return statistics.median(sample["seconds"] for sample in samples[1:])


def first_text(record):
    samples = record["measurement"].get("samples", [])
    if samples:
        return f"{samples[0]['seconds']:.6g}"
    if record["status"] == "timeout" and record["phase"] == "trial=0":
        return f">{record['timeout_seconds']:g}"
    return record["status"]


def warm_text(record):
    if "warm_lower_bound_seconds" in record:
        return f">{record['warm_lower_bound_seconds']:g} (one warm call)"
    seconds = warm(record)
    return f"{seconds:.6g}" if seconds is not None else "—"


def cutoff_assessment(records, names):
    """Assess warm comparisons separately from first-call censoring bounds."""
    lookup = {(record["case"], record["algorithm"]): record for record in records}
    rows = []
    for name in names:
        reference = lookup[name, "weyl"]
        parabolic = lookup[name, "experimental_parabolic"]
        features = cutoff_features(reference["measurement"])
        old, new = warm(reference), warm(parabolic)
        row = dict(case=name, **features, weyl_warm=old, parabolic_warm=new)
        reference_samples = reference["measurement"].get("samples", [])
        parabolic_samples = parabolic["measurement"].get("samples", [])
        row["first_ratio"] = (
            reference_samples[0]["seconds"] / parabolic_samples[0]["seconds"]
            if reference_samples and parabolic_samples
            else None
        )
        if old is not None and new is not None:
            row.update(kind="warm", ratio=old / new, tie=1 / 1.1 <= old / new <= 1.1)
        elif "warm_lower_bound_seconds" in parabolic and old is not None:
            row.update(
                kind="warm_upper_bound",
                ratio=old / parabolic["warm_lower_bound_seconds"],
            )
        elif (
            reference["status"] == "timeout"
            and reference["phase"] == "trial=0"
            and parabolic["measurement"].get("samples")
        ):
            row.update(
                kind="first_lower_bound",
                ratio=reference["timeout_seconds"]
                / parabolic["measurement"]["samples"][0]["seconds"],
            )
        elif (
            parabolic["status"] == "timeout"
            and parabolic["phase"] == "trial=0"
            and reference["measurement"].get("samples")
        ):
            row.update(
                kind="first_upper_bound",
                ratio=reference["measurement"]["samples"][0]["seconds"]
                / parabolic["timeout_seconds"],
            )
        else:
            row.update(kind="unresolved", ratio=None)
        rows.append(row)
    measured = [row for row in rows if row["kind"] == "warm"]
    decisive = [row for row in measured if not row["tie"]]
    weyl_scores = [row["score"] for row in decisive if row["ratio"] < 1]
    parabolic_scores = [row["score"] for row in decisive if row["ratio"] > 1]
    lower = max(weyl_scores, default=0)
    upper = min(parabolic_scores, default=math.inf)
    first_decisive = [
        row
        for row in rows
        if row["first_ratio"] is not None and not 1 / 1.1 <= row["first_ratio"] <= 1.1
    ]
    checks = []
    for threshold in (16, 20, 24, 32):
        wrong = [
            row["case"]
            for row in decisive
            if (row["score"] >= threshold) != (row["ratio"] > 1)
        ]
        penalties = [
            max(1, 1 / row["ratio"] if row["score"] >= threshold else row["ratio"])
            for row in measured
        ]
        checks.append(
            dict(
                threshold=threshold,
                decisive_misses=wrong,
                worst_warm_penalty=max(penalties, default=1),
                first_call_misses=[
                    row["case"]
                    for row in first_decisive
                    if (row["score"] >= threshold) != (row["first_ratio"] > 1)
                ],
            )
        )
    return dict(
        rows=rows,
        warm_pairs=len(measured),
        decisive_pairs=len(decisive),
        decisive_first_pairs=len(first_decisive),
        clear_weyl_max_score=lower,
        clear_parabolic_min_score=upper,
        thresholds=checks,
    )


def cutoff_table(assessment, lookup):
    """Render warm winners and first-call censoring bounds in separate terms."""
    lines = [
        "",
        "| Model | Score | Weyl warm (s) | Parabolic warm (s) | Comparison |",
        "|:--|--:|--:|--:|:--|",
    ]
    for row in assessment["rows"]:
        old, new = (
            lookup[row["case"], "weyl"],
            lookup[row["case"], "experimental_parabolic"],
        )
        if row["kind"] == "warm":
            comparison = (
                "within 10%"
                if row["tie"]
                else ("parabolic" if row["ratio"] > 1 else "Weyl")
            )
        elif row["kind"] == "first_lower_bound":
            comparison = f"first-call Weyl / parabolic > {row['ratio']:.3g}"
        elif row["kind"] == "first_upper_bound":
            comparison = f"first-call Weyl / parabolic < {row['ratio']:.3g}"
        elif row["kind"] == "warm_upper_bound":
            comparison = f"one warm call: Weyl / parabolic < {row['ratio']:.3g}"
        else:
            comparison = "unresolved"
        lines.append(
            f"| {label(row['case'], old['measurement'])} | {row['score']:.6g} | {warm_text(old)} | {warm_text(new)} | {comparison} |"
        )
    return lines


def render_report(result, previous):
    """Produce the timing table, cleanup comparison and cutoff audit together."""
    records = result["records"]
    lookup = {(row["case"], row["algorithm"]): row for row in records}
    environment = records[0]["measurement"]
    lines = [
        "# Experimental parabolic Chow benchmark",
        "",
        "**LLM-suggested experiment; the proof and code need careful checking.**",
        "The reduction and hypotheses are in [the documentation](../docs/src/experimental-chow.md).",
        "",
        f"Measured {result['date']} on {environment['cpu']}, Julia {environment['julia']}, Singular.jl {environment['singular']}, one Julia thread.",
        f"The Chow implementation is commit `{result['source_commit']}`; the worker adds the cutoff models listed in the raw data.",
        "Both algorithms ran serially in fresh processes, with the same generic warmup on 3-Kronecker $(1,1)$.",
        "The first timed call includes remaining case-specific compilation.",
        "Each warm time is the median of three further constructions, with Memoization results cleared and full GC before every sample.",
        "Compiled code and Singular's polynomial parent-ring caches persist; relations, Gröbner bases and experimental Betti numbers are recomputed.",
        "Startup, warmup, cache clearing, preceding GC and subsequent quotient-dimension checks are outside the timer.",
        "GC during construction is included. Both algorithms use `unsafe=true`; parabolic timings include the Betti check.",
        "The two Chern rows include a fresh moduli object, ring construction and all Chern numbers.",
        "",
        "Each phase has a three-minute maximum and a 20 GiB process-tree memory limit.",
        "The additional cutoff runs have a 30-second phase limit; shared rows reuse the longer main-table runs.",
        "For 5-Kronecker $(2,5)$, an additional probe allowed 180 seconds for the first call and 30 seconds for the next warm call. Its warm bound refers to one construction, not a three-sample median.",
        "On macOS the watchdog measures both resident memory and physical footprint and kills the process group if either exceeds the limit.",
        "",
        "## Main timings",
        "",
        "Times are in seconds. A bound in the first-call column is a timeout, not an exact time.",
        "",
        "| Model / computation | Weyl first | Parabolic first | Weyl warm | Parabolic warm | Weyl / parabolic peak RSS (GiB) |",
        "|:--|--:|--:|--:|--:|--:|",
    ]
    for name, title in BASE_LABELS.items():
        old, new = lookup[name, "weyl"], lookup[name, "experimental_parabolic"]
        lines.append(
            f"| {title} | {first_text(old)} | {first_text(new)} | {warm_text(old)} | {warm_text(new)} | {old['peak_rss_bytes'] / 2**30:.3f} / {new['peak_rss_bytes'] / 2**30:.3f} |"
        )
    lines += [
        "",
        "## Effect of the cleanup",
        "",
        "The earlier measurements used the same machine, dependency versions and three-warm-sample convention.",
        "This comparison measures the combined cleanup; it does not isolate the effect of individual edits.",
        "The benchmark worker was also moved into a function during cleanup, which can affect overhead on the smallest rows.",
        "",
        "| Model | Previous parabolic warm (s) | Current parabolic warm (s) | Improvement | Current cumulative Julia allocation (GB) |",
        "|:--|--:|--:|--:|--:|",
    ]
    for name in ("chain_552", "fork_552", "K3_34", "K3_35", "K3_23", "subspace7"):
        old = next(
            row
            for row in previous["records"]
            if row["case"] == name and row["algorithm"] == "experimental_parabolic"
        )
        old_warm = statistics.median(
            sample["seconds"] for sample in old["measurement"]["samples"][1:]
        )
        new = lookup[name, "experimental_parabolic"]
        new_warm = warm(new)
        allocation = (
            statistics.median(
                sample["julia_bytes"] for sample in new["measurement"]["samples"][1:]
            )
            / 1e9
        )
        lines.append(
            f"| {BASE_LABELS[name]} | {old_warm:.6g} | {new_warm:.6g} | {old_warm / new_warm:.3f}× | {allocation:.3f} |"
        )
    for name in ("chain_552", "fork_552", "K3_35"):
        old, new = lookup[name, "weyl"], lookup[name, "experimental_parabolic"]
        if old["status"] == "timeout" and old["phase"] == "trial=0":
            first = new["measurement"]["samples"][0]["seconds"]
            bound = math.floor(100 * old["timeout_seconds"] / first) / 100
            lines += [
                "",
                f"For {BASE_LABELS[name]}, Weyl exceeded {old['timeout_seconds']:g} seconds on its first call; the first-call speedup is therefore greater than {bound:g}×.",
            ]
    peak_rss = max(row["peak_rss_bytes"] for row in records) / 2**30
    peak_footprint = max(row["peak_footprint_bytes"] for row in records) / 2**30
    lines += [
        "",
        f"Across all runs, maximum resident memory was {peak_rss:.3f} GiB and maximum physical footprint was {peak_footprint:.3f} GiB.",
        "Julia allocation counts are cumulative and exclude native Singular allocations; they are not peak memory.",
        "",
        "## Cutoff assessment",
        "",
        "Let $W=\\prod_i d_i!$, let $F$ count the minimal forbidden vectors, and let $B$ count the minimal weighted monomials above the dimension in the ambient Chern ring, before imposing the normalization relation.",
        "The proposed score is $S=FW^2/B$, with parabolic chosen when $S\\geq32$.",
        "The $W^2$ term reflects the reference basis enumeration and explicit alternating sums; $B$ approximates the truncation overhead in the parabolic method.",
        "This is an empirical proxy for the current implementations, not a complexity theorem.",
        "",
    ]
    assessment = result["cutoff"]
    lines += [
        f"The original 29-model cutoff sample was rerun. There are {assessment['warm_pairs']} complete warm comparisons and {assessment['decisive_pairs']} outside a 10% tie window.",
        "Only complete, three-sample warm medians enter the threshold-fit table.",
        "First-call timeouts provide first-call bounds; the targeted warm timeout is shown separately. Rows with both algorithms censored remain unresolved.",
        "",
        f"The largest score among clear warm Weyl wins is {assessment['clear_weyl_max_score']:.6g}; the smallest among clear warm parabolic wins is {assessment['clear_parabolic_min_score']:.6g}.",
        "",
        "The first-call column below is assessed separately, using completed first calls outside the same 10% tie window.",
        "",
        "| Threshold | Clear warm mistakes | Worst warm slowdown | Clear first-call mistakes |",
        "|--:|:--|--:|:--|",
    ]
    for row in assessment["thresholds"]:
        names = ", ".join(f"`{name}`" for name in row["decisive_misses"]) or "none"
        first_names = (
            ", ".join(f"`{name}`" for name in row["first_call_misses"]) or "none"
        )
        lines.append(
            f"| {row['threshold']} | {names} | {row['worst_warm_penalty']:.3f}× | {first_names} |"
        )
    boundary_reference = lookup["K5_14", "weyl"]
    boundary_parabolic = lookup["K5_14", "experimental_parabolic"]
    reference_first = boundary_reference["measurement"]["samples"][0]["seconds"]
    parabolic_first = boundary_parabolic["measurement"]["samples"][0]["seconds"]
    lines += [
        "",
        f"At score 25.0435, 5-Kronecker $(1,4)$ takes {1000 * warm(boundary_reference):.3f} ms warm with Weyl and {1000 * warm(boundary_parabolic):.3f} ms with parabolic.",
        f"Its first timed calls go the other way: {reference_first:.4f} s with Weyl and {parabolic_first:.4f} s with parabolic.",
        "A threshold for warmed-up calculations need not give the best first-call choice under this particular generic warmup.",
    ]
    lines += cutoff_table(assessment, lookup)
    if result.get("boundary_repeats"):
        pairs = [(warm(boundary_reference), warm(boundary_parabolic))]
        for repeat in sorted({row["repeat"] for row in result["boundary_repeats"]}):
            values = {
                row["algorithm"]: warm(row)
                for row in result["boundary_repeats"]
                if row["repeat"] == repeat
            }
            pairs.append((values["weyl"], values["experimental_parabolic"]))
        ratios = [reference / parabolic for reference, parabolic in pairs]
        lines += [
            "",
            f"Including three additional fresh-process pairs with alternating algorithm order, parabolic wins {sum(ratio > 1.1 for ratio in ratios)} of {len(ratios)} clear warm comparisons on 5-Kronecker $(1,4)$.",
            f"The Weyl/parabolic warm ratios range from {min(ratios):.3f} to {max(ratios):.3f}.",
            "These are repeat measurements of one model, not additional cutoff models.",
        ]
    if result.get("holdout_cases"):
        holdout = result["holdout"]
        lines += [
            "",
            "## Preselected checks near the cutoff",
            "",
            f"After seeing the original sample, {len(result['holdout_cases'])} further three-vertex cases were selected by score and rank-vector diversity, before measuring their timings.",
            "All passed separate checks for nonemptiness, absence of properly semistable representations and ample stability.",
            "The candidate warm threshold was fixed at 24 before running this group; 32 was retained for comparison.",
            "These cases test the proposed revision; they are not included in the original-sample fit.",
            "",
        ]
        for check in holdout["thresholds"]:
            if check["threshold"] in (24, 32):
                misses = (
                    ", ".join(f"`{name}`" for name in check["decisive_misses"])
                    or "none"
                )
                lines += [
                    f"At threshold {check['threshold']}, clear warm mistakes are: {misses}. The worst measured warm selection penalty is {check['worst_warm_penalty']:.3f}×.",
                    "",
                ]
        lines += cutoff_table(holdout, lookup)
    diagnostic = lookup["K5_25", "experimental_parabolic"]
    if "warm_lower_bound_seconds" in diagnostic:
        reference = lookup["K5_25", "weyl"]
        score = cutoff_features(reference["measurement"])["score"]
        first = diagnostic["measurement"]["samples"][0]["seconds"]
        bound = diagnostic["warm_lower_bound_seconds"]
        lines += [
            "",
            "## What cutoff to use",
            "",
            f"The 5-Kronecker $(2,5)$ probe has score {score:.6g}. Its parabolic first call finished in {first:.3f} s, then a cache-cleared warm construction exceeded {bound:g} s. Weyl's three-sample warm median is {warm(reference):.3f} s.",
        ]
        lower_wins = [
            row
            for row in (result.get("holdout") or {}).get("rows", [])
            if row["kind"] == "warm" and row["ratio"] > 1.1 and row["score"] < score
        ]
        if lower_wins:
            witness = min(lower_wins, key=lambda row: row["score"])
            title = label(
                witness["case"], lookup[witness["case"], "weyl"]["measurement"]
            )
            lines += [
                f"The preselected {title} has the lower score {witness['score']:.6g} and favours parabolic by {witness['ratio']:.3f}× in warm medians.",
                "These observations do not support a sharp universal cutoff in this score.",
            ]
            measured = [
                row
                for row in assessment["rows"] + result["holdout"]["rows"]
                if row["kind"] == "warm" and not row["tie"]
            ]
            outside_agrees = all(
                (row["score"] >= 16 or row["ratio"] < 1)
                and (row["score"] <= 32 or row["ratio"] > 1)
                for row in measured
            )
            if outside_agrees:
                lines += [
                    "A practical provisional guide is Weyl below 16 and parabolic above 32, with scores 16–32 treated as a gray zone for direct measurement. Some rows remain censored, so these bands are provisional; first-call choices also depend strongly on compilation."
                ]
        lines += ["No automatic selection has been added to the library."]
    lines += [
        "",
        "For three-vertex rows, $Q(a,b,c)$ has adjacency matrix with upper-triangular entries $(a,b,c)$.",
        "The original 29-model rerun reuses the fitting sample; the six preselected new cases provide a separate, limited check.",
        "",
        "## Reproduce",
        "",
        "The [raw results](parabolic-chow-results.json) contain matrices, dimension vectors, stability parameters, individual samples, resource measurements, cutoff features and source hashes.",
        "Instantiate the repository environment with the recorded Julia and Singular versions and set a writable `JULIA_DEPOT_PATH`.",
        "Python 3.11 or newer is required. Run the benchmark separately from tests and other Julia jobs.",
        "",
        "```sh",
        "python3 benchmark/run_chow.py --julia /path/to/julia --timeout 180 --memory-limit-gib 20",
        "python3 benchmark/run_chow.py --julia /path/to/julia --timeout 30 --memory-limit-gib 20 \\",
        "  --cases "
        + " \\\n    ".join(
            " ".join(result["cutoff_cases"][start : start + 6])
            for start in range(0, len(result["cutoff_cases"]), 6)
        ),
        "python3 benchmark/run_chow.py --julia /path/to/julia --timeout 30 --cases "
        + " ".join(result.get("holdout_cases", [])),
        "python3 benchmark/run_chow.py --julia /path/to/julia --timeout 180 --warm-timeout 30 --cases K5_25 --algorithms experimental_parabolic",
        "python3 benchmark/summarize_chow.py",
        "/path/to/julia --project=. -e 'using Pkg; Pkg.test()'",
        "```",
        "",
        "The second command also reruns the shared small models with a 30-second limit; this does not affect their completed samples.",
        "Generated TOML samples and logs are kept in `benchmark/raw/`.",
        "Complete quotient dimensions agree between algorithms for every paired construction that finished.",
        "The main test suite also checks full ideal equality on its reference models and known Chern numbers.",
        "For the censored large Weyl cases, no reference ideal equality is claimed.",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--cutoff-catalog", type=Path)
    parser.add_argument("--source-commit")
    parser.add_argument("--holdout-catalog", type=Path)
    parser.add_argument("--boundary-results", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    raw = args.raw or root / "benchmark/raw"
    output = root / "benchmark/parabolic-chow-results.json"
    existing = json.loads(output.read_text())
    previous = (
        json.loads(args.previous.read_text())
        if args.previous
        else existing.get("previous", existing)
    )
    names = (
        [row["name"] for row in json.loads(args.cutoff_catalog.read_text())]
        if args.cutoff_catalog
        else existing["cutoff_cases"]
    )
    holdout_names = (
        [row["name"] for row in json.loads(args.holdout_catalog.read_text())]
        if args.holdout_catalog
        else existing.get("holdout_cases", [])
    )
    boundary_repeats = (
        json.loads(args.boundary_results.read_text()) if args.boundary_results else []
    )
    records = json.loads((raw / "status.json").read_text())
    for record in records:
        sample_path = raw / f"{record['case']}_{record['algorithm']}.toml"
        record["measurement"] = tomllib.loads(sample_path.read_text())
    lookup = {(row["case"], row["algorithm"]): row for row in records}
    for name in {row["case"] for row in records}:
        values = [
            sample["chow_total_dimension"]
            for algorithm in ("weyl", "experimental_parabolic")
            for sample in lookup[name, algorithm]["measurement"].get("samples", [])
            if "chow_total_dimension" in sample
        ]
        if values:
            assert len(set(values)) == 1, name
    result = dict(
        date=str(date.today()),
        source_commit=args.source_commit
        or subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        reference_commit=subprocess.check_output(
            ["git", "rev-parse", "main"], cwd=root, text=True
        ).strip(),
        source_sha256={
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                root / "src/Chow.jl",
                root / "src/ExperimentalChow.jl",
                root / "benchmark/chow_worker.jl",
            )
        },
        resolved_dependencies={
            name: [
                {
                    key: entry[key]
                    for key in ("uuid", "version", "git-tree-sha1")
                    if key in entry
                }
                for entry in entries
            ]
            for name, entries in tomllib.loads((root / "Manifest.toml").read_text())[
                "deps"
            ].items()
        },
        records=records,
        cutoff_cases=names,
        cutoff=cutoff_assessment(records, names),
        holdout_cases=holdout_names,
        holdout=cutoff_assessment(records, holdout_names) if holdout_names else None,
        boundary_repeats=boundary_repeats,
        previous=previous,
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    (root / "benchmark/parabolic-chow.md").write_text(render_report(result, previous))
    print(json.dumps(result["cutoff"]["thresholds"], indent=2))


if __name__ == "__main__":
    main()
