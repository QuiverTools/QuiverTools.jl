"""Serial, per-phase time-limited benchmark; no Python dependencies.

Usage: python3 benchmark/run_chow.py --julia /path/to/julia
Use a writable JULIA_DEPOT_PATH with an instantiated repository environment.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument("--julia", default="julia")
parser.add_argument("--timeout", type=float, default=120)
parser.add_argument("--startup-timeout", type=float, default=180)
parser.add_argument("--algorithms", nargs="+", default=["weyl", "experimental_parabolic"])
parser.add_argument("--cases", nargs="+", default=[
    "P3", "Gr24", "Gr36", "K3_23", "subspace7", "chain_552", "fork_552",
    "K3_34", "K3_35", "K3_23_chern", "subspace7_chern"])
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / "benchmark/raw"
out.mkdir(parents=True, exist_ok=True)
status_file = out / "status.json"
records = json.loads(status_file.read_text()) if status_file.exists() else []
for case in args.cases:
    for algorithm in args.algorithms:
        dest = out / f"{case}_{algorithm}.toml"
        phase = Path(str(dest) + ".phase")
        for path in (dest, phase):
            path.unlink(missing_ok=True)
        with open(out / f"{case}_{algorithm}.log", "w") as log:
            process = subprocess.Popen([
                args.julia, "--startup-file=no", "--project=.",
                "benchmark/chow_worker.jl", case, algorithm, str(dest)],
                cwd=root, stdout=log, stderr=subprocess.STDOUT,
                env={**os.environ, "JULIA_NUM_THREADS": "1"})
            started = changed = time.monotonic()
            last_phase = "startup"
            status = "complete"
            while process.poll() is None:
                current = phase.read_text() if phase.exists() else "startup"
                if current != last_phase:
                    last_phase, changed = current, time.monotonic()
                limit = args.startup_timeout if current == "startup" else args.timeout
                if time.monotonic() - changed > limit:
                    process.kill()
                    process.wait()
                    status = "timeout"
                    break
                time.sleep(0.1)
            if status == "complete" and process.returncode != 0:
                status = "error"
            if phase.exists():
                last_phase = phase.read_text()
            record = dict(case=case, algorithm=algorithm, status=status,
                          phase=last_phase, timeout_seconds=args.timeout,
                          elapsed_wall_seconds=time.monotonic() - started)
            records = [r for r in records if
                       (r["case"], r["algorithm"]) != (case, algorithm)]
            records.append(record)
            status_file.write_text(json.dumps(records, indent=2) + "\n")
            print(json.dumps(record), flush=True)
