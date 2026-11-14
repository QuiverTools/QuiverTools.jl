"""Run serial Chow benchmarks with per-phase time and process-tree memory limits.

Usage: python3 benchmark/run_chow.py --julia /path/to/julia
Memory monitoring supports macOS and Linux and uses only the standard library.
"""

import argparse
import ctypes
import errno
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


class ProcessMemory:
    """Read resident memory and physical footprint for a process and its children."""

    def __init__(self):
        if sys.platform == "darwin":
            self.library = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
            self.library.proc_pidinfo.argtypes = [
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_uint64,
                ctypes.c_void_p,
                ctypes.c_int,
            ]
            self.library.proc_pid_rusage.argtypes = [
                ctypes.c_int,
                ctypes.c_int,
                ctypes.c_void_p,
            ]
            self.library.proc_listchildpids.argtypes = [
                ctypes.c_int,
                ctypes.c_void_p,
                ctypes.c_int,
            ]
        elif not sys.platform.startswith("linux"):
            raise RuntimeError("memory monitoring requires macOS or Linux")

    def children(self, pid):
        """Return child PIDs; macOS returns their count, not the buffer's byte size."""
        if sys.platform != "darwin":
            path = Path(f"/proc/{pid}/task/{pid}/children")
            return list(map(int, path.read_text().split()))
        buffer = (ctypes.c_int * 4096)()
        count = self.library.proc_listchildpids(pid, buffer, ctypes.sizeof(buffer))
        if count < 0:
            raise OSError(ctypes.get_errno(), "cannot list child processes")
        if count >= len(buffer):
            raise RuntimeError("child-process buffer exhausted")
        return [child for child in buffer[:count] if child > 0]

    def memory(self, pid):
        """Return `(resident_bytes, footprint_bytes)` for one process."""
        if sys.platform != "darwin":
            for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    resident = int(line.split()[1]) * 1024
                    return resident, resident
            return 0, 0
        # PROC_PIDTASKINFO consists of twelve uint64_t fields; resident size is
        # the second. proc_pid_rusage v0 has sixteen UUID bytes followed by ten
        # uint64_t fields, with physical footprint at index seven.
        task = (ctypes.c_uint64 * 12)()
        usage = (ctypes.c_uint64 * 12)()
        if self.library.proc_pidinfo(
            pid, 4, 0, task, ctypes.sizeof(task)
        ) != ctypes.sizeof(task):
            raise OSError(ctypes.get_errno(), f"cannot read memory for PID {pid}")
        if self.library.proc_pid_rusage(pid, 0, usage) != 0:
            raise OSError(ctypes.get_errno(), f"cannot read footprint for PID {pid}")
        return task[1], usage[9]

    def sample(self, pid):
        """Sum memory over the process tree, tolerating children that have exited."""
        pids = [pid]
        resident = footprint = 0
        for current in pids:
            try:
                resident_bytes, footprint_bytes = self.memory(current)
                resident += resident_bytes
                footprint += footprint_bytes
                for child in self.children(current):
                    if child not in pids:
                        pids.append(child)
            except OSError as error:
                if error.errno not in (errno.ENOENT, errno.ESRCH):
                    raise
        return resident, footprint


def is_warm_trial(phase):
    """Recognize positive trial markers, tolerating an incomplete phase-file write."""
    name, _, value = phase.partition("=")
    return name == "trial" and value.isdecimal() and int(value) > 0


def run_case(
    command,
    log_path,
    phase_path,
    timeout,
    startup_timeout,
    memory_limit_gib,
    cwd,
    warm_timeout=None,
):
    """Run and reap a worker, preserving its final phase and resource measurements."""
    monitor = ProcessMemory()
    with log_path.open("w") as log:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env={**os.environ, "JULIA_NUM_THREADS": "1"},
        )
        started = changed = time.monotonic()
        last_phase, status = "startup", "complete"
        peak_resident = peak_footprint = 0
        try:
            while process.poll() is None:
                current = (
                    phase_path.read_text().strip() if phase_path.exists() else "startup"
                )
                current = current or last_phase
                if current != last_phase:
                    last_phase, changed = current, time.monotonic()
                try:
                    resident, footprint = monitor.sample(process.pid)
                except OSError:
                    if process.poll() is not None:
                        break
                    raise
                peak_resident = max(peak_resident, resident)
                peak_footprint = max(peak_footprint, footprint)
                limit = startup_timeout if current == "startup" else timeout
                warm_trial = is_warm_trial(current)
                if warm_trial and warm_timeout is not None:
                    limit = warm_timeout
                if max(resident, footprint) > memory_limit_gib * 2**30:
                    status = "memory_limit"
                    break
                if time.monotonic() - changed > limit:
                    status = "timeout"
                    break
                time.sleep(0.05)
        finally:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            process.wait()
        if status == "complete" and process.returncode != 0:
            status = "error"
        result = dict(
            status=status,
            phase=last_phase,
            timeout_seconds=limit if status == "timeout" else timeout,
            warm_timeout_seconds=warm_timeout,
            memory_limit_gib=memory_limit_gib,
            elapsed_wall_seconds=time.monotonic() - started,
            peak_rss_bytes=peak_resident,
            peak_footprint_bytes=peak_footprint,
        )
        if status == "timeout" and is_warm_trial(last_phase):
            result["warm_lower_bound_seconds"] = limit
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument(
        "--warm-timeout",
        type=float,
        help="optional shorter limit for warm constructions",
    )
    parser.add_argument("--startup-timeout", type=float, default=180)
    parser.add_argument("--memory-limit-gib", type=float, default=20)
    parser.add_argument(
        "--algorithms", nargs="+", default=["weyl", "experimental_parabolic"]
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=[
            "P3",
            "Gr24",
            "Gr36",
            "K3_23",
            "subspace7",
            "chain_552",
            "fork_552",
            "K3_34",
            "K3_35",
            "K3_23_chern",
            "subspace7_chern",
        ],
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / "benchmark/raw"
    out.mkdir(parents=True, exist_ok=True)
    status_file = out / "status.json"
    records = json.loads(status_file.read_text()) if status_file.exists() else []
    for case in args.cases:
        for algorithm in args.algorithms:
            stem = f"{case}_{algorithm}"
            destination = out / (stem + ".toml")
            phase = Path(str(destination) + ".phase")
            for path in (destination, phase):
                path.unlink(missing_ok=True)
            command = [
                args.julia,
                "--startup-file=no",
                "--project=.",
                "benchmark/chow_worker.jl",
                case,
                algorithm,
                str(destination),
            ]
            record = dict(
                case=case,
                algorithm=algorithm,
                **run_case(
                    command,
                    out / (stem + ".log"),
                    phase,
                    args.timeout,
                    args.startup_timeout,
                    args.memory_limit_gib,
                    root,
                    warm_timeout=args.warm_timeout,
                ),
            )
            records = [
                row
                for row in records
                if (row["case"], row["algorithm"]) != (case, algorithm)
            ]
            records.append(record)
            temporary = status_file.with_suffix(".tmp")
            temporary.write_text(json.dumps(records, indent=2) + "\n")
            temporary.replace(status_file)
            print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
