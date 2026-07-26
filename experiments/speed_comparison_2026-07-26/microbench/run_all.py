#!/usr/bin/env python3
"""Run the persistent-process microbenchmark for every verification matrix."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import time


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent


def run_one(binary: Path, matrix: str, target_seconds: float, cpu: int, env: dict[str, str]) -> dict:
    completed = subprocess.run(
        ["taskset", "-c", str(cpu), str(binary), matrix, str(target_seconds)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=True,
    )
    values: dict[str, int | float] = {}
    integer_keys = {
        "repetitions",
        "batches",
        "warmup_repetitions",
        "calibration_repetitions",
        "ess_count",
        "candidate_count",
        "checksum",
    }
    for line in completed.stdout.splitlines():
        key, raw = line.split("=", 1)
        values[key] = int(raw) if key in integer_keys else float(raw)
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-seconds", type=float, default=3.0)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument(
        "--output",
        type=Path,
        default=EXPERIMENT / "results" / "microbench_exact_pd.json",
    )
    args = parser.parse_args()

    matrices_path = EXPERIMENT / "input" / "verification_matrices.json"
    matrices = sorted(
        json.loads(matrices_path.read_text())["matrices"],
        key=lambda row: int(row["id"]),
    )
    binaries = {
        "current": (HERE / "build" / "benchmark_current").resolve(),
        "exact_pd_only": (HERE / "build" / "benchmark_exact_pd_only").resolve(),
    }

    env = os.environ.copy()
    libdir = "/tmp/fracessa-sysroot/usr/lib64"
    env["LD_LIBRARY_PATH"] = libdir + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")

    started = time.time()
    results = []
    for ordinal, matrix in enumerate(matrices, start=1):
        matrix_id = int(matrix["id"])
        expected_ess = int(matrix["number_ess"])
        matrix_arg = f'{int(matrix["dimension"])}#{matrix["matrix"]}'
        order = ("current", "exact_pd_only") if matrix_id & 1 else ("exact_pd_only", "current")
        measurements = {}
        for name in order:
            measurements[name] = run_one(
                binaries[name], matrix_arg, args.target_seconds, args.cpu, env
            )
            if measurements[name]["ess_count"] != expected_ess:
                raise RuntimeError(
                    f"matrix {matrix_id}: {name} returned {measurements[name]['ess_count']} ESS, "
                    f"expected {expected_ess}"
                )

        if measurements["current"]["candidate_count"] != measurements["exact_pd_only"]["candidate_count"]:
            raise RuntimeError(f"matrix {matrix_id}: candidate counts differ")

        current_ns = measurements["current"]["median_ns"]
        exact_ns = measurements["exact_pd_only"]["median_ns"]
        results.append(
            {
                "id": matrix_id,
                "dimension": int(matrix["dimension"]),
                "number_ess_expected": expected_ess,
                "current": measurements["current"],
                "exact_pd_only": measurements["exact_pd_only"],
                "ratio": exact_ns / current_ns,
            }
        )
        print(
            f"[{ordinal:02d}/{len(matrices)}] id={matrix_id:02d} "
            f"current={current_ns / 1000.0:.3f}us "
            f"exact={exact_ns / 1000.0:.3f}us "
            f"change={(exact_ns / current_ns - 1.0) * 100.0:+.2f}%",
            flush=True,
        )

    payload = {
        "metadata": {
            "target_seconds": args.target_seconds,
            "cpu_affinity": args.cpu,
            "hostname": platform.node(),
            "machine": platform.machine(),
            "processing_time": time.time() - started,
            "protocol": (
                "matrix parsed once; candidates enabled; analyzer constructed and destroyed per repetition; "
                "50ms warmup for fast matrices; up to 7 adaptive batches; steady_clock nanoseconds"
            ),
            "binaries": {name: str(path) for name, path in binaries.items()},
        },
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
