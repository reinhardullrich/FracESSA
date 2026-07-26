#!/usr/bin/env python3
"""Benchmark one FracESSA binary with the historical verification protocol."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


def iteration_count(dimension: int) -> int:
    if dimension < 10:
        return 20
    if dimension < 15:
        return 5
    return 1


def run_once(binary: Path, matrix_arg: str, cpu: int, env: dict[str, str]) -> tuple[int, float, int, int]:
    command = ["taskset", "-c", str(cpu), str(binary), "-c", "-t", matrix_arg]
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr.decode(errors='replace')}"
        )

    lines = completed.stdout.splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"missing ESS/timing output from {' '.join(command)}")

    return int(lines[0]), float(lines[1]), max(0, len(lines) - 3), len(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--matrices", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=2)
    args = parser.parse_args()

    binary = args.binary.resolve()
    matrices_path = args.matrices.resolve()
    data = json.loads(matrices_path.read_text())
    matrices = sorted(data["matrices"], key=lambda row: int(row["id"]))

    env = os.environ.copy()
    libdir = "/tmp/fracessa-sysroot/usr/lib64"
    env["LD_LIBRARY_PATH"] = libdir + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")

    started = time.time()
    results = []
    for ordinal, matrix in enumerate(matrices, start=1):
        matrix_id = int(matrix["id"])
        dimension = int(matrix["dimension"])
        expected = int(matrix["number_ess"])
        matrix_arg = f'{dimension}#{matrix["matrix"]}'
        samples = []
        candidate_count = None
        output_bytes = None

        for _ in range(iteration_count(dimension)):
            ess_count, elapsed, candidates, byte_count = run_once(binary, matrix_arg, args.cpu, env)
            if ess_count != expected:
                raise RuntimeError(
                    f"matrix {matrix_id}: expected {expected} ESS, binary returned {ess_count}"
                )
            samples.append(elapsed)
            candidate_count = candidates
            output_bytes = byte_count

        row = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": expected,
            "ess_count": expected,
            "iterations": len(samples),
            "timing": min(samples),
            "median_timing": statistics.median(samples),
            "samples": samples,
            "candidate_rows": candidate_count,
            "output_bytes": output_bytes,
        }
        results.append(row)
        print(
            f"[{ordinal:02d}/{len(matrices)}] id={matrix_id:02d} n={dimension:02d} "
            f"min={row['timing']:.6f}s median={row['median_timing']:.6f}s",
            flush=True,
        )

    payload = {
        "metadata": {
            "label": args.label,
            "binary": str(binary),
            "binary_sha256": subprocess.check_output(["sha256sum", str(binary)], text=True).split()[0],
            "matrices": str(matrices_path),
            "timestamp": time.time(),
            "processing_time": time.time() - started,
            "cpu_affinity": args.cpu,
            "hostname": platform.node(),
            "machine": platform.machine(),
            "protocol": "include candidates; safe parser; internal C++ timing; min of 20/5/1 runs for n<10/n<15/n>=15",
        },
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

