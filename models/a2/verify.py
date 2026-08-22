#!/usr/bin/env python3
"""Compare safe-only model a2 with the unchanged production baseline."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sqlite3
import statistics
import subprocess
import sys


MODEL_DIR = Path(__file__).resolve().parent
ROOT = MODEL_DIR.parents[1]
BUILD = MODEL_DIR / "build"
DATABASE = ROOT / "testdata" / "fracessa_testdata.sqlite3"
DATABASE_IDS = (9, 11, 22, 34, 35, 46, 51, 71, 2010, 2070, 2076, 2153, 2344, 2440, 2689, 2820, 2975)
DIRECT_CASES = (
    ("pair_ceiling", "3#0,0,0,0,0,0", None, None),
    ("larger_ceiling", "3#-1,2,0,-1,0,0", None, None),
    ("singular_inertia_bound", "3#-2,0,-2,-2,-2,-3", None, None),
    ("strictly_concave", "3#-1,0,0,-1,0,-1", None, None),
)


def load_cases() -> list[tuple[str, str, int | None, dict[str, int] | None]]:
    placeholders = ",".join("?" for _ in DATABASE_IDS)
    with sqlite3.connect(DATABASE) as connection:
        rows = connection.execute(
            f"SELECT matrix_id, dimension, matrix, ess_count, ess_structure "
            f"FROM matrices WHERE matrix_id IN ({placeholders}) ORDER BY matrix_id",
            DATABASE_IDS,
        ).fetchall()
    if len(rows) != len(DATABASE_IDS):
        raise RuntimeError("verification matrix is missing from the canonical database")
    return list(DIRECT_CASES) + [
        (f"matrix_{matrix_id}", f"{dimension}#{matrix}", ess_count, json.loads(ess_structure))
        for matrix_id, dimension, matrix, ess_count, ess_structure in rows
    ]


def worker(module_dir: Path, target_seconds: float) -> None:
    sys.path.insert(0, str(module_dir))
    import fracessa_core  # type: ignore[import-not-found]

    results = []
    for name, matrix, expected_ess_count, expected_ess_structure in load_cases():
        result = fracessa_core.compute_matrix("safe", matrix)
        elapsed = [result["elapsed_ns"]]
        iterations = max(1, math.ceil(target_seconds * 1_000_000_000 / max(1, elapsed[0])))
        for _ in range(1, iterations):
            repeated = fracessa_core.compute_matrix("safe", matrix)
            if repeated["ess_count"] != result["ess_count"] or repeated["ess_structure"] != result["ess_structure"]:
                raise RuntimeError(f"non-deterministic ESS output for {name}")
            elapsed.append(repeated["elapsed_ns"])
        results.append(
            {
                "name": name,
                "status": result["status"],
                "error_message": result["error_message"],
                "ess_count": result["ess_count"],
                "ess_structure": result["ess_structure"],
                "elapsed_ns": statistics.median_low(elapsed),
                "iterations": iterations,
                "expected_ess_count": expected_ess_count,
                "expected_ess_structure": expected_ess_structure,
            }
        )
    print(json.dumps(results, sort_keys=True))


def run(module_dir: Path, target_seconds: float) -> list[dict]:
    configured_python = sys.executable
    cache = BUILD / "CMakeCache.txt"
    if cache.exists():
        for line in cache.read_text().splitlines():
            if line.startswith("_Python3_EXECUTABLE:INTERNAL="):
                configured_python = line.split("=", 1)[1]
                break
    process = subprocess.run(
        [configured_python, str(MODEL_DIR / "verify.py"), "--worker", str(module_dir), "--target-seconds", str(target_seconds)],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(process.stdout)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--target-seconds", type=float, default=0.0)
    arguments = parser.parse_args()
    if arguments.target_seconds < 0:
        parser.error("--target-seconds must be non-negative")
    if arguments.worker:
        worker(arguments.worker, arguments.target_seconds)
        return

    baseline = run(BUILD / "baseline", arguments.target_seconds)
    model = run(BUILD / "a2", arguments.target_seconds)
    failed = False
    print("case\tbaseline_ns\ta2_ns\tratio\tess\tbaseline_n\ta2_n")
    for expected, actual in zip(baseline, model, strict=True):
        same_baseline = all(
            actual[key] == expected[key] for key in ("name", "status", "error_message", "ess_count", "ess_structure")
        )
        same_database = actual["expected_ess_count"] is None or (
            actual["ess_count"] == actual["expected_ess_count"]
            and actual["ess_structure"] == actual["expected_ess_structure"]
        )
        failed |= not same_baseline or not same_database
        ratio = actual["elapsed_ns"] / expected["elapsed_ns"] if expected["elapsed_ns"] else float("nan")
        print(
            f'{actual["name"]}\t{expected["elapsed_ns"]}\t{actual["elapsed_ns"]}\t{ratio:.3f}\t{actual["ess_count"]}'
            f'\t{expected["iterations"]}\t{actual["iterations"]}'
        )
    if failed:
        raise SystemExit("A2 ESS output differs from the production baseline")


if __name__ == "__main__":
    main()
