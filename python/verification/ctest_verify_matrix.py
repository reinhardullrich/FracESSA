#!/usr/bin/env python3
"""
Per-matrix correctness verifier for CTest.

Checks:
1. ESS count equals expected `number_ess` from verification_matrices.json.
2. Candidate rows equal baseline_candidates.csv for this matrix.
   Comparison includes deterministic candidate_id.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


def load_matrices(verification_dir: Path):
    matrices_file = verification_dir / "verification_matrices.json"
    with matrices_file.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    matrices = {int(m["id"]): m for m in data.get("matrices", [])}
    return matrices


def parse_bool(value):
    value_str = str(value).strip().lower()
    return value_str in {"1", "true", "yes"}


def normalize_reason(reason: str):
    """
    Canonicalize reason strings for correctness comparison.

    Policy:
    - Treat legacy `T_pd_rat` as `T_pd_frc`.
    - Treat all positive-definite variants as equivalent (`T_pd_*`)
      because dbl/frc branch choice can vary with numerical thresholds.
    """
    normalized = str(reason).strip()
    if normalized == "T_pd_rat":
        normalized = "T_pd_frc"
    if normalized in {"T_pd_dbl", "T_pd_frc"}:
        return "T_pd_*"
    return normalized


def normalize_candidate_row(row):
    payoff_double_raw = row["payoff_double"]
    try:
        payoff_double = float(payoff_double_raw)
    except (TypeError, ValueError):
        payoff_double = 0.0

    return {
        "candidate_id": int(row["candidate_id"]),
        "vector": str(row["vector"]),
        "support": int(row["support"]),
        "support_size": int(row["support_size"]),
        "extended_support": int(row["extended_support"]),
        "extended_support_size": int(row["extended_support_size"]),
        "multiplier": int(row["multiplier"]) if row["multiplier"] else None,
        "is_ess": parse_bool(row["is_ess"]),
        "reason_ess": normalize_reason(row["reason_ess"]),
        "payoff": str(row["payoff"]),
        "payoff_double": payoff_double,
    }


def parse_cli_candidates(output_lines):
    candidates = []
    if len(output_lines) < 3:
        return candidates

    for line in output_lines[3:]:
        if not line.strip():
            continue
        parts = line.split(";")
        if len(parts) < 11:
            continue
        payoff_double = 0.0
        if len(parts) > 10 and parts[10]:
            try:
                payoff_double = float(parts[10])
            except ValueError:
                payoff_double = 0.0

        candidates.append(
            {
                "candidate_id": int(parts[0]) if parts[0] else 0,
                "vector": parts[1],
                "support": int(parts[2]) if parts[2] else 0,
                "support_size": int(parts[3]) if parts[3] else 0,
                "extended_support": int(parts[4]) if parts[4] else 0,
                "extended_support_size": int(parts[5]) if parts[5] else 0,
                "multiplier": int(parts[6]) if parts[6] else None,
                "is_ess": parse_bool(parts[7]) if len(parts) > 7 else False,
                "reason_ess": normalize_reason(parts[8] if len(parts) > 8 else ""),
                "payoff": parts[9] if len(parts) > 9 else "0",
                "payoff_double": payoff_double,
            }
        )
    return candidates


def run_fracessa(fracessa_exe: Path, matrix_cli: str):
    cmd = [str(fracessa_exe), "-c", "-t", matrix_cli]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"fracessa failed ({result.returncode}): {result.stderr.strip()}")

    lines = [line.strip() for line in result.stdout.strip().splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("fracessa produced no output")

    try:
        ess_count = int(lines[0])
    except ValueError as exc:
        raise RuntimeError(f"failed parsing ESS count from '{lines[0]}'") from exc

    return ess_count, parse_cli_candidates(lines)


def load_baseline_candidates(verification_dir: Path, matrix_id: int):
    baseline_file = verification_dir / "baseline_candidates.csv"
    if not baseline_file.exists():
        raise RuntimeError(f"missing baseline CSV: {baseline_file}")

    rows = []
    with baseline_file.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if int(row["matrix_id"]) == matrix_id:
                rows.append(normalize_candidate_row(row))
    return rows


def comparable_key(candidate):
    return (
        candidate["candidate_id"],
        candidate["vector"],
        candidate["support"],
        candidate["support_size"],
        candidate["extended_support"],
        candidate["extended_support_size"],
        candidate["multiplier"],
        candidate["is_ess"],
        candidate["reason_ess"],
        candidate["payoff"],
        round(candidate["payoff_double"], 9),
    )


def main():
    parser = argparse.ArgumentParser(description="Per-matrix correctness check for CTest.")
    parser.add_argument("--matrix-id", type=int, required=True, help="Matrix ID from verification dataset.")
    parser.add_argument("--verification-dir", type=Path, required=True, help="Path to python/verification.")
    parser.add_argument("--fracessa-exe", type=Path, required=True, help="Path to fracessa executable.")
    args = parser.parse_args()

    matrices = load_matrices(args.verification_dir)
    matrix = matrices.get(args.matrix_id)
    if matrix is None:
        print(f"[ERROR] matrix id {args.matrix_id} not found in verification_matrices.json")
        return 1

    matrix_cli = f"{matrix['dimension']}#{matrix['matrix']}"
    expected_ess = int(matrix["number_ess"])

    try:
        ess_count, produced_candidates = run_fracessa(args.fracessa_exe, matrix_cli)
    except Exception as exc:
        print(f"[ERROR] matrix {args.matrix_id}: execution failed: {exc}")
        return 1

    if ess_count != expected_ess:
        print(
            f"[ERROR] matrix {args.matrix_id}: ESS mismatch "
            f"(expected {expected_ess}, got {ess_count})"
        )
        return 1

    produced_norm = [normalize_candidate_row(row) for row in produced_candidates]
    baseline_norm = load_baseline_candidates(args.verification_dir, args.matrix_id)

    produced_norm.sort(key=lambda row: row["candidate_id"])
    baseline_norm.sort(key=lambda row: row["candidate_id"])

    produced_keys = [comparable_key(row) for row in produced_norm]
    baseline_keys = [comparable_key(row) for row in baseline_norm]

    if produced_keys != baseline_keys:
        print(
            f"[ERROR] matrix {args.matrix_id}: candidate mismatch "
            f"(produced {len(produced_keys)}, baseline {len(baseline_keys)})"
        )
        produced_set = set(produced_keys)
        baseline_set = set(baseline_keys)
        extra = list(produced_set - baseline_set)
        missing = list(baseline_set - produced_set)
        if extra:
            print(f"  extra sample: {extra[0]}")
        if missing:
            print(f"  missing sample: {missing[0]}")
        return 1

    print(f"[OK] matrix {args.matrix_id}: ess and candidates match baseline")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
