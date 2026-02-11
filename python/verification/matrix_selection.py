#!/usr/bin/env python3
"""
Matrix selection helper for verification workflows.

Modes:
- --list-tests: print one line per matrix as "<id>,<is_fast_flag>"
- --fast-ids: print space-separated matrix IDs for fast correctness runs
- --full-ids: print space-separated matrix IDs for full correctness runs
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

FAST_EXCLUDED_IDS = {32, 34}


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def load_matrices(verification_dir: Path):
    matrices_file = verification_dir / "verification_matrices.json"
    data = load_json(matrices_file)
    matrices = data.get("matrices", [])
    matrices.sort(key=lambda m: int(m.get("id", 0)))
    return matrices


def select_fast_ids(matrices):
    all_ids = [int(matrix["id"]) for matrix in matrices]
    return [matrix_id for matrix_id in all_ids if matrix_id not in FAST_EXCLUDED_IDS]


def main():
    parser = argparse.ArgumentParser(description="Matrix selection helper.")
    parser.add_argument(
        "--verification-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Path to python/verification directory.",
    )
    parser.add_argument("--list-tests", action="store_true", help="Emit id,is_fast lines.")
    parser.add_argument("--fast-ids", action="store_true", help="Emit fast matrix ids.")
    parser.add_argument("--full-ids", action="store_true", help="Emit all matrix ids.")
    args = parser.parse_args()

    matrices = load_matrices(args.verification_dir)
    fast_ids = set(select_fast_ids(matrices))

    if args.list_tests:
        for matrix in matrices:
            matrix_id = int(matrix["id"])
            is_fast_flag = 1 if matrix_id in fast_ids else 0
            print(f"{matrix_id},{is_fast_flag}")
        return 0

    if args.full_ids:
        ids = [str(int(matrix["id"])) for matrix in matrices]
        print(" ".join(ids))
        return 0

    if args.fast_ids:
        ids = sorted(fast_ids)
        print(" ".join(str(matrix_id) for matrix_id in ids))
        return 0

    parser.error("One of --list-tests, --fast-ids, or --full-ids is required.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
