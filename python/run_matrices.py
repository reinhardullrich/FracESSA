#!/usr/bin/env python3
"""
Speed benchmarking script for FRACESSA Python wrapper runs.

Purpose:
- execute selected verification matrices,
- record timings/results,
- compare speed against baseline and latest run,
- write timestamped JSON/CSV outputs.

Correctness (ESS/candidate equality) is handled by CTest matrix verification tests.
"""

import argparse
import csv
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from fracessa_py import Fracessa, Matrix, FracessaError
from verification.matrix_selection import FAST_EXCLUDED_IDS


def process_matrix(matrix_data):
    """
    Process one matrix and keep the best timing across a small number of repeats
    for low-dimensional noisy cases.
    """
    matrix_id = matrix_data["id"]
    dimension = matrix_data["dimension"]
    number_ess = matrix_data["number_ess"]
    is_cs = matrix_data["is_cs"]
    matrix_str = matrix_data["matrix"]

    if dimension < 10:
        iterations = 20
    elif dimension < 15:
        iterations = 5
    else:
        iterations = 1

    best_result = None
    min_time = float("inf")

    try:
        fracessa = Fracessa()
        matrix = Matrix(matrix_str, dimension, is_circular=is_cs)

        for _ in range(iterations):
            result = fracessa.compute_ess(
                matrix=matrix,
                include_candidates=True,
                enable_logging=False,
                matrix_id=matrix_id,
            )

            if not result.success:
                return (
                    matrix_id,
                    dimension,
                    number_ess,
                    {
                        "id": matrix_id,
                        "dimension": dimension,
                        "number_ess_expected": number_ess,
                        "is_cs": is_cs,
                        "matrix": matrix_str,
                        "result": {
                            "success": False,
                            "error": result.error,
                            "timing": result.computation_time,
                        },
                    },
                    [],
                    {"success": False, "error": result.error},
                )

            if result.computation_time < min_time:
                min_time = result.computation_time
                best_result = result

        if best_result and best_result.success:
            candidates_data = []
            for candidate in best_result.candidates:
                candidate_dict = {
                    "candidate_id": candidate.candidate_id,
                    "vector": candidate.vector,
                    "support": candidate.support,
                    "support_size": candidate.support_size,
                    "extended_support": candidate.extended_support,
                    "extended_support_size": candidate.extended_support_size,
                    "shift_reference": candidate.shift_reference,
                    "is_ess": candidate.is_ess,
                    "reason_ess": str(candidate.reason_ess),
                    "payoff": candidate.payoff,
                    "payoff_double": candidate.payoff_double,
                }
                candidates_data.append(candidate_dict)

            computation_result = {
                "success": True,
                "ess_count": best_result.ess_count,
                "timing": min_time,
                "candidates": candidates_data,
            }
        else:
            computation_result = {
                "success": False,
                "error": "No valid result found",
                "timing": 0.0,
            }

    except FracessaError as exc:
        computation_result = {
            "success": False,
            "error": f"FRACESSA Error: {str(exc)}",
            "timing": 0.0,
        }
    except Exception as exc:
        computation_result = {
            "success": False,
            "error": f"Exception: {str(exc)}",
            "timing": 0.0,
        }

    candidates_list = []
    if computation_result["success"] and "candidates" in computation_result:
        for candidate in computation_result["candidates"]:
            candidate_with_matrix_id = {"matrix_id": matrix_id}
            candidate_with_matrix_id.update(candidate)
            candidates_list.append(candidate_with_matrix_id)

    result_without_candidates = {
        key: value for key, value in computation_result.items() if key != "candidates"
    }

    result_entry = {
        "id": matrix_id,
        "dimension": dimension,
        "number_ess_expected": number_ess,
        "is_cs": is_cs,
        "matrix": matrix_str,
        "result": result_without_candidates,
    }

    return (
        matrix_id,
        dimension,
        number_ess,
        result_entry,
        candidates_list,
        computation_result,
    )


def load_baseline_timings(baseline_results_file: Path):
    timings = {}
    if not baseline_results_file.exists():
        return timings

    try:
        with open(baseline_results_file, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        for entry in data.get("results", []):
            result = entry.get("result", {})
            if result.get("success"):
                timings[entry["id"]] = result.get("timing", 0.0)
    except Exception:
        pass

    return timings


def get_latest_timings(results_dir: Path):
    """Find latest previous run timings from python/results."""
    if not results_dir.exists():
        return {}, None

    files = sorted(results_dir.glob("fracessa_verification_result_*.json"), key=lambda x: x.name)
    if not files:
        return {}, None

    latest_file = files[-1]
    timings = {}
    try:
        with open(latest_file, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        for result in data.get("results", []):
            if result.get("result", {}).get("success"):
                timings[result["id"]] = result["result"].get("timing", 0.0)
        return timings, latest_file.name
    except Exception:
        return {}, None


def select_matrices(matrices, include_all: bool):
    """
    Selection policy:
    - full mode: all matrices
    - fast mode: all matrices except hardcoded excluded IDs
    """
    matrices_sorted = sorted(matrices, key=lambda m: m.get("id", 0))
    if include_all:
        return matrices_sorted

    selected = []
    for matrix in matrices_sorted:
        matrix_id = int(matrix.get("id", 0))
        if matrix_id not in FAST_EXCLUDED_IDS:
            selected.append(matrix)
    return selected


def format_comparison(current_time, ref_time):
    if ref_time <= 0:
        return f"(was {ref_time:.6f}s)"

    pct = ((current_time - ref_time) / ref_time) * 100
    if current_time > ref_time:
        return f"(was {ref_time:.6f}s, +{pct:.2f}%)"
    if current_time < ref_time:
        return f"(was {ref_time:.6f}s, {pct:.2f}%)"
    return f"(was {ref_time:.6f}s, 0.00%)"


def main():
    parser = argparse.ArgumentParser(description="Run speed benchmarks on verification matrices.")
    parser.add_argument("--full", action="store_true", help="Process all matrices.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    matrices_file = script_dir / "verification" / "verification_matrices.json"
    baseline_results_file = script_dir / "verification" / "baseline_result.json"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = script_dir / "results"
    results_dir.mkdir(exist_ok=True)
    results_file = results_dir / f"fracessa_verification_result_{timestamp}.json"
    candidates_file = results_dir / f"fracessa_verification_candidates_{timestamp}.csv"

    if not matrices_file.exists():
        print(f"Error: {matrices_file} not found")
        return 1

    try:
        import fracessa_py  # noqa: F401
        print("FRACESSA Python bindings loaded successfully")
    except ImportError as exc:
        print(f"Error: Cannot import fracessa_py: {exc}")
        return 1

    with open(matrices_file, "r", encoding="utf-8") as fh:
        matrices_data = json.load(fh)

    baseline_timings = load_baseline_timings(baseline_results_file)
    latest_timings, latest_filename = get_latest_timings(results_dir)

    if latest_filename:
        print(f"Comparing against Previous Run: {latest_filename}")

    matrices = select_matrices(
        matrices_data.get("matrices", []),
        include_all=args.full,
    )

    if args.full:
        print(f"Mode: FULL (all matrices). Processing {len(matrices)} matrices")
    else:
        excluded = ", ".join(str(matrix_id) for matrix_id in sorted(FAST_EXCLUDED_IDS))
        print(f"Mode: FAST (all except IDs: {excluded}). Processing {len(matrices)} matrices")

    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)

    print("\nStarting sequential processing...")
    start_time = time.time()

    for i, matrix_data in enumerate(matrices, start=1):
        matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result = process_matrix(matrix_data)

        results.append(result_entry)
        all_candidates.extend(candidates_list)

        if computation_result["success"]:
            actual_ess = computation_result["ess_count"]
            timing = computation_result["timing"]

            status = f"✅ {actual_ess} ESS in {timing:.6f}s"

            if matrix_id in baseline_timings:
                status += f" | Old: {format_comparison(timing, baseline_timings[matrix_id])}"
            if matrix_id in latest_timings:
                status += f" | Prev: {format_comparison(timing, latest_timings[matrix_id])}"
            if actual_ess != number_ess:
                status += f" ⚠️ (expected {number_ess})"

            print(f"[{i}/{total_matrices}] Matrix ID {matrix_id} (dim {dimension}): {status}")
            successful += 1
        else:
            error_msg = computation_result.get("error", "Unknown error")
            print(f"[{i}/{total_matrices}] Matrix ID {matrix_id}: ❌ {error_msg}")
            failed += 1

    elapsed_time = time.time() - start_time
    print(f"\nProcessing completed in {elapsed_time:.2f}s")

    results.sort(key=lambda x: x["id"])
    all_candidates.sort(key=lambda x: (x["matrix_id"], x["candidate_id"]))

    with open(results_file, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "metadata": {
                    "total_matrices": len(matrices),
                    "successful": successful,
                    "failed": failed,
                    "timestamp": time.time(),
                    "processing_time": elapsed_time,
                    "mode": "full" if args.full else "fast",
                    "fast_excluded_ids": sorted(FAST_EXCLUDED_IDS),
                },
                "results": results,
            },
            fh,
            indent=2,
        )

    csv_columns = [
        "matrix_id",
        "candidate_id",
        "vector",
        "support",
        "support_size",
        "extended_support",
        "extended_support_size",
        "shift_reference",
        "is_ess",
        "reason_ess",
        "payoff",
        "payoff_double",
    ]

    with open(candidates_file, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_columns)
        writer.writeheader()
        for candidate in all_candidates:
            row = candidate.copy()
            if "vector" in row and isinstance(row["vector"], list):
                row["vector"] = ",".join(str(v) for v in row["vector"])
            writer.writerow(row)

    print(f"json saved: {results_file}")
    print(f"csv saved: {candidates_file}")

    if failed == 0:
        print("✅ SPEED RUN COMPLETED")
        return 0

    print(f"❌ SPEED RUN FAILED ({failed} matrices failed)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
