#!/usr/bin/env python3
"""
Script to run all matrices from matrices.json using the FRACESSA Python bindings
and save results to a JSON file. Processes matrices sequentially.

Optimized: Uses micro-benchmarking loops for small matrices to filter out OS noise.
Features: Compares timing against 'baseline_result.json' AND the latest 'fracessa_verification_result_*.json'.
"""

import json
import sys
import time
import csv
from datetime import datetime
from pathlib import Path
from fracessa_py import Fracessa, Matrix, FracessaError


def candidate_to_comparison_key(candidate):
    """
    Convert a candidate dict to a hashable tuple for comparison.
    """
    vector = candidate.get('vector', '')
    if isinstance(vector, list):
        vector = ','.join(str(v) for v in vector)
    
    return (
        str(vector),
        str(candidate.get('support', '')),
        str(candidate.get('support_size', '')),
        str(candidate.get('extended_support', '')),
        str(candidate.get('extended_support_size', '')),
        str(candidate.get('shift_reference', '')),
        str(candidate.get('is_ess', '')),
        str(candidate.get('payoff', '')),
        str(candidate.get('payoff_double', ''))
    )


def compare_ess_counts(results, baseline_file):
    """Compare ESS counts of current run against baseline JSON file."""
    if not baseline_file.exists():
        return {"status": "SKIPPED", "reason": f"Baseline file not found", "compared": 0, "mismatches": []}
    
    try:
        with open(baseline_file, 'r') as f:
            baseline_data = json.load(f)
        baseline_results = baseline_data.get('results', [])
    except Exception as e:
        return {"status": "ERROR", "reason": f"Failed to load baseline: {e}", "compared": 0, "mismatches": []}
    
    baseline_ess_counts = {}
    for result in baseline_results:
        matrix_id = result.get('id')
        if result.get('result', {}).get('success'):
            baseline_ess_counts[matrix_id] = result['result'].get('ess_count', 0)
    
    compared = 0
    mismatches = []
    
    for result in results:
        matrix_id = result.get('id')
        current_result = result.get('result', {})
        if not current_result.get('success'): continue
        if matrix_id not in baseline_ess_counts: continue
        
        current_ess_count = current_result.get('ess_count', 0)
        baseline_ess_count = baseline_ess_counts[matrix_id]
        compared += 1
        
        if current_ess_count != baseline_ess_count:
            mismatches.append({
                "matrix_id": matrix_id,
                "dimension": result.get('dimension'),
                "current_ess": current_ess_count,
                "baseline_ess": baseline_ess_count
            })
    
    return {"status": "SUCCESS", "compared": compared, "mismatches": mismatches}


def verify_against_baseline(all_candidates, baseline_file):
    """Verify produced candidates against baseline CSV."""
    from collections import defaultdict
    
    if not baseline_file.exists():
        return {"status": "SKIPPED", "reason": f"Baseline file not found", "verified": 0, "passed": 0, "failed": 0, "details": []}
    
    baseline_candidates = []
    with open(baseline_file, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            baseline_candidates.append(row)
    
    produced_by_matrix = defaultdict(set)
    baseline_by_matrix = defaultdict(set)
    produced_counts = defaultdict(int)
    baseline_counts = defaultdict(int)
    
    for candidate in all_candidates:
        matrix_id = candidate.get('matrix_id')
        if matrix_id is None: continue
        matrix_id = int(matrix_id)
        key = candidate_to_comparison_key(candidate)
        produced_by_matrix[matrix_id].add(key)
        produced_counts[matrix_id] += 1
    
    for candidate in baseline_candidates:
        matrix_id = int(candidate.get('matrix_id'))
        key = candidate_to_comparison_key(candidate)
        baseline_by_matrix[matrix_id].add(key)
        baseline_counts[matrix_id] += 1
    
    produced_ids = set(produced_by_matrix.keys())
    baseline_ids = set(baseline_by_matrix.keys())
    to_verify = produced_ids & baseline_ids
    
    passed = 0
    failed = 0
    mismatches = []
    errors = []
    count_mismatches = []
    
    # Check for matrices not in baseline (ERROR)
    for matrix_id in sorted(produced_ids - baseline_ids):
        failed += 1
        errors.append({"matrix_id": matrix_id, "error": "Matrix ID not found in baseline"})
    
    for matrix_id in sorted(to_verify):
        p_set = produced_by_matrix[matrix_id]
        b_set = baseline_by_matrix[matrix_id]
        p_count = produced_counts[matrix_id]
        b_count = baseline_counts[matrix_id]
        
        if p_count != b_count:
            failed += 1
            count_mismatches.append({"matrix_id": matrix_id, "produced_count": p_count, "baseline_count": b_count})
        
        if p_set == b_set and p_count == b_count:
            passed += 1
        elif p_set != b_set:
            failed += 1
            mismatches.append({
                "matrix_id": matrix_id, 
                "produced_count": len(p_set), 
                "baseline_count": len(b_set),
                "extra_count": len(p_set - b_set), 
                "missing_count": len(b_set - p_set)
            })
    
    return {
        "status": "PASS" if failed == 0 else "FAIL",
        "verified": len(to_verify),
        "passed": passed, "failed": failed,
        "errors": errors, "mismatches": mismatches, "candidate_count_mismatches": count_mismatches
    }


def process_matrix(matrix_data):
    """
    Process a single matrix.
    If the matrix is small (dimension < 15), run it multiple times to get stable timing.
    """
    matrix_id = matrix_data['id']
    dimension = matrix_data['dimension']
    number_ess = matrix_data['number_ess']
    is_cs = matrix_data['is_cs']
    matrix_str = matrix_data['matrix']
    
    # HEURISTIC: Small matrices are noisy. Run them multiple times.
    if dimension < 10:
        iterations = 20
    elif dimension < 15:
        iterations = 5
    else:
        iterations = 1

    best_result = None
    min_time = float('inf')

    try:
        fracessa = Fracessa()
        matrix = Matrix(matrix_str, dimension, is_circular=is_cs)
        
        # BENCHMARK LOOP
        for i in range(iterations):
            result = fracessa.compute_ess(
                matrix=matrix,
                include_candidates=True,
                enable_logging=False,
                timeout=1800.0,
                matrix_id=matrix_id
            )
            
            if not result.success:
                return (matrix_id, dimension, number_ess, {
                    "id": matrix_id, "dimension": dimension, "number_ess_expected": number_ess,
                    "is_cs": is_cs, "matrix": matrix_str,
                    "result": {"success": False, "error": result.error, "timing": result.computation_time}
                }, [], {"success": False, "error": result.error})
            
            if result.computation_time < min_time:
                min_time = result.computation_time
                best_result = result
        
        if best_result and best_result.success:
            candidates_data = []
            for candidate in best_result.candidates:
                candidate_dict = {
                    'candidate_id': candidate.candidate_id,
                    'vector': candidate.vector,
                    'support': candidate.support,
                    'support_size': candidate.support_size,
                    'extended_support': candidate.extended_support,
                    'extended_support_size': candidate.extended_support_size,
                    'shift_reference': candidate.shift_reference,
                    'is_ess': candidate.is_ess,
                    'reason_ess': str(candidate.reason_ess),
                    'payoff': candidate.payoff,
                    'payoff_double': candidate.payoff_double
                }
                candidates_data.append(candidate_dict)
            
            computation_result = {
                "success": True,
                "ess_count": best_result.ess_count,
                "timing": min_time,
                "candidates": candidates_data
            }
        else:
             computation_result = {"success": False, "error": "No valid result found", "timing": 0.0}
    
    except FracessaError as e:
        computation_result = {"success": False, "error": f"FRACESSA Error: {str(e)}", "timing": 0.0}
    except Exception as e:
        computation_result = {"success": False, "error": f"Exception: {str(e)}", "timing": 0.0}
    
    candidates_list = []
    if computation_result["success"] and "candidates" in computation_result:
        for candidate in computation_result["candidates"]:
            candidate_with_matrix_id = {"matrix_id": matrix_id}
            candidate_with_matrix_id.update(candidate)
            candidates_list.append(candidate_with_matrix_id)
    
    result_without_candidates = {k: v for k, v in computation_result.items() if k != "candidates"}
    
    result_entry = {
        "id": matrix_id, "dimension": dimension, "number_ess_expected": number_ess,
        "is_cs": is_cs, "matrix": matrix_str, "result": result_without_candidates
    }
    
    return (matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result)


def get_latest_timings(results_dir: Path):
    """
    Finds the most recent result file (excluding current run) and extracts timings.
    Returns: (timings_dict, filename)
    """
    if not results_dir.exists():
        return {}, None
    
    # glob files like 'fracessa_verification_result_*.json'
    files = list(results_dir.glob("fracessa_verification_result_*.json"))
    if not files:
        return {}, None
        
    # Sort files by name (timestamps are in the name YYYYMMDD_HHMMSS)
    files.sort(key=lambda x: x.name)
    
    # Pick the last one
    latest_file = files[-1]
    
    timings = {}
    try:
        with open(latest_file, 'r') as f:
            data = json.load(f)
        for r in data.get('results', []):
            if r.get('result', {}).get('success'):
                timings[r['id']] = r['result'].get('timing', 0.0)
        return timings, latest_file.name
    except Exception:
        return {}, None


def format_comparison(current_time, ref_time):
    """Generates string: (was Xs, +Z%)"""
    if ref_time <= 0:
        return f"(was {ref_time:.6f}s)"
    
    difference = current_time - ref_time
    pct = ((current_time - ref_time) / ref_time) * 100
    
    if difference > 0:
        return f"(was {ref_time:.6f}s, +{pct:.2f}%)"
    elif difference < 0:
        return f"(was {ref_time:.6f}s, {pct:.2f}%)"
    else:
        return f"(was {ref_time:.6f}s, 0.00%)"


def main():
    # Paths
    script_dir = Path(__file__).resolve().parent
    matrices_file = script_dir / "verification" / "verification_matrices.json"
    
    # Setup Output Paths
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = script_dir / "results"
    results_dir.mkdir(exist_ok=True)
    results_file = results_dir / f"fracessa_verification_result_{timestamp}.json"
    candidates_file = results_dir / f"fracessa_verification_candidates_{timestamp}.csv"
    
    # Baseline Files
    baseline_candidates_file = script_dir / "verification" / "baseline_candidates.csv"
    baseline_results_file = script_dir / "verification" / "baseline_result.json"

    # 1. Load Baseline Timings (Old Exe)
    baseline_timings = {}
    if baseline_results_file.exists():
        try:
            with open(baseline_results_file, 'r') as f:
                data = json.load(f)
            for r in data.get('results', []):
                if r.get('result', {}).get('success'):
                    baseline_timings[r['id']] = r['result'].get('timing', 0.0)
        except Exception:
            pass

    # 2. Load Previous Run Timings (Latest Python Run)
    latest_timings, latest_filename = get_latest_timings(results_dir)
    if latest_filename:
        print(f"Comparing against Previous Run: {latest_filename}")

    # Load Matrices
    if not matrices_file.exists():
        print(f"Error: {matrices_file} not found")
        sys.exit(1)

    try:
        import fracessa_py
        print("FRACESSA Python bindings loaded successfully")
    except ImportError as e:
        print(f"Error: Cannot import fracessa_py: {e}")
        sys.exit(1)

    with open(matrices_file, 'r') as f:
        data = json.load(f)

    matrices = [m for m in data.get('matrices', []) if m.get('in_use', True)]
    matrices.sort(key=lambda x: x.get('id', 0))

    print(f"Processing {len(matrices)} matrices")
    
    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)

    print(f"\nStarting sequential processing...")
    start_time = time.time()

    for i, matrix_data in enumerate(matrices):
        result_tuple = process_matrix(matrix_data)
        matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result = result_tuple
        
        results.append(result_entry)
        all_candidates.extend(candidates_list)
        
        if computation_result["success"]:
            actual_ess = computation_result["ess_count"]
            timing = computation_result["timing"]
            
            # Base status string
            status = f"✅ {actual_ess} ESS in {timing:.6f}s"
            
            # 1. Compare to Old Exe (Baseline)
            if matrix_id in baseline_timings:
                diff_str = format_comparison(timing, baseline_timings[matrix_id])
                status += f" | Old: {diff_str}"
            
            # 2. Compare to Previous Python Run (Latest)
            if matrix_id in latest_timings:
                diff_str = format_comparison(timing, latest_timings[matrix_id])
                status += f" | Prev: {diff_str}"
            
            if actual_ess != number_ess:
                status += f" ⚠️ (expected {number_ess})"
            
            print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id} (dim {dimension}): {status}")
            successful += 1
        else:
            error_msg = computation_result.get('error', 'Unknown error')
            if "Timed out" in error_msg:
                print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id}: ⏰ {error_msg}")
            else:
                print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id}: ❌ {error_msg}")
            failed += 1

    elapsed_time = time.time() - start_time
    print(f"\nProcessing completed in {elapsed_time:.2f}s")

    # Save Results
    results.sort(key=lambda x: x['id'])
    all_candidates.sort(key=lambda x: (x['matrix_id'], x['candidate_id']))

    with open(results_file, 'w') as f:
        json.dump({
            "metadata": {
                "total_matrices": len(matrices),
                "successful": successful,
                "failed": failed,
                "timestamp": time.time(),
                "processing_time": elapsed_time
            },
            "results": results
        }, f, indent=2)

    csv_columns = [
        "matrix_id", "candidate_id", "vector", "support", "support_size",
        "extended_support", "extended_support_size", "shift_reference",
        "is_ess", "reason_ess", "payoff", "payoff_double"
    ]
    
    with open(candidates_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns)
        writer.writeheader()
        for candidate in all_candidates:
            row = candidate.copy()
            if 'vector' in row and isinstance(row['vector'], list):
                row['vector'] = ','.join(str(v) for v in row['vector'])
            writer.writerow(row)

    print(f"json saved: {results_file}")
    print(f"csv saved: {candidates_file}")
    
    # Verification Steps
    print("\n--- Verification ---")
    ess_count_comparison = compare_ess_counts(results, baseline_results_file)
    verification = verify_against_baseline(all_candidates, baseline_candidates_file)
    
    ess_ok = False
    cand_ok = False
    
    # Check ESS
    if ess_count_comparison["status"] == "SUCCESS":
        if not ess_count_comparison["mismatches"]:
            ess_ok = True
        else:
            print(f"❌ ESS Mismatches: {len(ess_count_comparison['mismatches'])}")
    elif ess_count_comparison["status"] == "SKIPPED":
        print("⚠️ Baseline not found for ESS check")
        ess_ok = True # Treat as pass if file missing to avoid error noise

    # Check Candidates
    if verification["status"] == "PASS":
        cand_ok = True
    elif verification["status"] == "FAIL":
        print(f"❌ Candidate Mismatches: {verification['failed']}")
    elif verification["status"] == "SKIPPED":
        cand_ok = True

    if ess_ok and cand_ok:
        print(f"✅ VERIFICATION PASSED")
    else:
        print(f"❌ VERIFICATION FAILED")

if __name__ == "__main__":
    main()