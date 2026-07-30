#!/usr/bin/env python3
"""
Script to regenerate baseline files by running the modified executable (REF_linux64_modified.exe)
on all matrices from verification_matrices.json.

This script:
1. Loads matrices from verification_matrices.json
2. Runs the modified executable (python/the_old_exe/REF_linux64_modified.exe) with -c and -t flags
3. Parses the output and collapses circular candidates to bracelet representatives
4. Writes both baseline files:
   - baseline_result.json
   - baseline_candidates.csv
"""

import json
import csv
import subprocess
import sys
import time
import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def map_reason_ess_to_string(reason_ess_value: str) -> str:
    """
    Normalize reason_ess values to canonical string values.
    
    Old enum values (numbers):
    1 = true_pure_ess -> "T_pure_ess"
    2 = true_posdef_double -> "T_pd_dbl"
    3 = true_posdef_rational -> "T_pd_frc"
    4 = true_copositive -> "T_copos"
    5 = false_not_posdef_and_kay_0_1 -> "F_not_pd_kay_0_1"
    6 = false_not_partial_copositive -> "F_not_part_copos"
    7 = false_not_copositive -> "F_not_copos"
    """
    normalized = str(reason_ess_value).strip()

    numeric_mapping = {
        '1': 'T_pure_ess',
        '2': 'T_pd_dbl',
        '3': 'T_pd_frc',
        '4': 'T_copos',
        '5': 'F_not_pd_kay_0_1',
        '6': 'F_not_part_copos',
        '7': 'F_not_copos'
    }
    legacy_string_mapping = {
        'T_pd_double': 'T_pd_dbl',
        'T_pd_rat': 'T_pd_frc',
    }

    canonical_values = set(numeric_mapping.values())
    if normalized in canonical_values:
        return normalized
    if normalized in legacy_string_mapping:
        return legacy_string_mapping[normalized]
    return numeric_mapping.get(normalized, normalized)


def rotate_support_right(support: int, dimension: int) -> int:
    mask = (1 << dimension) - 1
    return ((support >> 1) | ((support & 1) << (dimension - 1))) & mask


def reflect_support(support: int, dimension: int) -> int:
    reflected = 0
    for bit in range(dimension):
        reflected |= ((support >> bit) & 1) << (dimension - 1 - bit)
    return reflected


def dihedral_orbit(support: int, dimension: int) -> set[int]:
    orbit = set()
    current = support
    while current not in orbit:
        orbit.add(current)
        current = rotate_support_right(current, dimension)

    current = reflect_support(support, dimension)
    while current not in orbit:
        orbit.add(current)
        current = rotate_support_right(current, dimension)
    return orbit


def compress_candidates(candidates: List[Dict], dimension: int, is_cs: bool) -> List[Dict]:
    if not is_cs:
        result = []
        for candidate_id, candidate in enumerate(candidates, start=1):
            representative = candidate.copy()
            representative['candidate_id'] = candidate_id
            representative['multiplier'] = None
            result.append(representative)
        return result

    groups = {}
    for candidate in candidates:
        orbit = dihedral_orbit(candidate['support'], dimension)
        canonical_support = min(orbit)
        group = groups.setdefault(
            canonical_support,
            {'multiplier': len(orbit), 'candidates': []},
        )
        group['candidates'].append(candidate)

    result = []
    for canonical_support, group in groups.items():
        matches = [
            candidate for candidate in group['candidates']
            if candidate['support'] == canonical_support
        ]
        if not matches:
            raise ValueError(
                f"legacy output has no canonical row for support {canonical_support}"
            )
        representative = min(matches, key=lambda candidate: candidate['candidate_id']).copy()
        representative['multiplier'] = group['multiplier']
        result.append(representative)

    result.sort(key=lambda candidate: (candidate['support_size'], candidate['support']))
    for candidate_id, representative in enumerate(result, start=1):
        representative['candidate_id'] = candidate_id
    return result


def full_matrix_to_upper_triangular(matrix_str: str, dimension: int) -> str:
    """
    Convert a full matrix (n*n elements in row-major order) to upper triangular format.
    
    For an n×n matrix, upper triangular includes elements where i <= j.
    In row-major order: row i starts at index i*n, and we take elements from column i to n-1.
    
    Args:
        matrix_str: Comma-separated string of n*n elements (row-major order)
        dimension: Matrix dimension n
        
    Returns:
        Comma-separated string of upper triangular elements
    """
    elements = matrix_str.split(',')
    if len(elements) != dimension * dimension:
        raise ValueError(f"Expected {dimension * dimension} elements, got {len(elements)}")
    
    upper_tri_elements = []
    for i in range(dimension):
        # Row i: take elements from column i to n-1
        row_start = i * dimension
        for j in range(i, dimension):
            upper_tri_elements.append(elements[row_start + j])
    
    return ','.join(upper_tri_elements)


def find_old_executable(script_dir: Path) -> Path:
    """Find the modified executable path."""
    # Try relative to script directory (python/verification/)
    # The modified executable is in the_old_exe subdirectory of verification
    old_exe_path = script_dir / "the_old_exe" / "REF_linux64_modified.exe"
    
    if old_exe_path.exists():
        return old_exe_path
    
    # Try from project root
    project_root = script_dir.parent.parent
    old_exe_path = project_root / "python" / "verification" / "the_old_exe" / "REF_linux64_modified.exe"
    if old_exe_path.exists():
        return old_exe_path
    
    raise FileNotFoundError(f"Modified executable not found. Tried: {old_exe_path}")


def run_old_executable(exe_path: Path, matrix_cli_string: str, matrix_id: int = -1) -> Tuple[bool, int, List[Dict], float, Optional[str]]:
    """
    Run the modified executable and parse output.
    The modified executable supports -t flag for C++ timing output.
    Note: matrix_id is kept for tracking but not passed to the executable (it doesn't support -m flag).
    
    Args:
        exe_path: Path to the executable
        matrix_cli_string: Matrix in CLI format (e.g., "2#0,1,1,0")
        matrix_id: Matrix ID (kept for tracking, not passed to executable)
    
    Returns:
        (success, ess_count, candidates, timing, error_message)
    """
    # Modified executable supports -t flag for timing (but not -m flag for matrix ID)
    cmd = [str(exe_path), "-c", "-t", matrix_cli_string]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            return (False, 0, [], 0.0, f"Command failed with return code {result.returncode}: {result.stderr}")
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        if not lines or not lines[0].strip():
            return (False, 0, [], 0.0, "No output from executable")
        
        try:
            ess_count = int(lines[0].strip())
        except ValueError:
            return (False, 0, [], 0.0, f"Could not parse ESS count from: {lines[0]}")
        
        # Parse timing from line 1 (C++ timing output)
        timing = 0.0
        if len(lines) > 1:
            try:
                timing = float(lines[1].strip())
            except (ValueError, IndexError):
                # If timing line is missing or invalid, use 0.0
                timing = 0.0
        
        # Parse candidates
        candidates = []
        # Modified executable output format with -t flag:
        # Line 0: ESS count
        # Line 1: Timing (float)
        # Line 2: CSV header
        # Lines 3+: CSV data rows
        if len(lines) > 3:
            # Skip line 0 (ESS count), line 1 (timing), and line 2 (header), process remaining lines
            candidate_lines = lines[3:]
            
            for line in candidate_lines:
                if line.strip():
                    parts = line.split(';')
                    if len(parts) >= 11:
                        # Parse is_ess: could be "1", "True", "true", etc.
                        is_ess_str = parts[7].strip().lower() if len(parts) > 7 else "0"
                        is_ess = is_ess_str in ("1", "true", "yes")
                        
                        # Map old reason_ess number to new string format
                        reason_ess_raw = parts[8] if len(parts) > 8 else ''
                        reason_ess_string = map_reason_ess_to_string(reason_ess_raw)
                        
                        candidate_data = {
                            'candidate_id': int(parts[0]) if parts[0] else 0,
                            'vector': parts[1] if parts[1] else '',  # Keep as string (comma-separated)
                            'support': int(parts[2]) if parts[2] else 0,
                            'support_size': int(parts[3]) if parts[3] else 0,
                            'extended_support': int(parts[4]) if parts[4] else 0,
                            'extended_support_size': int(parts[5]) if parts[5] else 0,
                            'is_ess': is_ess,
                            'reason_ess': reason_ess_string,
                            'payoff': parts[9] if len(parts) > 9 else '0',
                            'payoff_double': float(parts[10]) if len(parts) > 10 and parts[10] else 0.0
                        }
                        candidates.append(candidate_data)
        
        return (True, ess_count, candidates, timing, None)
        
    except Exception as e:
        return (False, 0, [], 0.0, f"Exception: {str(e)}")


def process_matrix(matrix_data: Dict, old_exe_path: Path) -> Tuple[Dict, List[Dict]]:
    """
    Process a single matrix.
    
    Returns:
        (result_entry, candidates_list)
    """
    matrix_id = matrix_data['id']
    dimension = matrix_data['dimension']
    number_ess = matrix_data['number_ess']
    is_cs = matrix_data['is_cs']
    # Use matrix_old if available (contains original format), otherwise fall back to matrix
    matrix_str = matrix_data.get('matrix_old', matrix_data['matrix'])
    
    # Build CLI string: dimension#matrix_string
    # For circular symmetric: matrix_str has n/2 elements (compact format)
    # For non-circular symmetric: matrix_str has n*n elements (full matrix in row-major order from matrix_old)
    cli_string = f"{dimension}#{matrix_str}"
    
    # Run old executable
    success, ess_count, candidates, timing, error = run_old_executable(old_exe_path, cli_string, matrix_id)

    if success:
        candidates = compress_candidates(candidates, dimension, is_cs)
        # For output, use the new format from 'matrix' field if available
        # (matrix_old was used for CLI string, but output should use the new format)
        if 'matrix' in matrix_data:
            matrix_str_for_output = matrix_data['matrix']
        else:
            # Fallback: convert matrix_old to appropriate format
            if not is_cs:
                matrix_str_for_output = full_matrix_to_upper_triangular(matrix_str, dimension)
            else:
                matrix_str_for_output = matrix_str
        
        result_entry = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": number_ess,
            "is_cs": is_cs,
            "matrix": matrix_str_for_output,
            "result": {
                "success": True,
                "ess_count": ess_count,
                "timing": timing
            }
        }
        
        # Add matrix_id to each candidate
        candidates_list = []
        for candidate in candidates:
            candidate_with_matrix_id = {"matrix_id": matrix_id}
            candidate_with_matrix_id.update(candidate)
            candidates_list.append(candidate_with_matrix_id)
        
        return (result_entry, candidates_list)
    else:
        # For output, use the new format from 'matrix' field if available
        if 'matrix' in matrix_data:
            matrix_str_for_output = matrix_data['matrix']
        else:
            # Fallback: convert matrix_old to appropriate format
            if not is_cs:
                matrix_str_for_output = full_matrix_to_upper_triangular(matrix_str, dimension)
            else:
                matrix_str_for_output = matrix_str
        
        result_entry = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": number_ess,
            "is_cs": is_cs,
            "matrix": matrix_str_for_output,
            "result": {
                "success": False,
                "error": error,
                "timing": timing
            }
        }
        return (result_entry, [])


def main():
    """Main function."""
    script_dir = Path(__file__).resolve().parent
    matrices_file = script_dir / "verification_matrices.json"
    baseline_results_file = script_dir / "baseline_result.json"
    baseline_candidates_file = script_dir / "baseline_candidates.csv"
    
    # Find modified executable
    try:
        old_exe_path = find_old_executable(script_dir)
        print(f"Using modified executable: {old_exe_path}")
        
        # Make executable if needed
        os.chmod(old_exe_path, 0o755)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Load matrices
    if not matrices_file.exists():
        print(f"Error: {matrices_file} not found")
        sys.exit(1)
    
    print(f"Loading matrices from {matrices_file}")
    with open(matrices_file, 'r') as f:
        data = json.load(f)
    
    matrices = data.get('matrices', [])
    
    # Baseline regeneration always processes all matrices in verification_matrices.json.
    total_count = len(matrices)
    print(f"Found {total_count} matrices total")
    print(f"Processing all {total_count} matrices")
    print()
    
    # Process matrices sequentially
    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)
    
    start_time = time.time()
    
    for i, matrix_data in enumerate(matrices):
        matrix_id = matrix_data['id']
        dimension = matrix_data['dimension']
        
        print(f"[{i+1}/{total_matrices}] Processing Matrix ID {matrix_id} (dim {dimension})...", end=" ", flush=True)
        
        result_entry, candidates_list = process_matrix(matrix_data, old_exe_path)
        results.append(result_entry)
        all_candidates.extend(candidates_list)
        
        if result_entry["result"]["success"]:
            actual_ess = result_entry["result"]["ess_count"]
            timing = result_entry["result"]["timing"]
            expected_ess = matrix_data['number_ess']
            status = f"✅ {actual_ess} ESS in {timing:.6f}s"
            if actual_ess != expected_ess:
                status += f" ⚠️ (expected {expected_ess})"
            print(status)
            successful += 1
        else:
            error_msg = result_entry["result"].get('error', 'Unknown error')
            print(f"❌ {error_msg}")
            failed += 1
    
    elapsed_time = time.time() - start_time
    print()
    print(f"Processing completed in {elapsed_time:.2f}s")
    
    # Sort results by matrix ID for consistent output
    results.sort(key=lambda x: x['id'])
    all_candidates.sort(key=lambda x: (x['matrix_id'], x['candidate_id']))
    
    # Write JSON baseline file
    print(f"\nWriting results to {baseline_results_file}")
    with open(baseline_results_file, 'w') as f:
        json.dump({
            "metadata": {
                "total_matrices": len(matrices),
                "successful": successful,
                "failed": failed,
                "timestamp": time.time(),
                "processing_time": elapsed_time,
                "num_processes": 1,  # Sequential processing
                "fracessa_settings": {
                    "include_candidates": True,
                    "enable_logging": False,
                    "exact_arithmetic": False,
                    "full_support_search": False,
                    "use_cpp_timing": True
                }
            },
            "results": results
        }, f, indent=2)
    
    # Write CSV baseline file
    print(f"Writing candidates to {baseline_candidates_file}")
    csv_columns = [
        "matrix_id", "candidate_id", "vector", "support", "support_size",
        "extended_support", "extended_support_size", "multiplier",
        "is_ess", "reason_ess", "payoff", "payoff_double"
    ]
    
    with open(baseline_candidates_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns, lineterminator='\n')
        writer.writeheader()
        for candidate in all_candidates:
            # Convert is_ess boolean to string representation for CSV
            # reason_ess is already converted to string format by map_reason_ess_to_string
            row = candidate.copy()
            row['is_ess'] = 'True' if row.get('is_ess', False) else 'False'
            # Ensure reason_ess is a string (should already be converted, but double-check)
            if 'reason_ess' in row:
                row['reason_ess'] = map_reason_ess_to_string(str(row['reason_ess']))
            writer.writerow(row)
    
    print("\n✅ Baseline files created successfully!")
    print(f"  JSON: {baseline_results_file}")
    print(f"  CSV:  {baseline_candidates_file}")
    print(f"\nSummary: {successful} successful, {failed} failed out of {len(matrices)} matrices")
    print(f"Total candidates: {len(all_candidates)}")
    
    if failed > 0:
        print("\nFailed matrices:")
        for result in results:
            if not result["result"]["success"]:
                print(f"  ID {result['id']}: {result['result']['error']}")


if __name__ == "__main__":
    main()
