#!/usr/bin/env python3
"""
Black-box parser checks against the built `fracessa` executable.

These tests validate CLI parser behavior through the real binary, complementing
unit tests on matrix_parser internals.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


def run_case(fracessa_exe: Path, args: list[str], timeout: float = 30.0) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(fracessa_exe), *args],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def first_non_empty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def assert_success_with_ess_output(
    fracessa_exe: Path,
    args: list[str],
    case_name: str,
) -> subprocess.CompletedProcess:
    result = run_case(fracessa_exe, args)
    if result.returncode != 0:
        raise AssertionError(
            f"{case_name}: expected success, got rc={result.returncode}, stderr={result.stderr.strip()}"
        )

    first_line = first_non_empty_line(result.stdout)
    if not re.fullmatch(r"\d+", first_line):
        raise AssertionError(
            f"{case_name}: expected integer ESS count on first line, got '{first_line}'"
        )
    return result


def assert_unsafe_warning(result: subprocess.CompletedProcess, case_name: str) -> None:
    warning = result.stderr.lower()
    required = (
        "heuristic numerical filtering",
        "miss exact candidates and ess results",
    )
    for text in required:
        if text not in warning:
            raise AssertionError(f"{case_name}: unsafe warning is missing '{text}'")


def assert_failure_with_stderr(
    fracessa_exe: Path,
    args: list[str],
    expected_substring: str,
    case_name: str,
) -> None:
    result = run_case(fracessa_exe, args)
    if result.returncode == 0:
        raise AssertionError(f"{case_name}: expected non-zero exit code")
    if expected_substring not in result.stderr:
        raise AssertionError(
            f"{case_name}: expected stderr to contain '{expected_substring}', got '{result.stderr.strip()}'"
        )


def assert_candidate_header_matches_rows(fracessa_exe: Path) -> None:
    result = run_case(fracessa_exe, ["--candidates", "2#0,1,0"])
    if result.returncode != 0:
        raise AssertionError(
            f"candidate_columns: expected success, got rc={result.returncode}, stderr={result.stderr.strip()}"
        )

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        raise AssertionError(f"candidate_columns: expected count, header, and candidate row; got {lines}")

    header_count = len(lines[1].split(";"))
    for row in lines[2:]:
        row_count = len(row.split(";"))
        if row_count != header_count:
            raise AssertionError(
                f"candidate_columns: header has {header_count} fields, candidate row has {row_count}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Black-box parser checks for fracessa CLI.")
    parser.add_argument("--fracessa-exe", type=Path, required=True, help="Path to built fracessa executable.")
    args = parser.parse_args()

    fracessa_exe = args.fracessa_exe
    if not fracessa_exe.exists():
        print(f"[ERROR] missing executable: {fracessa_exe}")
        return 1

    verification_file = Path(__file__).resolve().parents[2] / "python/verification/verification_matrices.json"
    with verification_file.open("r", encoding="utf-8") as fh:
        matrix_46 = next(matrix for matrix in json.load(fh)["matrices"] if matrix["id"] == 46)
    unsafe_counterexample = f"{matrix_46['dimension']}#{matrix_46['matrix']}"

    # Matrix 46 has one ESS under bounded-error/exact analysis, but unsafe filtering misses it.
    default_result = assert_success_with_ess_output(
        fracessa_exe, [unsafe_counterexample], "default_bounded_success"
    )
    explicit_unsafe_result = assert_success_with_ess_output(
        fracessa_exe, ["--unsafe", unsafe_counterexample], "explicit_unsafe_success"
    )
    if first_non_empty_line(default_result.stdout) != "1":
        raise AssertionError(f"bounded-error mode missed matrix 46 ESS: {default_result.stdout.strip()}")
    if first_non_empty_line(explicit_unsafe_result.stdout) != "0":
        raise AssertionError(f"unsafe mode unexpectedly retained matrix 46 ESS: {explicit_unsafe_result.stdout.strip()}")
    if "unsafe" in default_result.stderr.lower():
        raise AssertionError("default bounded-error mode unexpectedly printed an unsafe warning")
    assert_unsafe_warning(explicit_unsafe_result, "explicit_unsafe_warning")

    # Numerical mode and parser behavior
    exact_result = assert_success_with_ess_output(
        fracessa_exe, ["--exact", "2#0,1,0"], "exact_success"
    )
    if "unsafe numerical mode" in exact_result.stderr.lower():
        raise AssertionError("exact mode unexpectedly printed the unsafe warning")

    # Affine normalization restores the exact result for both historical cases.
    for case_name, matrix in (
        ("normalized_scale", "2#0,1/100000000000000000000,0"),
        ("normalized_translation", "2#100000000000000000000,100000000000000000001,100000000000000000000"),
    ):
        result = assert_success_with_ess_output(fracessa_exe, [matrix], case_name)
        if first_non_empty_line(result.stdout) != "1":
            raise AssertionError(f"{case_name}: expected one ESS, got {result.stdout.strip()}")

    # Other success paths
    assert_success_with_ess_output(fracessa_exe, ["5#1,3"], "circular_success")
    assert_candidate_header_matches_rows(fracessa_exe)

    assert_failure_with_stderr(
        fracessa_exe,
        ["--exact", "--unsafe", "2#0,1,0"],
        "cannot be used together",
        "exact_unsafe_conflict",
    )

    # Parser failure paths
    assert_failure_with_stderr(
        fracessa_exe,
        ["2,0,1,0"],
        "does not include '#'",
        "missing_hash_rejected",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["2#0#1"],
        "Invalid character",
        "multiple_hash_rejected",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["--unsafe", "64#1"],
        "supports dimensions in [1, 63]",
        "numerical_unsafe_uses_dimension_guard",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["2#1/0,0,1"],
        "denominator cannot be zero",
        "zero_denominator_rejected",
    )

    print("[OK] parser black-box checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
