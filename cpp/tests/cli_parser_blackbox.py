#!/usr/bin/env python3
"""
Black-box parser checks against the built `fracessa` executable.

These tests validate CLI parser behavior through the real binary, complementing
unit tests on matrix_parser internals.
"""

from __future__ import annotations

import argparse
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
    result = run_case(fracessa_exe, ["--candidates", "safe", "2#0,1,0"])
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

    # Search-method and parser behavior
    safe_result = assert_success_with_ess_output(
        fracessa_exe, ["safe", "2#0,1,0"], "safe_success"
    )
    if "warning" in safe_result.stderr.lower():
        raise AssertionError("safe search unexpectedly printed a warning")
    assert_success_with_ess_output(fracessa_exe, ["test", "2#0,1,0"], "test_success")

    for removed_or_unknown_method in ("verified", "exact", "unsafe", "unknown"):
        assert_failure_with_stderr(
            fracessa_exe,
            [removed_or_unknown_method, "2#0,1,0"],
            "Unknown search method",
            f"{removed_or_unknown_method}_method_rejected",
        )
    assert_failure_with_stderr(fracessa_exe, ["--mode", "safe", "2#0,1,0"], "Unknown argument: --mode", "mode_flag_removed")
    assert_failure_with_stderr(fracessa_exe, ["2#0,1,0"], "matrix", "missing_method_rejected")

    # Fast search bypasses floating-point rejection for a game whose exact-to-double conversion triggers any input check.
    for case_name, matrix in (
        ("small_difference", "2#0,1/100000000000000000000,0"),
        ("collapsed_values", "2#100000000000000000000,100000000000000000001,100000000000000000000"),
    ):
        safe_result = assert_success_with_ess_output(
            fracessa_exe, ["safe", matrix], f"{case_name}_safe"
        )
        fast_result = assert_success_with_ess_output(
            fracessa_exe, ["fast", matrix], f"{case_name}_fast"
        )
        if first_non_empty_line(safe_result.stdout) != "1":
            raise AssertionError(f"{case_name}: safe search missed the ESS")
        if first_non_empty_line(fast_result.stdout) != "1":
            raise AssertionError(f"{case_name}: fast search did not fall back to safe analysis")
        if "warning" in fast_result.stderr.lower():
            raise AssertionError(f"{case_name}: fast search unexpectedly printed a warning")

    timing_result = assert_success_with_ess_output(
        fracessa_exe, ["--timing", "fast", "2#1,1000000000,1"], "safe_fallback_output"
    )
    timing_lines = [line.strip() for line in timing_result.stdout.splitlines() if line.strip()]
    if len(timing_lines) != 3 or timing_lines[2] != "precision_span":
        raise AssertionError(f"safe_fallback_output: expected precision_span on line 3, got {timing_lines}")

    # Fast and test must preserve all three exact ESS that their historical floating-point rejection tests lost.
    historical_false_rejections = (
        (
            "small_pivot",
            ["--fullsupport", "--candidates"],
            "1",
            "3#-3,1,2,-1000000000001/3000000000000,-1999999999999/3000000000000,-4000000000001/3000000000000",
        ),
        (
            "positive_probability",
            ["--fullsupport", "--candidates"],
            "1",
            "3#1/50000000,4001/200000000000,5000000001/100000000,1/50000000,10000000002001/200000000000,0",
        ),
        (
            "outside_payoff",
            ["--candidates"],
            "2",
            "4#1/5,40000000001/200000000000,501/10,6040000000001/400000000000,1/5,10020000000001/200000000000,"
            "10048000000001/400000000000,0,10040000000001/400000000000,0",
        ),
    )
    for case_name, options, expected_ess_count, matrix in historical_false_rejections:
        safe_result = assert_success_with_ess_output(fracessa_exe, [*options, "safe", matrix], f"{case_name}_safe")
        if first_non_empty_line(safe_result.stdout) != expected_ess_count:
            raise AssertionError(f"{case_name}: safe search returned the wrong ESS count")

        for method in ("fast", "test"):
            result = assert_success_with_ess_output(fracessa_exe, [*options, method, matrix], f"{case_name}_{method}")
            if result.stdout != safe_result.stdout:
                raise AssertionError(f"{case_name}: {method} search differs from safe search")

    # Other success paths
    assert_success_with_ess_output(fracessa_exe, ["safe", "5#1,3"], "circular_success")
    cyclic_off = assert_success_with_ess_output(fracessa_exe, ["safe", "8#7,11,7,13"], "cyclic_filter_off")
    cyclic_on = assert_success_with_ess_output(
        fracessa_exe,
        ["--cyclic-symmetry-filter", "safe", "8#7,11,7,13"],
        "cyclic_filter_on",
    )
    if first_non_empty_line(cyclic_on.stdout) != first_non_empty_line(cyclic_off.stdout):
        raise AssertionError("cyclic symmetry filter changed the represented ESS count")
    assert_success_with_ess_output(
        fracessa_exe,
        ["--matrixid", "9223372036854775807", "safe", "2#0,1,0"],
        "signed_64_bit_matrix_id",
    )
    assert_candidate_header_matches_rows(fracessa_exe)

    # Parser failure paths
    assert_failure_with_stderr(
        fracessa_exe,
        ["safe", "2,0,1,0"],
        "does not include '#'",
        "missing_hash_rejected",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["safe", "2#0#1"],
        "Multiple '#'",
        "multiple_hash_rejected",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["fast", "64#1"],
        "supports dimensions in [1, 63]",
        "fast_uses_dimension_guard",
    )
    assert_failure_with_stderr(
        fracessa_exe,
        ["safe", "2#1/0,0,1"],
        "denominator cannot be zero",
        "zero_denominator_rejected",
    )

    print("[OK] parser black-box checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
