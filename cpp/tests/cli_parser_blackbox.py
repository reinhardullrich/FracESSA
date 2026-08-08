#!/usr/bin/env python3
"""
Black-box parser checks against the built `fracessa` executable.

These tests validate CLI parser behavior through the real binary, complementing
unit tests on matrix_parser internals.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def run_case(fracessa_exe: Path, args: list[str], timeout: float = 30.0) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(fracessa_exe), *args],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


SUMMARY_FIELDS = [
    "run_id",
    "matrix_id",
    "status",
    "candidate_count",
    "ess_count",
    "candidate_structure",
    "ess_structure",
    "elapsed_ns",
    "safe_fallback",
    "error_message",
]


def output_lines(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if line.strip()]


def parse_summary(text: str, case_name: str) -> dict:
    lines = output_lines(text)
    if not lines:
        raise AssertionError(f"{case_name}: missing JSON summary")
    try:
        summary = json.loads(lines[0])
    except json.JSONDecodeError as error:
        raise AssertionError(f"{case_name}: invalid JSON summary: {lines[0]!r}") from error
    if list(summary) != SUMMARY_FIELDS:
        raise AssertionError(f"{case_name}: wrong summary fields: {list(summary)}")
    return summary


def comparable_output(result: subprocess.CompletedProcess, case_name: str) -> tuple[dict, list[str]]:
    summary = parse_summary(result.stdout, case_name)
    comparable = {
        key: summary[key]
        for key in (
            "matrix_id",
            "status",
            "candidate_count",
            "ess_count",
            "candidate_structure",
            "ess_structure",
            "error_message",
        )
    }
    return comparable, output_lines(result.stdout)[1:]


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

    summary = parse_summary(result.stdout, case_name)
    if summary["status"] != 0:
        raise AssertionError(f"{case_name}: expected status 0, got {summary['status']}")
    for field in ("candidate_count", "ess_count", "elapsed_ns"):
        if type(summary[field]) is not int or summary[field] < 0:
            raise AssertionError(f"{case_name}: invalid {field}: {summary[field]!r}")
    if summary["run_id"] is not None or summary["error_message"]:
        raise AssertionError(f"{case_name}: invalid success summary: {summary}")
    if "--candidates" not in args and len(output_lines(result.stdout)) != 1:
        raise AssertionError(f"{case_name}: expected exactly one output line")
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
    summary = parse_summary(result.stdout, case_name)
    if summary["status"] == 0 or not summary["error_message"]:
        raise AssertionError(f"{case_name}: error summary does not describe the failure: {summary}")
    if len(output_lines(result.stdout)) != 1:
        raise AssertionError(f"{case_name}: expected exactly one error-summary line")


def assert_candidate_header_matches_rows(fracessa_exe: Path) -> None:
    result = run_case(fracessa_exe, ["--candidates", "safe", "1#1/3"])
    if result.returncode != 0:
        raise AssertionError(
            f"candidate_columns: expected success, got rc={result.returncode}, stderr={result.stderr.strip()}"
        )

    lines = output_lines(result.stdout)
    if len(lines) < 3:
        raise AssertionError(f"candidate_columns: expected summary, header, and candidate row; got {lines}")
    summary = parse_summary(result.stdout, "candidate_columns")
    if summary["candidate_count"] != 1 or summary["candidate_structure"] != {"1": 1}:
        raise AssertionError(f"candidate_columns: wrong candidate summary: {summary}")

    header_count = len(lines[1].split(";"))
    for row in lines[2:]:
        row_count = len(row.split(";"))
        if row_count != header_count:
            raise AssertionError(
                f"candidate_columns: header has {header_count} fields, candidate row has {row_count}"
            )
    payoff_double = lines[2].rsplit(";", 1)[1]
    if payoff_double != "0.33333333333333331":
        raise AssertionError(f"candidate_columns: expected round-trip binary64 output, got {payoff_double}")


def assert_logged_candidates_match_sorted_output(fracessa_exe: Path) -> None:
    matrix = "13#7,18,18,10,7,10"
    with tempfile.TemporaryDirectory() as directory:
        result = subprocess.run(
            [str(fracessa_exe.resolve()), "--candidates", "--log", "safe", matrix],
            capture_output=True,
            text=True,
            timeout=30.0,
            cwd=directory,
        )
        if result.returncode != 0:
            raise AssertionError(
                f"logged_candidate_ids: expected success, got rc={result.returncode}, stderr={result.stderr.strip()}"
            )

        lines = output_lines(result.stdout)
        summary = parse_summary(result.stdout, "logged_candidate_ids")
        if summary["ess_count"] != 143 or len(lines[2:]) != 7:
            raise AssertionError("logged_candidate_ids: wrong ESS total or representative count")

        expected_rows = lines[1:]
        expected_ids = list(range(1, len(expected_rows)))
        actual_ids = [int(row.split(";", 1)[0]) for row in expected_rows[1:]]
        if actual_ids != expected_ids:
            raise AssertionError(f"logged_candidate_ids: output IDs are not sequential: {actual_ids}")

        log_path = Path(directory) / "log" / "fracessa.log"
        marker = "] [info] "
        log_payloads = [
            line.split(marker, 1)[1]
            for line in log_path.read_text(encoding="utf-8").splitlines()
            if marker in line
        ]
        header_positions = [index for index, payload in enumerate(log_payloads) if payload == expected_rows[0]]
        if not header_positions:
            raise AssertionError("logged_candidate_ids: candidate header missing from log")
        start = header_positions[-1]
        logged_rows = log_payloads[start:start + len(expected_rows)]
        if logged_rows != expected_rows:
            raise AssertionError("logged_candidate_ids: logged rows differ from sorted CLI candidate output")


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

    for removed_or_unknown_method in ("verified", "exact", "unsafe", "unknown", 'bad"method'):
        assert_failure_with_stderr(
            fracessa_exe,
            [removed_or_unknown_method, "2#0,1,0"],
            "Unknown search method",
            f"{removed_or_unknown_method}_method_rejected",
        )
    assert_failure_with_stderr(fracessa_exe, ["--mode", "safe", "2#0,1,0"], "Unknown argument: --mode", "mode_flag_removed")
    assert_failure_with_stderr(fracessa_exe, ["--timing", "safe", "2#0,1,0"], "Unknown argument: --timing", "timing_flag_removed")
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
        if parse_summary(safe_result.stdout, case_name)["ess_count"] != 1:
            raise AssertionError(f"{case_name}: safe search missed the ESS")
        if parse_summary(fast_result.stdout, case_name)["ess_count"] != 1:
            raise AssertionError(f"{case_name}: fast search did not fall back to safe analysis")
        if "warning" in fast_result.stderr.lower():
            raise AssertionError(f"{case_name}: fast search unexpectedly printed a warning")

    fallback_result = assert_success_with_ess_output(
        fracessa_exe, ["fast", "2#1,1000000000,1"], "safe_fallback_output"
    )
    fallback_summary = parse_summary(fallback_result.stdout, "safe_fallback_output")
    if fallback_summary["safe_fallback"] != "precision_span":
        raise AssertionError(f"safe_fallback_output: expected precision_span, got {fallback_summary}")

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
        if parse_summary(safe_result.stdout, case_name)["ess_count"] != int(expected_ess_count):
            raise AssertionError(f"{case_name}: safe search returned the wrong ESS count")

        for method in ("fast", "test"):
            result = assert_success_with_ess_output(fracessa_exe, [*options, method, matrix], f"{case_name}_{method}")
            if comparable_output(result, case_name) != comparable_output(safe_result, case_name):
                raise AssertionError(f"{case_name}: {method} search differs from safe search")

    # Other success paths
    assert_success_with_ess_output(fracessa_exe, ["safe", "5#1,3"], "circular_success")
    assert_logged_candidates_match_sorted_output(fracessa_exe)
    assert_failure_with_stderr(
        fracessa_exe,
        ["--cyclic-symmetry-filter", "safe", "5#1,3"],
        "Unknown argument: --cyclic-symmetry-filter",
        "cyclic_filter_flag_removed",
    )
    matrix_id_result = assert_success_with_ess_output(
        fracessa_exe,
        ["--matrixid", "9223372036854775807", "safe", "2#0,1,0"],
        "signed_64_bit_matrix_id",
    )
    if parse_summary(matrix_id_result.stdout, "signed_64_bit_matrix_id")["matrix_id"] != 9223372036854775807:
        raise AssertionError("signed_64_bit_matrix_id: summary lost the matrix ID")
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
        ["fast", "65#1"],
        "supports dimensions in [1, 64]",
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
