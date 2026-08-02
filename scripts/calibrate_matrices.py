#!/usr/bin/env python3
"""Fill missing matrix timing calibrations and safe exact baselines."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from multiprocessing import get_context
from multiprocessing.connection import Connection
import os
from pathlib import Path
import sqlite3
from statistics import median
import sys
from time import monotonic


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = REPOSITORY_ROOT / "testdata" / "fracessa_testdata.sqlite3"
DEFAULT_MODULE_DIR = REPOSITORY_ROOT / "cpp" / "build-benchmark"
DEFAULT_CUTOFF_SECONDS = 1.0
STARTUP_TIMEOUT_SECONDS = 30.0


def _compute_matrix(
    module_dir: str,
    cpu_id: int,
    method: str,
    matrix_id: int,
    matrix: str,
    include_baseline: bool,
    single_sample: bool,
    cutoff_ns: int,
    messages: Connection,
) -> None:
    """Run one matrix in a disposable process so the parent can enforce a cutoff."""

    os.sched_setaffinity(0, {cpu_id})
    sys.path.insert(0, module_dir)
    import fracessa_core

    messages.send(("ready", None))
    if include_baseline:
        result = fracessa_core.compute_matrix(
            method=method,
            matrix=matrix,
            include_candidates=True,
            full_support=False,
            enable_logging=False,
            matrix_id=matrix_id,
        )
        messages.send(("baseline", result))
        if result["status"] != 0 or single_sample:
            return

    for sample_index in range(1 if single_sample else 3):
        result = fracessa_core.compute_matrix(
            method=method,
            matrix=matrix,
            include_candidates=False,
            full_support=False,
            enable_logging=False,
            matrix_id=matrix_id,
        )
        messages.send(("sample", result))
        if result["status"] != 0 or (sample_index == 0 and int(result["elapsed_ns"]) >= cutoff_ns):
            return


def _terminate(process) -> None:
    process.terminate()
    process.join(timeout=1)
    if process.is_alive():
        process.kill()
        process.join()


def _receive(messages: Connection, process, timeout_seconds: float, matrix_id: int):
    if not messages.poll(timeout_seconds):
        if process.is_alive():
            return None
        raise RuntimeError(f"matrix {matrix_id}: calibration process exited without a result")
    try:
        return messages.recv()
    except EOFError:
        process.join(timeout=0)
        return "worker_exit", process.exitcode


def _check_result(result: dict, matrix_id: int) -> int:
    if result["status"] != 0 or not result["success"]:
        raise RuntimeError(f"matrix {matrix_id}: {result['error_message'] or 'native computation failed'}")
    elapsed_ns = int(result["elapsed_ns"])
    if elapsed_ns <= 0:
        raise RuntimeError(f"matrix {matrix_id}: native timing must be positive")
    return elapsed_ns


def _exact_baseline(result: dict, matrix_id: int) -> tuple[tuple, list[tuple]]:
    candidates = result["candidates"]
    candidate_structure: Counter[int] = Counter()
    ess_structure: Counter[int] = Counter()
    candidate_rows = []
    weighted_candidates = 0
    weighted_ess = 0

    if [int(row["candidate_id"]) for row in candidates] != list(range(1, len(candidates) + 1)):
        raise RuntimeError(f"matrix {matrix_id}: candidate IDs are not consecutive")
    if len({int(row["support"]) for row in candidates}) != len(candidates):
        raise RuntimeError(f"matrix {matrix_id}: duplicate candidate supports")

    for row in candidates:
        multiplier = row["multiplier"]
        weight = 1 if multiplier is None else int(multiplier)
        if weight < 1:
            raise RuntimeError(f"matrix {matrix_id}: invalid candidate multiplier")
        support_size = int(row["support_size"])
        is_ess = bool(row["is_ess"])
        weighted_candidates += weight
        candidate_structure[support_size] += weight
        if is_ess:
            weighted_ess += weight
            ess_structure[support_size] += weight
        candidate_rows.append(
            (
                matrix_id,
                int(row["candidate_id"]),
                str(row["vector"]),
                int(row["support"]),
                support_size,
                int(row["extended_support"]),
                int(row["extended_support_size"]),
                None if multiplier is None else int(multiplier),
                int(is_ess),
                str(row["stability"]),
                str(row["payoff"]),
                str(float(row["payoff_dbl"])),
            )
        )

    ess_count = int(result["ess_count"])
    if weighted_ess != ess_count:
        raise RuntimeError(f"matrix {matrix_id}: weighted ESS count {weighted_ess} != native count {ess_count}")
    summary = (
        weighted_candidates,
        ess_count,
        json.dumps(dict(sorted(candidate_structure.items())), separators=(",", ":")),
        json.dumps(dict(sorted(ess_structure.items())), separators=(",", ":")),
    )
    return summary, candidate_rows


def _run_one(
    context,
    module_dir: Path,
    cpu_id: int,
    cutoff_seconds: float,
    method: str,
    matrix_id: int,
    matrix: str,
    needs_baseline: bool,
    expected_ess: int | None,
    single_sample: bool,
) -> tuple[int, tuple[tuple, list[tuple]] | None, str | None]:
    receive_messages, send_messages = context.Pipe(duplex=False)
    cutoff_ns = int(cutoff_seconds * 1_000_000_000)
    process = context.Process(
        target=_compute_matrix,
        args=(str(module_dir), cpu_id, method, matrix_id, matrix, needs_baseline, single_sample, cutoff_ns, send_messages),
    )
    process.start()
    send_messages.close()
    samples = []
    baseline = None
    timed_out = False
    failure = None
    try:
        message = _receive(receive_messages, process, STARTUP_TIMEOUT_SECONDS, matrix_id)
        if message is None or message[0] != "ready":
            raise RuntimeError(f"matrix {matrix_id}: calibration process did not become ready")

        if needs_baseline:
            message = _receive(receive_messages, process, cutoff_seconds, matrix_id)
            if message is None:
                timed_out = True
            elif message[0] == "worker_exit":
                failure = f"worker_exit={message[1]}"
            else:
                kind, result = message
                if kind != "baseline":
                    raise RuntimeError(f"matrix {matrix_id}: expected exact baseline result")
                _check_result(result, matrix_id)
                baseline = _exact_baseline(result, matrix_id)
                expected_ess = int(result["ess_count"])
                if single_sample:
                    samples.append(int(result["elapsed_ns"]))

        expected_samples = 1
        while not timed_out and failure is None and len(samples) < expected_samples:
            message = _receive(receive_messages, process, cutoff_seconds, matrix_id)
            if message is None:
                timed_out = True
                break
            if message[0] == "worker_exit":
                failure = f"worker_exit={message[1]}"
                break
            kind, result = message
            if kind != "sample":
                raise RuntimeError(f"matrix {matrix_id}: expected calibration sample")
            elapsed_ns = _check_result(result, matrix_id)
            observed_ess = int(result["ess_count"])
            if expected_ess is not None and observed_ess != expected_ess:
                raise RuntimeError(f"matrix {matrix_id}: ESS count {observed_ess} != expected {expected_ess}")
            samples.append(elapsed_ns)
            if len(samples) == 1 and elapsed_ns < cutoff_ns:
                expected_samples = 3
    finally:
        if process.is_alive():
            if timed_out:
                _terminate(process)
            else:
                process.join(timeout=1)
        if process.is_alive():
            _terminate(process)
        receive_messages.close()

    return (-1 if timed_out or failure else int(median(samples))), baseline, ("timeout" if timed_out else failure)


def _store_result(
    connection: sqlite3.Connection,
    method: str,
    matrix_id: int,
    calibration_ns: int,
    baseline: tuple[tuple, list[tuple]] | None,
    retry_timeout: bool,
) -> None:
    calibration_column = f"{method}_calibration_ns"
    with connection:
        if baseline is not None:
            summary, candidate_rows = baseline
            existing_rows = connection.execute("SELECT COUNT(*) FROM candidates WHERE matrix_id = ?", (matrix_id,)).fetchone()[0]
            if existing_rows:
                raise RuntimeError(f"matrix {matrix_id}: cannot add a baseline over existing candidate rows")
            connection.executemany(
                """INSERT INTO candidates (
                       matrix_id, candidate_id, vector, support, support_size, extended_support,
                       extended_support_size, multiplier, is_ess, reason_ess, payoff, payoff_double
                   ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                candidate_rows,
            )
            calibration_condition = f"{calibration_column} = -1" if retry_timeout else f"{calibration_column} IS NULL"
            cursor = connection.execute(
                f"""UPDATE matrices
                    SET candidate_count = ?, ess_count = ?, candidate_structure = ?, ess_structure = ?, {calibration_column} = ?
                    WHERE matrix_id = ? AND candidate_count IS NULL AND {calibration_condition}""",
                (*summary, calibration_ns, matrix_id),
            )
        else:
            calibration_condition = f"{calibration_column} = -1" if retry_timeout else f"{calibration_column} IS NULL"
            cursor = connection.execute(
                f"UPDATE matrices SET {calibration_column} = ? WHERE matrix_id = ? AND {calibration_condition}",
                (calibration_ns, matrix_id),
            )
        if cursor.rowcount != 1:
            raise RuntimeError(f"matrix {matrix_id}: database row changed during calibration")


def calibrate(
    method: str,
    module_dir: Path,
    database: Path,
    cpu_id: int,
    cutoff_seconds: float,
    retry_timeouts: bool,
) -> None:
    available_cpus = set(os.sched_getaffinity(0))
    if cpu_id not in available_cpus:
        raise ValueError(f"CPU {cpu_id} is unavailable; choose one of {sorted(available_cpus)}")
    if cutoff_seconds <= 0:
        raise ValueError("cutoff must be positive")
    if retry_timeouts and method != "safe":
        raise ValueError("timeout retries are supported only for safe calibration")
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    if not list(module_dir.glob("fracessa_core*.so")):
        raise ValueError(f"fracessa_core is missing from {module_dir}")

    calibration_column = f"{method}_calibration_ns"
    context = get_context("spawn")
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        inconsistent = connection.execute(
            """SELECT m.matrix_id
               FROM matrices m
               WHERE m.candidate_count IS NULL
                 AND EXISTS (SELECT 1 FROM candidates c WHERE c.matrix_id = m.matrix_id)
               LIMIT 1"""
        ).fetchone()
        if inconsistent:
            raise RuntimeError(f"matrix {inconsistent[0]} has candidate rows but no baseline summary")
        calibration_condition = f"{calibration_column} = -1" if retry_timeouts else f"{calibration_column} IS NULL"
        rows = connection.execute(
            f"""SELECT matrix_id, dimension, matrix, candidate_count, ess_count
                FROM matrices
                WHERE {calibration_condition}
                ORDER BY dimension, matrix_id"""
        ).fetchall()
        started = monotonic()
        for position, (matrix_id, dimension, values, candidate_count, ess_count) in enumerate(rows, start=1):
            needs_baseline = method == "safe" and candidate_count is None
            calibration_ns, baseline, failure = _run_one(
                context,
                module_dir,
                cpu_id,
                cutoff_seconds,
                method,
                matrix_id,
                f"{dimension}#{values}",
                needs_baseline,
                ess_count,
                retry_timeouts,
            )
            _store_result(connection, method, matrix_id, calibration_ns, baseline, retry_timeouts)
            calibration = failure or ("timeout" if calibration_ns == -1 else f"{calibration_ns / 1_000:.3f} us")
            baseline_status = " baseline=filled" if baseline is not None else ""
            print(
                f"[{position}/{len(rows)}] matrix={matrix_id} dimension={dimension} calibration={calibration}{baseline_status} "
                f"total_minutes={(monotonic() - started) / 60:.2f}",
                flush=True,
            )
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=("fast", "safe"))
    parser.add_argument("--module-dir", type=Path, default=DEFAULT_MODULE_DIR)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--cutoff-seconds", type=float, default=DEFAULT_CUTOFF_SECONDS)
    parser.add_argument("--retry-timeouts", action="store_true")
    arguments = parser.parse_args()
    calibrate(
        arguments.method,
        arguments.module_dir.resolve(),
        arguments.database.resolve(),
        arguments.cpu,
        arguments.cutoff_seconds,
        arguments.retry_timeouts,
    )


if __name__ == "__main__":
    main()
