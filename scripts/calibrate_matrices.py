#!/usr/bin/env python3
"""Fill matrix timing calibrations and candidate results."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
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
DEFAULT_AUDIT_CUTOFF_SECONDS = 120.0
AUDIT_TARGET_NS = 1_000_000_000
STARTUP_TIMEOUT_SECONDS = 30.0
SAFE_FALLBACKS = {"precision_span", "equilibration_invalid", "equilibration_non_convergence"}


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


def _compute_audit_matrix(
    module_dir: str,
    cpu_id: int,
    method: str,
    matrix_id: int,
    matrix: str,
    iterations: int,
    messages: Connection,
) -> None:
    """Measure one matrix repeatedly and return complete candidates with the first sample."""

    os.sched_setaffinity(0, {cpu_id})
    sys.path.insert(0, module_dir)
    import fracessa_core

    messages.send(("ready", None))
    samples = []
    baseline_result = None
    safe_fallback = None
    expected_ess = None
    for sample_index in range(iterations):
        result = fracessa_core.compute_matrix(
            method=method,
            matrix=matrix,
            include_candidates=sample_index == 0,
            full_support=False,
            enable_logging=False,
            matrix_id=matrix_id,
        )
        try:
            elapsed_ns = _check_result(result, matrix_id)
            current_fallback = _safe_fallback(result, matrix_id)
            current_ess = int(result["ess_count"])
            if sample_index == 0:
                baseline_result = result
                safe_fallback = current_fallback
                expected_ess = current_ess
            elif current_fallback != safe_fallback:
                raise RuntimeError("safe fallback changed between audit samples")
            elif current_ess != expected_ess:
                raise RuntimeError(f"ESS count changed from {expected_ess} to {current_ess} between audit samples")
            samples.append(elapsed_ns)
        except RuntimeError as error:
            messages.send(("failure", str(error)))
            return
    messages.send(("audit", (baseline_result, int(median(samples)), safe_fallback, len(samples))))


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
    if result["status"] != 0:
        raise RuntimeError(f"matrix {matrix_id}: {result['error_message'] or 'native computation failed'}")
    elapsed_ns = int(result["elapsed_ns"])
    if elapsed_ns <= 0:
        raise RuntimeError(f"matrix {matrix_id}: native timing must be positive")
    return elapsed_ns


def _safe_fallback(result: dict, matrix_id: int) -> str | None:
    fallback = result.get("safe_fallback")
    if fallback is not None and fallback not in SAFE_FALLBACKS:
        raise RuntimeError(f"matrix {matrix_id}: unknown safe fallback {fallback!r}")
    return fallback


def _candidate_result(result: dict, matrix_id: int) -> tuple[tuple, list[tuple]]:
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


def _audit_iterations(calibration_ns: int | None) -> int:
    if calibration_ns is None or calibration_ns <= 0:
        return 1
    return max(1, (AUDIT_TARGET_NS + calibration_ns - 1) // calibration_ns)


def _audit_cutoff_seconds(calibration_ns: int | None, minimum_cutoff_seconds: float) -> float:
    if calibration_ns is None or calibration_ns <= 0:
        return minimum_cutoff_seconds
    return max(minimum_cutoff_seconds, 1.2 * calibration_ns / 1_000_000_000)


def _run_audit_one(
    context,
    module_dir: Path,
    cpu_id: int,
    method: str,
    matrix_id: int,
    matrix: str,
    calibration_ns: int | None,
    minimum_cutoff_seconds: float,
) -> tuple[int, tuple[tuple, list[tuple]] | None, str, str | None, int, int, float]:
    iterations = _audit_iterations(calibration_ns)
    cutoff_seconds = _audit_cutoff_seconds(calibration_ns, minimum_cutoff_seconds)
    receive_messages, send_messages = context.Pipe(duplex=False)
    process = context.Process(
        target=_compute_audit_matrix,
        args=(str(module_dir), cpu_id, method, matrix_id, matrix, iterations, send_messages),
    )
    process.start()
    send_messages.close()
    baseline = None
    safe_fallback = None
    status = "completed"
    completed_samples = 0
    try:
        message = _receive(receive_messages, process, STARTUP_TIMEOUT_SECONDS, matrix_id)
        if message is None or message[0] != "ready":
            status = "startup_timeout"
        else:
            message = _receive(receive_messages, process, cutoff_seconds, matrix_id)
            if message is None:
                status = "timeout"
            elif message[0] == "worker_exit":
                status = f"worker_exit_{message[1]}"
            elif message[0] == "failure":
                status = f"failure: {message[1]}"
            elif message[0] != "audit":
                status = f"unexpected_message_{message[0]}"
            else:
                baseline_result, measured_ns, safe_fallback, completed_samples = message[1]
                baseline = _candidate_result(baseline_result, matrix_id)
    finally:
        if process.is_alive():
            if status != "completed" or completed_samples < iterations:
                _terminate(process)
            else:
                process.join(timeout=1)
        if process.is_alive():
            _terminate(process)
        receive_messages.close()

    if status != "completed" or completed_samples != iterations:
        measured_ns = -1
    return measured_ns, baseline, status, safe_fallback, iterations, completed_samples, cutoff_seconds


def _load_stored_baseline(connection: sqlite3.Connection, matrix_id: int) -> tuple[tuple, list[tuple]] | None:
    summary = connection.execute(
        "SELECT candidate_count, ess_count, candidate_structure, ess_structure FROM matrices WHERE matrix_id = ?",
        (matrix_id,),
    ).fetchone()
    candidate_rows = connection.execute(
        """SELECT matrix_id, candidate_id, vector, support, support_size, extended_support, extended_support_size,
                  multiplier, is_ess, reason_ess, payoff, payoff_double
           FROM candidates WHERE matrix_id = ? ORDER BY candidate_id""",
        (matrix_id,),
    ).fetchall()
    if summary[0] is None:
        if candidate_rows:
            raise RuntimeError(f"matrix {matrix_id}: candidate rows exist without a baseline summary")
        return None
    if not candidate_rows and int(summary[0]) != 0:
        raise RuntimeError(f"matrix {matrix_id}: baseline summary exists without candidate rows")
    return tuple(summary), [tuple(row) for row in candidate_rows]


def _baseline_differences(
    stored: tuple[tuple, list[tuple]], observed: tuple[tuple, list[tuple]], method: str
) -> list[str]:
    differences = []
    stored_summary, stored_candidates = stored
    observed_summary, observed_candidates = observed
    labels = ("candidate_count", "ess_count", "candidate_structure", "ess_structure")
    for index, label in enumerate(labels):
        stored_value = stored_summary[index]
        observed_value = observed_summary[index]
        if label.endswith("_structure"):
            stored_value = json.loads(stored_value)
            observed_value = json.loads(observed_value)
        if stored_value != observed_value:
            differences.append(f"{method} {label} stored={stored_value!r} observed={observed_value!r}")
    if stored_candidates != observed_candidates:
        first_difference = next(
            (
                index + 1
                for index, (stored_row, observed_row) in enumerate(zip(stored_candidates, observed_candidates))
                if stored_row != observed_row
            ),
            min(len(stored_candidates), len(observed_candidates)) + 1,
        )
        differences.append(
            f"{method} candidates differ at candidate_id={first_difference} "
            f"(stored_rows={len(stored_candidates)}, observed_rows={len(observed_candidates)})"
        )
    return differences


def _insert_baseline(
    connection: sqlite3.Connection, matrix_id: int, baseline: tuple[tuple, list[tuple]]
) -> None:
    summary, candidate_rows = baseline
    connection.executemany(
        """INSERT INTO candidates (
               matrix_id, candidate_id, vector, support, support_size, extended_support,
               extended_support_size, multiplier, is_ess, reason_ess, payoff, payoff_double
           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        candidate_rows,
    )
    connection.execute(
        """UPDATE matrices
           SET candidate_count = ?, ess_count = ?, candidate_structure = ?, ess_structure = ?
           WHERE matrix_id = ? AND candidate_count IS NULL""",
        (*summary, matrix_id),
    )


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
) -> tuple[int, tuple[tuple, list[tuple]] | None, str | None, str | None]:
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
    safe_fallback = None
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
                    raise RuntimeError(f"matrix {matrix_id}: expected candidate result")
                _check_result(result, matrix_id)
                safe_fallback = _safe_fallback(result, matrix_id)
                baseline = _candidate_result(result, matrix_id)
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
            current_safe_fallback = _safe_fallback(result, matrix_id)
            if samples and current_safe_fallback != safe_fallback:
                raise RuntimeError(f"matrix {matrix_id}: safe fallback changed between calibration samples")
            safe_fallback = current_safe_fallback
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

    return (-1 if timed_out or failure else int(median(samples))), baseline, ("timeout" if timed_out else failure), safe_fallback


def _store_result(
    connection: sqlite3.Connection,
    method: str,
    matrix_id: int,
    calibration_ns: int,
    baseline: tuple[tuple, list[tuple]] | None,
    retry_timeout: bool,
    safe_fallback: str | None,
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
                    SET candidate_count = ?, ess_count = ?, candidate_structure = ?, ess_structure = ?, {calibration_column} = ?,
                        safe_fallback = ?
                    WHERE matrix_id = ? AND candidate_count IS NULL AND {calibration_condition}""",
                (*summary, calibration_ns, safe_fallback, matrix_id),
            )
        else:
            calibration_condition = f"{calibration_column} = -1" if retry_timeout else f"{calibration_column} IS NULL"
            cursor = connection.execute(
                f"UPDATE matrices SET {calibration_column} = ?, safe_fallback = ? WHERE matrix_id = ? AND {calibration_condition}",
                (calibration_ns, safe_fallback, matrix_id),
            )
        if cursor.rowcount != 1:
            raise RuntimeError(f"matrix {matrix_id}: database row changed during calibration")


def calibrate(
    method: str,
    module_dir: Path,
    database: Path,
    cpu_ids: list[int],
    cutoff_seconds: float,
    retry_timeouts: bool,
) -> None:
    available_cpus = set(os.sched_getaffinity(0))
    if not cpu_ids or len(cpu_ids) != len(set(cpu_ids)):
        raise ValueError("provide at least one CPU without duplicates")
    unavailable_cpus = sorted(set(cpu_ids) - available_cpus)
    if unavailable_cpus:
        raise ValueError(f"CPUs {unavailable_cpus} are unavailable; choose from {sorted(available_cpus)}")
    if cutoff_seconds <= 0:
        raise ValueError("cutoff must be positive")
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    if not list(module_dir.glob("fracessa_core*.so")):
        raise ValueError(f"fracessa_core is missing from {module_dir}")

    sys.path.insert(0, str(module_dir))
    import fracessa_core

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
        row_iterator = iter(rows)
        worker_count = min(len(cpu_ids), len(rows))
        if worker_count == 0:
            return

        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            in_flight = {}

            def submit(row, cpu_id: int) -> None:
                matrix_id, dimension, values, candidate_count, ess_count = row
                matrix = f"{dimension}#{values}"
                classified_fallback = fracessa_core.classify_safe_fallback(matrix)
                if classified_fallback is not None and classified_fallback not in SAFE_FALLBACKS:
                    raise RuntimeError(f"matrix {matrix_id}: unknown classified safe fallback {classified_fallback!r}")
                future = executor.submit(
                    _run_one,
                    context,
                    module_dir,
                    cpu_id,
                    cutoff_seconds,
                    method,
                    matrix_id,
                    matrix,
                    candidate_count is None,
                    ess_count,
                    retry_timeouts,
                )
                in_flight[future] = (cpu_id, matrix_id, dimension, classified_fallback)

            for cpu_id in cpu_ids[:worker_count]:
                submit(next(row_iterator), cpu_id)

            position = 0
            while in_flight:
                completed, _ = wait(in_flight, return_when=FIRST_COMPLETED)
                for future in completed:
                    cpu_id, matrix_id, dimension, classified_fallback = in_flight.pop(future)
                    calibration_ns, baseline, failure, observed_fallback = future.result()
                    if method == "fast" and observed_fallback is not None and observed_fallback != classified_fallback:
                        raise RuntimeError(
                            f"matrix {matrix_id}: computed safe fallback {observed_fallback!r} != "
                            f"classified {classified_fallback!r}"
                        )
                    _store_result(connection, method, matrix_id, calibration_ns, baseline, retry_timeouts, classified_fallback)
                    position += 1
                    calibration = failure or ("timeout" if calibration_ns == -1 else f"{calibration_ns / 1_000:.3f} us")
                    baseline_status = " baseline=filled" if baseline is not None else ""
                    print(
                        f"[{position}/{len(rows)}] matrix={matrix_id} dimension={dimension} cpu={cpu_id} "
                        f"calibration={calibration}{baseline_status} safe_fallback={classified_fallback or 'null'} "
                        f"total_minutes={(monotonic() - started) / 60:.2f}",
                        flush=True,
                    )

                    next_row = next(row_iterator, None)
                    if next_row is not None:
                        submit(next_row, cpu_id)
    finally:
        connection.close()


def _audit_matrix(
    context,
    module_dir: Path,
    cpu_id: int,
    minimum_cutoff_seconds: float,
    row: tuple,
    stored_baseline: tuple[tuple, list[tuple]] | None,
    classified_fallback: str | None,
    initial_messages: list[str],
) -> tuple:
    (
        matrix_id,
        dimension,
        values,
        old_fast_calibration_ns,
        old_safe_calibration_ns,
        _,
        calibration_comment,
    ) = row
    matrix = f"{dimension}#{values}"
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    messages = list(initial_messages)

    try:
        fast_result = _run_audit_one(
            context,
            module_dir,
            cpu_id,
            "fast",
            matrix_id,
            matrix,
            old_fast_calibration_ns,
            minimum_cutoff_seconds,
        )
    except Exception as error:
        fast_iterations = _audit_iterations(old_fast_calibration_ns)
        fast_cutoff_seconds = _audit_cutoff_seconds(old_fast_calibration_ns, minimum_cutoff_seconds)
        fast_result = (-1, None, f"failure: {error}", None, fast_iterations, 0, fast_cutoff_seconds)

    (
        fast_calibration_ns,
        fast_baseline,
        fast_status,
        observed_fallback,
        fast_iterations,
        fast_completed_samples,
        fast_cutoff_seconds,
    ) = fast_result
    if observed_fallback != classified_fallback and fast_baseline is not None:
        messages.append(f"fast safe_fallback classified={classified_fallback!r} observed={observed_fallback!r}")
    safe_fallback = observed_fallback if fast_baseline is not None else classified_fallback
    if fast_baseline is not None and stored_baseline is not None:
        messages.extend(_baseline_differences(stored_baseline, fast_baseline, "fast"))

    event = {
        "timestamp": timestamp,
        "cpu_id": cpu_id,
        "target_seconds": AUDIT_TARGET_NS / 1_000_000_000,
        "fast": {
            "status": fast_status,
            "cutoff_seconds": fast_cutoff_seconds,
            "requested_iterations": fast_iterations,
            "completed_iterations": fast_completed_samples,
            "calibration_ns": None if fast_calibration_ns == -1 else fast_calibration_ns,
            "safe_fallback": safe_fallback,
        },
    }

    safe_baseline = None
    baseline_source = None
    baseline_to_store = None
    if fast_calibration_ns == -1:
        safe_calibration_ns = -1
        event["safe"] = {
            "status": "skipped_fast_timeout_or_failure",
            "cutoff_seconds": None,
            "requested_iterations": 0,
            "completed_iterations": 0,
            "calibration_ns": None,
        }
        messages.append(
            f"fast {fast_status} after cutoff {fast_cutoff_seconds:g} seconds; safe skipped and both calibrations set to -1"
        )
    elif safe_fallback is not None:
        safe_calibration_ns = fast_calibration_ns
        safe_baseline = fast_baseline
        event["safe"] = {
            "status": "copied_from_fast_safe_fallback",
            "cutoff_seconds": fast_cutoff_seconds,
            "requested_iterations": fast_iterations,
            "completed_iterations": fast_completed_samples,
            "calibration_ns": safe_calibration_ns,
        }
    else:
        try:
            safe_result = _run_audit_one(
                context,
                module_dir,
                cpu_id,
                "safe",
                matrix_id,
                matrix,
                old_safe_calibration_ns,
                minimum_cutoff_seconds,
            )
        except Exception as error:
            safe_iterations = _audit_iterations(old_safe_calibration_ns)
            safe_cutoff_seconds = _audit_cutoff_seconds(old_safe_calibration_ns, minimum_cutoff_seconds)
            safe_result = (-1, None, f"failure: {error}", None, safe_iterations, 0, safe_cutoff_seconds)
        (
            safe_calibration_ns,
            safe_baseline,
            safe_status,
            _,
            safe_iterations,
            safe_completed_samples,
            safe_cutoff_seconds,
        ) = safe_result
        event["safe"] = {
            "status": safe_status,
            "cutoff_seconds": safe_cutoff_seconds,
            "requested_iterations": safe_iterations,
            "completed_iterations": safe_completed_samples,
            "calibration_ns": None if safe_calibration_ns == -1 else safe_calibration_ns,
        }
        if safe_calibration_ns == -1:
            messages.append(f"safe {safe_status} after cutoff {safe_cutoff_seconds:g} seconds; calibration set to -1")
        if safe_baseline is not None and stored_baseline is not None:
            messages.extend(_baseline_differences(stored_baseline, safe_baseline, "safe"))
        if fast_baseline is not None and safe_baseline is not None:
            messages.extend(_baseline_differences(safe_baseline, fast_baseline, "fast versus safe"))

    if stored_baseline is None:
        if safe_baseline is not None:
            baseline_to_store = safe_baseline
            baseline_source = "safe" if safe_fallback is None else "fast_safe_fallback"
        else:
            baseline_to_store = fast_baseline
            baseline_source = (
                "fast_safe_fallback" if fast_baseline is not None and safe_fallback is not None
                else "fast_unverified" if fast_baseline is not None
                else None
            )
        if baseline_to_store is not None:
            messages.append(f"filled missing candidate data from {baseline_source}")

    event["messages"] = messages
    try:
        history = json.loads(calibration_comment)
    except (TypeError, json.JSONDecodeError) as error:
        raise RuntimeError(f"matrix {matrix_id}: invalid calibration_comment JSON: {error}") from error
    if not isinstance(history, list):
        raise RuntimeError(f"matrix {matrix_id}: calibration_comment must be a JSON array")
    history.append(event)
    return (
        fast_calibration_ns,
        safe_calibration_ns,
        safe_fallback,
        timestamp,
        json.dumps(history, indent=2),
        baseline_to_store,
        event,
        messages,
    )


def audit(
    module_dir: Path,
    database: Path,
    cpu_ids: list[int],
    minimum_cutoff_seconds: float,
    refresh_all: bool,
) -> None:
    available_cpus = set(os.sched_getaffinity(0))
    if not cpu_ids or len(cpu_ids) != len(set(cpu_ids)):
        raise ValueError("provide at least one CPU without duplicates")
    unavailable_cpus = sorted(set(cpu_ids) - available_cpus)
    if unavailable_cpus:
        raise ValueError(f"CPUs {unavailable_cpus} are unavailable; choose from {sorted(available_cpus)}")
    if minimum_cutoff_seconds <= 0:
        raise ValueError("cutoff must be positive")
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    if not list(module_dir.glob("fracessa_core*.so")):
        raise ValueError(f"fracessa_core is missing from {module_dir}")

    sys.path.insert(0, str(module_dir))
    import fracessa_core

    context = get_context("spawn")
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        condition = "1" if refresh_all else "calibration_timestamp IS NULL"
        rows = connection.execute(
            f"""SELECT matrix_id, dimension, matrix, fast_calibration_ns, safe_calibration_ns,
                       safe_fallback, calibration_comment
                FROM matrices
                WHERE {condition}
                ORDER BY dimension, matrix_id"""
        ).fetchall()
        started = monotonic()
        row_iterator = iter(enumerate(rows, start=1))
        worker_count = min(len(cpu_ids), len(rows))
        if worker_count == 0:
            return

        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            in_flight = {}

            def submit(position_and_row: tuple[int, tuple], cpu_id: int) -> None:
                ordered_position, row = position_and_row
                matrix_id, dimension, values, _, _, stored_fallback, _ = row
                matrix = f"{dimension}#{values}"
                stored_baseline = _load_stored_baseline(connection, matrix_id)
                initial_messages = []
                try:
                    classified_fallback = fracessa_core.classify_safe_fallback(matrix)
                    if classified_fallback is not None and classified_fallback not in SAFE_FALLBACKS:
                        raise RuntimeError(f"unknown classified safe fallback {classified_fallback!r}")
                except Exception as error:
                    classified_fallback = stored_fallback
                    initial_messages.append(f"safe fallback classification failed: {error}")
                future = executor.submit(
                    _audit_matrix,
                    context,
                    module_dir,
                    cpu_id,
                    minimum_cutoff_seconds,
                    row,
                    stored_baseline,
                    classified_fallback,
                    initial_messages,
                )
                in_flight[future] = (cpu_id, ordered_position, matrix_id, dimension, stored_baseline)

            for cpu_id in cpu_ids[:worker_count]:
                submit(next(row_iterator), cpu_id)

            completed_count = 0
            while in_flight:
                completed, _ = wait(in_flight, return_when=FIRST_COMPLETED)
                for future in completed:
                    cpu_id, ordered_position, matrix_id, dimension, stored_baseline = in_flight.pop(future)
                    (
                        fast_calibration_ns,
                        safe_calibration_ns,
                        safe_fallback,
                        timestamp,
                        calibration_comment,
                        baseline_to_store,
                        event,
                        messages,
                    ) = future.result()
                    with connection:
                        if stored_baseline is None and baseline_to_store is not None:
                            _insert_baseline(connection, matrix_id, baseline_to_store)
                        connection.execute(
                            """UPDATE matrices
                               SET fast_calibration_ns = ?, safe_calibration_ns = ?, safe_fallback = ?,
                                   calibration_timestamp = ?, calibration_comment = ?
                               WHERE matrix_id = ?""",
                            (
                                fast_calibration_ns,
                                safe_calibration_ns,
                                safe_fallback,
                                timestamp,
                                calibration_comment,
                                matrix_id,
                            ),
                        )
                    completed_count += 1
                    print(
                        f"[{completed_count}/{len(rows)} order={ordered_position}] matrix={matrix_id} dimension={dimension} "
                        f"cpu={cpu_id} fast={event['fast']['status']}:{fast_calibration_ns} "
                        f"safe={event['safe']['status']}:{safe_calibration_ns} "
                        f"messages={len(messages)} total_minutes={(monotonic() - started) / 60:.2f}",
                        flush=True,
                    )

                    next_row = next(row_iterator, None)
                    if next_row is not None:
                        submit(next_row, cpu_id)
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=("fast", "safe", "audit"))
    parser.add_argument("--module-dir", type=Path, default=DEFAULT_MODULE_DIR)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--cpu", type=int, action="append", metavar="ID", help="CPU ID; repeat to run matrices in parallel")
    parser.add_argument("--cutoff-seconds", type=float)
    parser.add_argument("--retry-timeouts", action="store_true")
    parser.add_argument("--refresh-all", action="store_true", help="audit every matrix even if it already has an audit timestamp")
    arguments = parser.parse_args()
    cpu_ids = arguments.cpu or [2]
    if arguments.method == "audit":
        if arguments.retry_timeouts:
            parser.error("--retry-timeouts does not apply to audit")
        audit(
            arguments.module_dir.resolve(),
            arguments.database.resolve(),
            cpu_ids,
            DEFAULT_AUDIT_CUTOFF_SECONDS if arguments.cutoff_seconds is None else arguments.cutoff_seconds,
            arguments.refresh_all,
        )
        return
    if arguments.refresh_all:
        parser.error("--refresh-all applies only to audit")
    calibrate(
        arguments.method,
        arguments.module_dir.resolve(),
        arguments.database.resolve(),
        cpu_ids,
        DEFAULT_CUTOFF_SECONDS if arguments.cutoff_seconds is None else arguments.cutoff_seconds,
        arguments.retry_timeouts,
    )


if __name__ == "__main__":
    main()
