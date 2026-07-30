"""Record sequential FracESSA timings in the canonical SQLite database.

The runner measures one build per invocation and one matrix at a time. Reuse a
session name across invocations to compare builds without loading incompatible
native extensions into the same Python interpreter.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
import hashlib
from importlib import import_module, invalidate_caches
import math
import os
from pathlib import Path
import platform
import sqlite3
from statistics import median
import subprocess
import sys
from time import perf_counter_ns


DEFAULT_DATABASE = (
    Path(__file__).resolve().parents[2] / "testdata" / "fracessa_testdata.sqlite3"
)
_CLI_TO_NS = {
    "ns": Decimal(1),
    "us": Decimal(1_000),
    "ms": Decimal(1_000_000),
    "s": Decimal(1_000_000_000),
}


def _sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _machine_name() -> str:
    """Return a compact machine and platform identifier."""

    return f"{platform.node()} ({platform.system()} {platform.machine()})"


def _new_session() -> str:
    """Return a sortable UTC benchmark-session identifier."""

    stamp = datetime.now(tz=timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"timing_{stamp}"


def _pin_cpu(cpu_id: int) -> set[int]:
    """Pin this process to ``cpu_id`` and return its previous CPU affinity."""

    if not hasattr(os, "sched_getaffinity") or not hasattr(os, "sched_setaffinity"):
        raise RuntimeError("CPU affinity requires Linux os.sched_setaffinity()")
    available = set(os.sched_getaffinity(0))
    if cpu_id not in available:
        raise ValueError(
            f"CPU {cpu_id} is unavailable; choose one of {sorted(available)}"
        )
    os.sched_setaffinity(0, {cpu_id})
    return available


def _restore_affinity(cpu_ids: set[int]) -> None:
    """Restore the process CPU affinity saved by :func:`_pin_cpu`."""

    os.sched_setaffinity(0, cpu_ids)


def _load_matrices(
    connection: sqlite3.Connection,
    matrix_ids: list[int] | None,
    size_class: str,
) -> list[tuple[int, int, str, int]]:
    """Load selected matrix IDs, dimensions, values, and expected ESS counts."""

    if matrix_ids:
        requested = list(dict.fromkeys(matrix_ids))
        placeholders = ",".join("?" for _ in requested)
        rows = connection.execute(
            f"""SELECT matrix_id, dimension, matrix, ess_count
                FROM matrices WHERE matrix_id IN ({placeholders})
                ORDER BY matrix_id""",
            requested,
        ).fetchall()
        missing = sorted(set(requested) - {row[0] for row in rows})
        if missing:
            raise ValueError(f"unknown matrix IDs: {missing}")
        return rows

    if size_class == "all":
        rows = connection.execute(
            """SELECT matrix_id, dimension, matrix, ess_count
               FROM matrices ORDER BY matrix_id"""
        ).fetchall()
    else:
        rows = connection.execute(
            """SELECT matrix_id, dimension, matrix, ess_count
               FROM matrices WHERE size_class = ? ORDER BY matrix_id""",
            (size_class,),
        ).fetchall()
    if not rows:
        raise ValueError("matrix selection is empty")
    return rows


def _pybind_arguments(
    matrix: str,
    matrix_id: int,
    mode: str,
    supports_safe_mode: bool,
) -> dict:
    """Return native keyword arguments for one named numerical mode."""

    if mode == "safe" and not supports_safe_mode:
        raise ValueError("this Pybind build does not expose safe mode")
    arguments = {
        "matrix": matrix,
        "include_candidates": False,
        "exact": mode == "exact",
        "full_support": False,
        "enable_logging": False,
        "matrix_id": matrix_id,
    }
    if supports_safe_mode:
        arguments["unsafe"] = mode == "unsafe"
    return arguments


def _pybind_runner(module_dir: Path | None) -> tuple[Callable, Path]:
    """Load one native extension and return its sequential timing callable."""

    if module_dir is None:
        from .core import load_native_module

        native = load_native_module()
    else:
        module_dir = module_dir.resolve()
        if not module_dir.is_dir():
            raise ValueError(f"Pybind module directory does not exist: {module_dir}")
        sys.path.insert(0, str(module_dir))
        invalidate_caches()
        native = import_module("fracessa_core")
    module_path = Path(native.__file__).resolve()
    if module_dir is not None and not module_path.is_relative_to(module_dir):
        raise RuntimeError(
            f"loaded {module_path}, not a fracessa_core module from {module_dir}"
        )
    supports_safe_mode = "unsafe:" in (native.compute_matrix.__doc__ or "")

    def run(matrix_id: int, matrix: str, mode: str) -> tuple[int, int]:
        """Run one matrix through the loaded extension."""

        result = native.compute_matrix(
            **_pybind_arguments(matrix, matrix_id, mode, supports_safe_mode)
        )
        if result["status"] != 0:
            raise RuntimeError(result["error_message"] or "native computation failed")
        return int(result["ess_count"]), int(result["elapsed_ns"])

    return run, module_path


def _parse_cli_output(stdout: str, unit: str) -> tuple[int, int]:
    """Parse ESS count and elapsed time from the CLI's first two output lines."""

    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError("CLI did not print ESS count and timing on two lines")
    try:
        ess_count = int(lines[0])
        raw_time = Decimal(lines[1])
    except (ValueError, InvalidOperation) as error:
        raise RuntimeError(f"invalid CLI timing output: {lines[:2]}") from error
    if ess_count < 0 or not raw_time.is_finite() or raw_time < 0:
        raise RuntimeError(f"invalid CLI timing output: {lines[:2]}")
    elapsed_ns = int(
        (raw_time * _CLI_TO_NS[unit]).to_integral_value(rounding=ROUND_HALF_UP)
    )
    return ess_count, elapsed_ns


def _cli_runner(
    executable: Path,
    unit: str,
    safe_default: bool,
    unsafe_default: bool,
) -> tuple[Callable, Path]:
    """Return a sequential callable for a CLI-only FracESSA build."""

    executable = executable.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"CLI executable is missing or not executable: {executable}")

    def run(matrix_id: int, matrix: str, mode: str) -> tuple[int, int]:
        """Run one matrix through the selected CLI executable."""

        del matrix_id
        if mode == "safe":
            if not safe_default:
                raise ValueError("safe CLI mode requires --safe-default")
            mode_arguments = []
        elif mode == "unsafe":
            mode_arguments = [] if unsafe_default else ["-u"]
        else:
            mode_arguments = ["-e"]

        completed = subprocess.run(
            [str(executable), "-t", *mode_arguments, matrix],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(f"CLI exited with {completed.returncode}: {detail}")
        return _parse_cli_output(completed.stdout, unit)

    return run, executable


def _measure_target(
    runner: Callable,
    matrix_id: int,
    matrix: str,
    mode: str,
    target_ns: int,
) -> tuple[int, int, int, int]:
    """Measure one matrix for about ``target_ns`` and return its native average.

    A pilot taking at least the target is the sole measurement. Faster calls
    use the pilot as warmup, then divide the target by one measured call's wall
    time to choose the batch size.
    """

    pilot_started = perf_counter_ns()
    pilot_ess, pilot_elapsed_ns = runner(matrix_id, matrix, mode)
    pilot_wall_ns = max(1, perf_counter_ns() - pilot_started)
    if pilot_elapsed_ns < 0:
        raise RuntimeError("native timing cannot be negative")
    if pilot_wall_ns >= target_ns:
        return pilot_ess, pilot_elapsed_ns, 1, pilot_wall_ns

    measured_started = perf_counter_ns()
    ess_count, elapsed_ns = runner(matrix_id, matrix, mode)
    first_wall_ns = max(1, perf_counter_ns() - measured_started)
    if ess_count != pilot_ess:
        raise RuntimeError("ESS count changed between timing iterations")
    if elapsed_ns < 0:
        raise RuntimeError("native timing cannot be negative")

    iterations = max(1, (target_ns + first_wall_ns - 1) // first_wall_ns)
    total_elapsed_ns = elapsed_ns
    completed = 1
    while True:
        while completed < iterations:
            current_ess, current_elapsed_ns = runner(matrix_id, matrix, mode)
            if current_ess != ess_count:
                raise RuntimeError("ESS count changed between timing iterations")
            if current_elapsed_ns < 0:
                raise RuntimeError("native timing cannot be negative")
            total_elapsed_ns += current_elapsed_ns
            completed += 1

        measured_wall_ns = max(1, perf_counter_ns() - measured_started)
        if measured_wall_ns >= target_ns:
            break
        average_wall_ns = max(1, measured_wall_ns // iterations)
        remaining_ns = target_ns - measured_wall_ns
        iterations += (remaining_ns + average_wall_ns - 1) // average_wall_ns

    average_ns = (total_elapsed_ns + iterations // 2) // iterations
    return ess_count, average_ns, iterations, measured_wall_ns


def _validate_run(arguments: argparse.Namespace) -> None:
    """Reject inconsistent benchmark command-line options."""

    if not math.isfinite(arguments.target_seconds) or arguments.target_seconds <= 0:
        raise ValueError("--target-seconds must be finite and positive")
    if round(arguments.target_seconds * 1_000_000_000) < 1:
        raise ValueError("--target-seconds is below one nanosecond")
    if len(arguments.mode) != len(set(arguments.mode)):
        raise ValueError("each --mode may be specified only once")
    if arguments.safe_default and arguments.unsafe_default:
        raise ValueError("--safe-default and --unsafe-default are mutually exclusive")
    if arguments.backend == "pybind":
        if arguments.executable is not None:
            raise ValueError("--executable is only valid with --backend cli")
        if arguments.safe_default or arguments.unsafe_default:
            raise ValueError("CLI default-mode flags are invalid with --backend pybind")
    else:
        if arguments.executable is None:
            raise ValueError("--executable is required with --backend cli")
        if arguments.module_dir is not None:
            raise ValueError("--module-dir is only valid with --backend pybind")
    for name in ("build_label", "source_ref", "revision"):
        if not getattr(arguments, name).strip():
            raise ValueError(f"--{name.replace('_', '-')} cannot be empty")


def _run(arguments: argparse.Namespace) -> int:
    """Execute one sequential benchmark build and persist every sample."""

    _validate_run(arguments)
    session = arguments.session or _new_session()
    target_ns = round(arguments.target_seconds * 1_000_000_000)
    database = arguments.database.resolve()
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")

    connection = sqlite3.connect(database)
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        matrices = _load_matrices(
            connection, arguments.matrix_id, arguments.size_class
        )
        existing = connection.execute(
            "SELECT 1 FROM timings WHERE session = ? AND build_label = ? LIMIT 1",
            (session, arguments.build_label),
        ).fetchone()
        if existing:
            raise ValueError(
                f"session {session!r} already contains build {arguments.build_label!r}"
            )

        if arguments.backend == "pybind":
            runner, binary = _pybind_runner(arguments.module_dir)
        else:
            runner, binary = _cli_runner(
                arguments.executable,
                arguments.cli_unit,
                arguments.safe_default,
                arguments.unsafe_default,
            )

        recorded_at = datetime.now(tz=timezone.utc).isoformat()
        binary_sha256 = _sha256(binary)
        machine = _machine_name()
        print(f"session={session}")
        previous_affinity = _pin_cpu(arguments.cpu)
        try:
            for mode in arguments.mode:
                for matrix_id, dimension, values, expected_ess in matrices:
                    matrix = f"{dimension}#{values}"
                    ess_count, elapsed_ns, iterations, measured_wall_ns = (
                        _measure_target(runner, matrix_id, matrix, mode, target_ns)
                    )
                    connection.execute(
                        """INSERT INTO timings (
                               session, recorded_at, machine, cpu_id, comment,
                               build_label, source_ref, revision, binary_sha256,
                               backend, mode, matrix_id, target_ns, iterations,
                               measured_wall_ns, elapsed_ns, ess_count
                           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                        (
                            session,
                            recorded_at,
                            machine,
                            arguments.cpu,
                            arguments.comment,
                            arguments.build_label,
                            arguments.source_ref,
                            arguments.revision,
                            binary_sha256,
                            arguments.backend,
                            mode,
                            matrix_id,
                            target_ns,
                            iterations,
                            measured_wall_ns,
                            elapsed_ns,
                            ess_count,
                        ),
                    )
                    connection.commit()
                    status = "ok" if ess_count == expected_ess else "mismatch"
                    print(
                        f"{mode} matrix={matrix_id} iterations={iterations} "
                        f"average_ns={elapsed_ns} "
                        f"measured_s={measured_wall_ns / 1_000_000_000:.6f} "
                        f"ess={ess_count} expected={expected_ess} {status}",
                        flush=True,
                    )
        finally:
            _restore_affinity(previous_affinity)
    finally:
        connection.close()
    return 0


def _report(arguments: argparse.Namespace) -> int:
    """Print average timings and correctness for one stored session."""

    database = arguments.database.resolve()
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    connection = sqlite3.connect(database)
    try:
        session = arguments.session
        if session == "latest":
            row = connection.execute(
                "SELECT session FROM timings ORDER BY recorded_at DESC LIMIT 1"
            ).fetchone()
            if row is None:
                raise ValueError("the timing table is empty")
            session = row[0]
        rows = connection.execute(
            """SELECT t.build_label, t.backend, t.source_ref, t.revision,
                      t.binary_sha256, t.machine, t.cpu_id, t.comment, t.mode,
                      t.matrix_id, m.is_cs, m.dimension, t.target_ns,
                      t.iterations, t.measured_wall_ns, t.elapsed_ns,
                      t.ess_count, m.ess_count
               FROM timings AS t
               JOIN matrices AS m USING (matrix_id)
               WHERE t.session = ?
               ORDER BY t.build_label, t.mode, t.matrix_id""",
            (session,),
        ).fetchall()
    finally:
        connection.close()
    if not rows:
        raise ValueError(f"unknown timing session: {session}")

    print(f"session={session}")
    builds = {}
    measurements = []
    averages = {}
    for row in rows:
        build, backend, source_ref, revision, digest, machine, cpu, comment = row[:8]
        (
            mode,
            matrix_id,
            is_cs,
            dimension,
            target_ns,
            iterations,
            measured_wall_ns,
            elapsed_ns,
            ess_count,
            expected_ess,
        ) = row[8:]
        builds.setdefault(
            build, (backend, source_ref, revision, digest, machine, cpu, comment)
        )
        measurements.append(
            (
                build,
                mode,
                matrix_id,
                is_cs,
                dimension,
                target_ns,
                iterations,
                measured_wall_ns,
                elapsed_ns,
                ess_count,
                expected_ess,
            )
        )
        averages[build, mode, matrix_id] = elapsed_ns

    for build, metadata in builds.items():
        backend, source_ref, revision, digest, machine, cpu, comment = metadata
        print(
            f"build={build} backend={backend} source_ref={source_ref} "
            f"revision={revision} sha256={digest} machine={machine} "
            f"cpu={cpu} comment={comment!r}"
        )

    print(
        "build\tmode\tmatrix_id\tis_cs\tdimension\ttarget_s\titerations\t"
        "measured_s\taverage_ns\tess\texpected\tgamma_lower_bound\tstatus"
    )
    for measurement in sorted(measurements):
        (
            build,
            mode,
            matrix_id,
            is_cs,
            dimension,
            target_ns,
            iterations,
            measured_wall_ns,
            elapsed_ns,
            ess_count,
            expected_ess,
        ) = measurement
        status = "ok" if ess_count == expected_ess else "mismatch"
        gamma_lower_bound = expected_ess ** (1 / dimension)
        print(
            f"{build}\t{mode}\t{matrix_id}\t{is_cs}\t{dimension}\t"
            f"{target_ns / 1_000_000_000:g}\t"
            f"{iterations}\t{measured_wall_ns / 1_000_000_000:.6f}\t"
            f"{elapsed_ns}\t{ess_count}\t{expected_ess}\t"
            f"{gamma_lower_bound:.6f}\t{status}"
        )

    if arguments.baseline:
        if arguments.baseline not in builds:
            raise ValueError(f"unknown baseline build: {arguments.baseline}")
        ratios = defaultdict(list)
        for (build, mode, matrix_id), value in averages.items():
            baseline = averages.get((arguments.baseline, mode, matrix_id))
            if build != arguments.baseline and baseline:
                ratios[(build, mode)].append(value / baseline)
        print("build\tmode\tshared_matrices\tmedian_ratio_to_baseline")
        for (build, mode), values in sorted(ratios.items()):
            print(f"{build}\t{mode}\t{len(values)}\t{median(values):.6f}")
    return 0


def _parser() -> argparse.ArgumentParser:
    """Build the command-line parser for timing runs and reports."""

    parser = argparse.ArgumentParser(
        description="Record sequential FracESSA timings in SQLite."
    )
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="time one build")
    run.set_defaults(handler=_run)
    run.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    run.add_argument("--backend", choices=("pybind", "cli"), required=True)
    run.add_argument("--module-dir", type=Path)
    run.add_argument("--executable", type=Path)
    run.add_argument("--build-label", required=True)
    run.add_argument("--source-ref", required=True)
    run.add_argument("--revision", required=True)
    run.add_argument(
        "--mode", action="append", choices=("safe", "unsafe", "exact"), required=True
    )
    run.add_argument("--session")
    run.add_argument("--comment", default="")
    run.add_argument("--cpu", type=int, required=True)
    run.add_argument("--target-seconds", type=float, default=1.0)
    run.add_argument(
        "--size-class", choices=("small", "medium", "large", "all"), default="small"
    )
    run.add_argument("--matrix-id", type=int, action="append")
    run.add_argument("--cli-unit", choices=tuple(_CLI_TO_NS), default="ns")
    run.add_argument(
        "--safe-default",
        action="store_true",
        help="declare that a CLI with no mode flag is safe",
    )
    run.add_argument(
        "--unsafe-default",
        action="store_true",
        help="declare that a legacy CLI with no mode flag is unsafe",
    )

    report = commands.add_parser("report", help="show one stored session")
    report.set_defaults(handler=_report)
    report.add_argument("session", nargs="?", default="latest")
    report.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    report.add_argument("--baseline", help="build label used for timing ratios")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the timing command-line interface and return its process status."""

    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        return arguments.handler(arguments)
    except (ImportError, OSError, RuntimeError, sqlite3.Error, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
