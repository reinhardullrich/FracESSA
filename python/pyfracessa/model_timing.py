"""Benchmark production and experimental FracESSA models with persistent workers."""

from __future__ import annotations

import sys

if not __package__ and sys.path and sys.path[0].endswith(("/pyfracessa", "\\pyfracessa")):
    sys.path.pop(0)

import argparse
from collections import deque
from collections.abc import Sequence
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import queue
import re
import sqlite3
import subprocess
import threading
from time import monotonic
from typing import TextIO

if __package__:
    from . import timing
else:  # Model workers use the Python interpreter recorded by their own CMake build.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import timing  # type: ignore[no-redef]


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATABASE = ROOT / "testdata" / "fracessa_testdata.sqlite3"
DEFAULT_RESULTS_DATABASE = ROOT / "experiments" / "model_timings.sqlite3"
STARTUP_TIMEOUT_SECONDS = 30.0
RESULT_SCHEMA = """
CREATE TABLE IF NOT EXISTS model_timings (
    session TEXT NOT NULL CHECK(length(session) > 0),
    model_id TEXT NOT NULL CHECK(length(model_id) > 0 AND model_id = lower(model_id)),
    method TEXT NOT NULL CHECK(method IN ('fast', 'safe')),
    binary_sha256 TEXT NOT NULL CHECK(length(binary_sha256) = 64),
    binary_path TEXT NOT NULL CHECK(length(binary_path) > 0),
    python_path TEXT NOT NULL CHECK(length(python_path) > 0),
    recorded_at TEXT NOT NULL CHECK(length(recorded_at) > 0),
    machine TEXT NOT NULL CHECK(length(machine) > 0),
    cpu_id INTEGER NOT NULL CHECK(cpu_id >= 0),
    comment TEXT NOT NULL DEFAULT '',
    matrix_id INTEGER NOT NULL CHECK(matrix_id > 0),
    status TEXT NOT NULL CHECK(status IN ('ok', 'mismatch', 'timeout', 'error')),
    target_ns INTEGER NOT NULL CHECK(target_ns > 0),
    timeout_ns INTEGER NOT NULL CHECK(timeout_ns > 0),
    iterations INTEGER CHECK(iterations IS NULL OR iterations > 0),
    measured_wall_ns INTEGER CHECK(measured_wall_ns IS NULL OR measured_wall_ns > 0),
    elapsed_ns INTEGER CHECK(elapsed_ns IS NULL OR elapsed_ns > 0),
    ess_count INTEGER CHECK(ess_count IS NULL OR ess_count >= 0),
    safe_fallback TEXT CHECK(
        safe_fallback IS NULL OR safe_fallback IN (
            'precision_span', 'equilibration', 'equilibration_invalid', 'equilibration_non_convergence'
        )
    ),
    message TEXT,
    PRIMARY KEY(session, model_id, method, matrix_id),
    CHECK(
        (status IN ('ok', 'mismatch') AND iterations IS NOT NULL AND measured_wall_ns IS NOT NULL
            AND elapsed_ns IS NOT NULL AND ess_count IS NOT NULL)
        OR (status IN ('timeout', 'error') AND iterations IS NULL AND elapsed_ns IS NULL AND ess_count IS NULL)
    )
) STRICT;
"""
RESULT_UPSERT = """INSERT INTO model_timings (
       session, model_id, method, binary_sha256, binary_path, python_path, recorded_at, machine, cpu_id, comment,
       matrix_id, status, target_ns, timeout_ns, iterations, measured_wall_ns, elapsed_ns, ess_count, safe_fallback, message
   ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
   ON CONFLICT(session, model_id, method, matrix_id) DO UPDATE SET
       binary_sha256 = excluded.binary_sha256,
       binary_path = excluded.binary_path,
       python_path = excluded.python_path,
       recorded_at = excluded.recorded_at,
       machine = excluded.machine,
       cpu_id = excluded.cpu_id,
       comment = excluded.comment,
       status = excluded.status,
       target_ns = excluded.target_ns,
       timeout_ns = excluded.timeout_ns,
       iterations = excluded.iterations,
       measured_wall_ns = excluded.measured_wall_ns,
       elapsed_ns = excluded.elapsed_ns,
       ess_count = excluded.ess_count,
       safe_fallback = excluded.safe_fallback,
       message = excluded.message"""


@dataclass(frozen=True)
class Model:
    model_id: str
    method: str
    module_dir: Path
    python: Path
    binary: Path
    binary_sha256: str


@dataclass
class Worker:
    cpu_id: int
    generation: int
    process: subprocess.Popen
    error_log: TextIO
    error_log_path: Path
    state: str = "starting"
    task: dict | None = None
    started_at: float | None = None
    deadline: float | None = None


def _configured_python(cache: Path) -> Path:
    """Return the Python interpreter used to compile one Pybind model."""

    if not cache.is_file():
        raise ValueError(f"CMake cache does not exist: {cache}")
    values = {}
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = value
    executable = values.get("_Python3_EXECUTABLE:INTERNAL") or values.get("Python3_EXECUTABLE:FILEPATH")
    if not executable:
        raise ValueError(f"CMake cache does not record its Python interpreter: {cache}")
    python = Path(executable).resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise ValueError(f"configured Python interpreter is unavailable: {python}")
    return python


def _parse_model(value: str) -> tuple[str, str]:
    """Parse ``MODEL:METHOD`` without introducing a model registry."""

    try:
        model_id, method = value.rsplit(":", 1)
    except ValueError as error:
        raise ValueError(f"model must use MODEL:METHOD syntax: {value!r}") from error
    if model_id != "production" and re.fullmatch(r"a[1-9][0-9]*", model_id) is None:
        raise ValueError(f"unknown model {model_id!r}; use 'production' or an aN directory under models/")
    if method not in ("fast", "safe"):
        raise ValueError(f"unknown method {method!r}; use 'fast' or 'safe'")
    return model_id, method


def _resolve_model(value: str, root: Path = ROOT) -> Model:
    """Resolve a conventional production or ``models/aN`` Pybind build."""

    model_id, method = _parse_model(value)
    if model_id == "production":
        build = root / "cpp" / "build-benchmark"
        module_dir = build
    else:
        build = root / "models" / model_id / "build"
        module_dir = build / model_id
    candidates = sorted((*module_dir.glob("fracessa_core*.so"), *module_dir.glob("fracessa_core*.pyd")))
    if len(candidates) != 1:
        raise ValueError(f"expected one fracessa_core extension in {module_dir}, found {len(candidates)}")
    binary = candidates[0].resolve()
    return Model(
        model_id=model_id,
        method=method,
        module_dir=module_dir.resolve(),
        python=_configured_python(build / "CMakeCache.txt"),
        binary=binary,
        binary_sha256=timing._sha256(binary),
    )


def _initialize_results_database(path: Path) -> sqlite3.Connection:
    """Open the local results database and ensure its one table exists."""

    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA busy_timeout = 5000")
    connection.executescript(RESULT_SCHEMA)
    return connection


def _worker_main(module_dir: Path, cpu_id: int) -> int:
    """Load one model once and serve timed matrix requests on standard input."""

    try:
        os.sched_setaffinity(0, {cpu_id})
        runner, binary = timing._pybind_runner(module_dir)
        print(json.dumps({"type": "ready", "binary": str(binary)}), flush=True)
    except Exception as error:
        print(json.dumps({"type": "fatal", "message": str(error)}), flush=True)
        return 1

    for line in sys.stdin:
        request = json.loads(line)
        if request.get("type") == "stop":
            return 0
        matrix_id = int(request["matrix_id"])
        print(json.dumps({"type": "started", "matrix_id": matrix_id}), flush=True)
        try:
            ess_count, elapsed_ns, iterations, measured_wall_ns, safe_fallback = timing._measure_target(
                runner,
                matrix_id,
                str(request["matrix"]),
                str(request["method"]),
                int(request["target_ns"]),
                int(request["calibration_ns"]),
            )
            response = {
                "type": "result",
                "matrix_id": matrix_id,
                "ess_count": ess_count,
                "elapsed_ns": elapsed_ns,
                "iterations": iterations,
                "measured_wall_ns": measured_wall_ns,
                "safe_fallback": safe_fallback,
            }
        except Exception as error:
            response = {"type": "error", "matrix_id": matrix_id, "message": str(error)}
        print(json.dumps(response, separators=(",", ":")), flush=True)
    return 0


def _read_worker_output(cpu_id: int, generation: int, stream, messages: queue.Queue) -> None:
    """Forward one worker's complete JSON lines to the scheduler."""

    for line in stream:
        messages.put((cpu_id, generation, line))
    messages.put((cpu_id, generation, None))


def _spawn_worker(model: Model, cpu_id: int, session: str, generation: int, messages: queue.Queue) -> Worker:
    """Start one persistent model process pinned by its worker entry point."""

    log_dir = ROOT / ".local-tmp" / "model-timing"
    log_dir.mkdir(parents=True, exist_ok=True)
    safe_session = re.sub(r"[^A-Za-z0-9_.-]", "_", session)
    error_log_path = log_dir / f"{safe_session}-{model.model_id}-{model.method}-cpu{cpu_id}.log"
    error_log = error_log_path.open("a", encoding="utf-8")
    process = subprocess.Popen(
        [
            str(model.python),
            str(Path(__file__).resolve()),
            "_worker",
            "--module-dir",
            str(model.module_dir),
            "--cpu",
            str(cpu_id),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=error_log,
        text=True,
        bufsize=1,
        cwd=ROOT,
    )
    assert process.stdout is not None
    threading.Thread(
        target=_read_worker_output,
        args=(cpu_id, generation, process.stdout, messages),
        daemon=True,
    ).start()
    return Worker(
        cpu_id=cpu_id,
        generation=generation,
        process=process,
        error_log=error_log,
        error_log_path=error_log_path,
        deadline=monotonic() + STARTUP_TIMEOUT_SECONDS,
    )


def _close_worker(worker: Worker, graceful: bool) -> None:
    """Stop one worker and close its local diagnostic stream."""

    if worker.process.poll() is None and graceful and worker.process.stdin is not None:
        try:
            worker.process.stdin.write('{"type":"stop"}\n')
            worker.process.stdin.flush()
        except BrokenPipeError:
            pass
    if worker.process.poll() is None:
        try:
            worker.process.wait(timeout=1.0 if graceful else 0.0)
        except subprocess.TimeoutExpired:
            worker.process.terminate()
            try:
                worker.process.wait(timeout=1.0)
            except subprocess.TimeoutExpired:
                worker.process.kill()
                worker.process.wait()
    if worker.process.stdin is not None:
        worker.process.stdin.close()
    if worker.process.stdout is not None:
        worker.process.stdout.close()
    worker.error_log.close()


def _store_result(
    connection: sqlite3.Connection,
    session: str,
    model: Model,
    worker: Worker,
    task: dict,
    status: str,
    target_ns: int,
    timeout_ns: int,
    comment: str,
    result: dict | None = None,
    message: str | None = None,
) -> None:
    """Persist one completed, failed, or timed-out model run."""

    result = result or {}
    connection.execute(
        RESULT_UPSERT,
        (
            session,
            model.model_id,
            model.method,
            model.binary_sha256,
            str(model.binary),
            str(model.python),
            datetime.now(tz=timezone.utc).isoformat(),
            timing._machine_name(),
            worker.cpu_id,
            comment,
            task["matrix_id"],
            status,
            target_ns,
            timeout_ns,
            result.get("iterations"),
            result.get("measured_wall_ns"),
            result.get("elapsed_ns"),
            result.get("ess_count"),
            result.get("safe_fallback"),
            message,
        ),
    )
    connection.commit()


def _existing_matrix_ids(connection: sqlite3.Connection, session: str, model: Model) -> set[int]:
    """Return rows already recorded for this exact model binary."""

    rows = connection.execute(
        "SELECT matrix_id, binary_sha256 FROM model_timings WHERE session = ? AND model_id = ? AND method = ?",
        (session, model.model_id, model.method),
    ).fetchall()
    hashes = {row[1] for row in rows}
    if hashes and hashes != {model.binary_sha256}:
        raise ValueError(
            f"session {session!r} already contains a different {model.model_id}:{model.method} binary"
        )
    return {row[0] for row in rows}


def _run_model(
    connection: sqlite3.Connection,
    session: str,
    model: Model,
    matrices: list[tuple],
    cpu_ids: list[int],
    target_ns: int,
    timeout_ns: int,
    comment: str,
    rerun: bool,
) -> None:
    """Run one model over selected matrices with one persistent process per CPU."""

    calibration_index = 4 if model.method == "fast" else 5
    missing = [row[0] for row in matrices if row[calibration_index] is None]
    if missing:
        raise ValueError(f"matrix IDs have no {model.method} calibration: {missing}")
    completed = set() if rerun else _existing_matrix_ids(connection, session, model)
    tasks = deque(
        {
            "matrix_id": matrix_id,
            "dimension": dimension,
            "matrix": f"{dimension}#{values}",
            "expected_ess": expected_ess,
            "calibration_ns": fast_calibration_ns if model.method == "fast" else safe_calibration_ns,
        }
        for matrix_id, dimension, values, expected_ess, fast_calibration_ns, safe_calibration_ns in matrices
        if matrix_id not in completed
    )
    total = len(tasks)
    print(
        f"model={model.model_id} method={model.method} matrices={total} binary_sha256={model.binary_sha256} "
        f"python={model.python} module={model.module_dir}",
        flush=True,
    )
    if not tasks:
        return

    messages: queue.Queue = queue.Queue()
    generation = 0
    workers: dict[int, Worker] = {}
    completed_count = 0
    timeout_seconds = timeout_ns / 1_000_000_000

    def spawn(cpu_id: int) -> None:
        nonlocal generation
        generation += 1
        workers[cpu_id] = _spawn_worker(model, cpu_id, session, generation, messages)

    def remove(cpu_id: int, graceful: bool, restart: bool) -> None:
        worker = workers.pop(cpu_id)
        _close_worker(worker, graceful)
        if restart and tasks:
            spawn(cpu_id)

    def record(worker: Worker, status: str, result: dict | None = None, message: str | None = None) -> None:
        nonlocal completed_count
        assert worker.task is not None
        _store_result(
            connection, session, model, worker, worker.task, status, target_ns, timeout_ns, comment, result, message
        )
        completed_count += 1
        if total <= 50 or completed_count % 25 == 0 or completed_count == total or status != "ok":
            elapsed = None if result is None else result.get("elapsed_ns")
            print(
                f"[{completed_count}/{total}] model={model.model_id} matrix={worker.task['matrix_id']} "
                f"dimension={worker.task['dimension']} cpu={worker.cpu_id} status={status} elapsed_ns={elapsed}",
                flush=True,
            )
        worker.task = None
        worker.started_at = None

    for cpu_id in cpu_ids[: min(len(cpu_ids), total)]:
        spawn(cpu_id)

    try:
        while workers:
            now = monotonic()
            deadlines = [worker.deadline for worker in workers.values() if worker.deadline is not None]
            wait_seconds = 0.25 if not deadlines else max(0.0, min(0.25, min(deadlines) - now))
            try:
                cpu_id, message_generation, line = messages.get(timeout=wait_seconds)
            except queue.Empty:
                line = ""
                cpu_id = -1
                message_generation = -1

            worker = workers.get(cpu_id)
            if worker is not None and worker.generation == message_generation:
                if line is None:
                    if worker.state == "starting":
                        raise RuntimeError(f"{model.model_id} worker exited during startup; see {worker.error_log_path}")
                    if worker.task is not None:
                        record(worker, "error", message=f"worker exited; see {worker.error_log_path}")
                    remove(cpu_id, graceful=False, restart=True)
                else:
                    try:
                        message = json.loads(line)
                    except json.JSONDecodeError:
                        if worker.state == "starting":
                            raise RuntimeError(f"{model.model_id} worker printed invalid startup output; see {worker.error_log_path}")
                        if worker.task is not None:
                            record(worker, "error", message=f"invalid worker output; see {worker.error_log_path}")
                        remove(cpu_id, graceful=False, restart=True)
                    else:
                        message_type = message.get("type")
                        if message_type == "ready" and worker.state == "starting":
                            if Path(message["binary"]).resolve() != model.binary:
                                raise RuntimeError(f"worker loaded {message['binary']}, expected {model.binary}")
                            worker.state = "idle"
                            worker.deadline = None
                        elif message_type == "started" and worker.state == "assigned":
                            if int(message["matrix_id"]) != worker.task["matrix_id"]:
                                raise RuntimeError("worker started the wrong matrix")
                            worker.state = "running"
                            worker.started_at = monotonic()
                            worker.deadline = worker.started_at + timeout_seconds
                        elif message_type in ("result", "error") and worker.state == "running":
                            if int(message["matrix_id"]) != worker.task["matrix_id"]:
                                raise RuntimeError("worker returned the wrong matrix")
                            if monotonic() > worker.deadline:
                                record(
                                    worker,
                                    "timeout",
                                    result={"measured_wall_ns": max(1, int(message.get("measured_wall_ns", timeout_ns)))},
                                    message=f"exceeded {timeout_seconds:g} seconds",
                                )
                            elif message_type == "error":
                                record(worker, "error", message=str(message["message"]))
                            else:
                                status = "ok" if int(message["ess_count"]) == worker.task["expected_ess"] else "mismatch"
                                record(worker, status, result=message)
                            worker.state = "idle"
                            worker.deadline = None
                        elif message_type == "fatal" and worker.state == "starting":
                            raise RuntimeError(f"{model.model_id} worker failed: {message['message']}")
                        else:
                            raise RuntimeError(f"unexpected worker message in state {worker.state}: {message!r}")

            now = monotonic()
            for current_cpu, current in list(workers.items()):
                if current.deadline is None or now < current.deadline:
                    continue
                if current.state == "running":
                    measured_wall_ns = max(1, round((now - current.started_at) * 1_000_000_000))
                    record(
                        current,
                        "timeout",
                        result={"measured_wall_ns": measured_wall_ns},
                        message=f"exceeded {timeout_seconds:g} seconds",
                    )
                    remove(current_cpu, graceful=False, restart=True)
                elif current.state == "assigned":
                    record(current, "error", message="worker did not acknowledge the matrix")
                    remove(current_cpu, graceful=False, restart=True)
                else:
                    raise RuntimeError(f"worker on CPU {current_cpu} did not start; see {current.error_log_path}")

            for current_cpu, current in list(workers.items()):
                if current.state != "idle":
                    continue
                if not tasks:
                    remove(current_cpu, graceful=True, restart=False)
                    continue
                task = tasks.popleft()
                request = {
                    "type": "run",
                    "matrix_id": task["matrix_id"],
                    "matrix": task["matrix"],
                    "method": model.method,
                    "target_ns": target_ns,
                    "calibration_ns": task["calibration_ns"],
                }
                assert current.process.stdin is not None
                try:
                    current.process.stdin.write(json.dumps(request, separators=(",", ":")) + "\n")
                    current.process.stdin.flush()
                except BrokenPipeError:
                    tasks.appendleft(task)
                    remove(current_cpu, graceful=False, restart=True)
                    continue
                current.task = task
                current.state = "assigned"
                current.deadline = monotonic() + STARTUP_TIMEOUT_SECONDS
    finally:
        for worker in list(workers.values()):
            _close_worker(worker, graceful=False)


def _run(arguments: argparse.Namespace) -> int:
    """Run selected models and store resumable per-matrix results."""

    if not math.isfinite(arguments.target_seconds) or arguments.target_seconds <= 0:
        raise ValueError("--target-seconds must be finite and positive")
    if not math.isfinite(arguments.timeout_seconds) or arguments.timeout_seconds <= 0:
        raise ValueError("--timeout-seconds must be finite and positive")
    if len(arguments.cpus) != len(set(arguments.cpus)):
        raise ValueError("--cpus must not contain duplicates")
    if arguments.parent_cpu in arguments.cpus:
        raise ValueError("--parent-cpu must not also be a worker CPU")
    available = set(os.sched_getaffinity(0))
    unavailable = sorted(({arguments.parent_cpu} | set(arguments.cpus)) - available)
    if unavailable:
        raise ValueError(f"CPU IDs are unavailable: {unavailable}")

    models = [_resolve_model(value) for value in arguments.models]
    model_keys = [(model.model_id, model.method) for model in models]
    if len(model_keys) != len(set(model_keys)):
        raise ValueError("each MODEL:METHOD may be specified only once")
    target_ns = round(arguments.target_seconds * 1_000_000_000)
    timeout_ns = round(arguments.timeout_seconds * 1_000_000_000)
    if target_ns < 1:
        raise ValueError("--target-seconds is below one nanosecond")
    if timeout_ns < 1:
        raise ValueError("--timeout-seconds is below one nanosecond")
    session = arguments.session or timing._new_session()
    print(f"session={session}")

    database = arguments.database.resolve()
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    with closing(sqlite3.connect(database)) as corpus:
        matrices = [
            row for row in timing._load_matrices(corpus, arguments.matrix_id, arguments.size_class) if row[1] != 2
        ]
    if not matrices:
        raise ValueError("matrix selection is empty after excluding dimension 2")

    previous_affinity = timing._pin_cpu(arguments.parent_cpu)
    results = _initialize_results_database(arguments.results_database.resolve())
    try:
        for model in models:
            _run_model(
                results,
                session,
                model,
                matrices,
                arguments.cpus,
                target_ns,
                timeout_ns,
                arguments.comment,
                arguments.rerun,
            )
    finally:
        results.close()
        timing._restore_affinity(previous_affinity)
    return 0


def _report(arguments: argparse.Namespace) -> int:
    """Report statuses and median speedups for one stored model session."""

    results_path = arguments.results_database.resolve()
    if not results_path.is_file():
        raise ValueError(f"results database does not exist: {results_path}")
    with closing(sqlite3.connect(results_path)) as connection:
        session = arguments.session
        if session == "latest":
            row = connection.execute(
                "SELECT session FROM model_timings ORDER BY recorded_at DESC LIMIT 1"
            ).fetchone()
            if row is None:
                raise ValueError("the model timing table is empty")
            session = row[0]
        rows = connection.execute(
            """SELECT model_id, method, binary_sha256, matrix_id, status, elapsed_ns, ess_count, safe_fallback
               FROM model_timings WHERE session = ? ORDER BY model_id, method, matrix_id""",
            (session,),
        ).fetchall()
    if not rows:
        raise ValueError(f"unknown session: {session}")

    print(f"session={session}")
    summaries = {}
    for model_id, method, digest, _matrix_id, status, _elapsed, _ess, _fallback in rows:
        summary = summaries.setdefault((model_id, method, digest), {"ok": 0, "mismatch": 0, "timeout": 0, "error": 0})
        summary[status] += 1
    print("model\tmethod\tbinary_sha256\tok\tmismatch\ttimeout\terror")
    for (model_id, method, digest), statuses in sorted(summaries.items()):
        print(
            f"{model_id}\t{method}\t{digest}\t{statuses['ok']}\t{statuses['mismatch']}\t"
            f"{statuses['timeout']}\t{statuses['error']}"
        )

    if arguments.baseline:
        baseline = _parse_model(arguments.baseline)
        usable = {
            (model_id, method, matrix_id): elapsed_ns
            for model_id, method, _digest, matrix_id, status, elapsed_ns, _ess, fallback in rows
            if status == "ok" and not (method == "fast" and fallback is not None)
        }
        baseline_times = {
            matrix_id: elapsed_ns
            for (model_id, method, matrix_id), elapsed_ns in usable.items()
            if (model_id, method) == baseline
        }
        print("model\tmethod\tshared_matrices\tmedian_speedup_vs_baseline")
        for model_id, method in sorted({(row[0], row[1]) for row in rows}):
            if (model_id, method) == baseline:
                continue
            speedups = [
                baseline_times[matrix_id] / elapsed_ns
                for (current_model, current_method, matrix_id), elapsed_ns in usable.items()
                if (current_model, current_method) == (model_id, method) and matrix_id in baseline_times
            ]
            if speedups:
                print(f"{model_id}\t{method}\t{len(speedups)}\t{timing.median(speedups):.6f}")
    return 0


def _parser() -> argparse.ArgumentParser:
    """Build the multi-model timing command-line parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="benchmark one or more MODEL:METHOD builds")
    run.set_defaults(handler=_run)
    run.add_argument("models", nargs="+")
    run.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    run.add_argument("--results-database", type=Path, default=DEFAULT_RESULTS_DATABASE)
    run.add_argument("--session")
    run.add_argument("--comment", default="")
    run.add_argument("--parent-cpu", type=int, required=True)
    run.add_argument("--cpus", type=int, nargs="+", required=True)
    run.add_argument("--target-seconds", type=float, default=0.5)
    run.add_argument("--timeout-seconds", type=float, required=True)
    run.add_argument("--size-class", choices=("small", "medium", "large", "super_large", "all"), default="all")
    run.add_argument("--matrix-id", type=int, action="append")
    run.add_argument("--rerun", action="store_true")

    report = commands.add_parser("report", help="summarize one model timing session")
    report.set_defaults(handler=_report)
    report.add_argument("session", nargs="?", default="latest")
    report.add_argument("--results-database", type=Path, default=DEFAULT_RESULTS_DATABASE)
    report.add_argument("--baseline", help="MODEL:METHOD used as the speedup baseline")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the public command or the private persistent-worker entry point."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "_worker":
        worker_parser = argparse.ArgumentParser(add_help=False)
        worker_parser.add_argument("_worker")
        worker_parser.add_argument("--module-dir", type=Path, required=True)
        worker_parser.add_argument("--cpu", type=int, required=True)
        worker_arguments = worker_parser.parse_args(arguments)
        return _worker_main(worker_arguments.module_dir, worker_arguments.cpu)

    parser = _parser()
    parsed = parser.parse_args(arguments)
    try:
        return parsed.handler(parsed)
    except (ImportError, OSError, RuntimeError, sqlite3.Error, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
