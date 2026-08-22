from contextlib import redirect_stdout
from io import StringIO
import json
import os
from pathlib import Path
from types import SimpleNamespace
import queue
import sqlite3
import sys
import tempfile
import unittest
from unittest import mock

from pyfracessa import model_timing


class ModelTimingTests(unittest.TestCase):
    def test_conventional_models_need_no_registry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for build, module in (
                (root / "cpp" / "build-benchmark", root / "cpp" / "build-benchmark"),
                (root / "models" / "a3" / "build", root / "models" / "a3" / "build" / "a3"),
            ):
                module.mkdir(parents=True)
                (module / "fracessa_core.test.so").write_bytes(b"model")
                (build / "CMakeCache.txt").write_text(f"_Python3_EXECUTABLE:INTERNAL={sys.executable}\n")

            production = model_timing._resolve_model("production:fast", root)
            a3 = model_timing._resolve_model("a3:safe", root)

        self.assertEqual((production.model_id, production.method), ("production", "fast"))
        self.assertEqual((a3.model_id, a3.method), ("a3", "safe"))
        self.assertEqual(len(a3.binary_sha256), 64)

    def test_results_resume_only_the_same_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            connection = model_timing._initialize_results_database(Path(directory) / "results.sqlite3")
            model = model_timing.Model("a1", "safe", Path("module"), Path(sys.executable), Path("binary"), "a" * 64)
            task = {"matrix_id": 7}
            model_timing._store_result(
                connection,
                "session",
                model,
                SimpleNamespace(cpu_id=2),
                task,
                "ok",
                500,
                1_000,
                "",
                {"iterations": 2, "measured_wall_ns": 600, "elapsed_ns": 100, "ess_count": 3},
            )
            self.assertEqual(model_timing._existing_matrix_ids(connection, "session", model), {7})
            changed = model_timing.Model("a1", "safe", Path("module"), Path(sys.executable), Path("binary"), "b" * 64)
            with self.assertRaisesRegex(ValueError, "different a1:safe binary"):
                model_timing._existing_matrix_ids(connection, "session", changed)
            connection.close()

    def test_report_excludes_fast_fallbacks_from_comparisons(self):
        with tempfile.TemporaryDirectory() as directory:
            results_path = Path(directory) / "results.sqlite3"
            connection = model_timing._initialize_results_database(results_path)
            rows = []
            for model, method, matrix_id, elapsed_ns, fallback in (
                ("production", "fast", 1, 100, None),
                ("production", "fast", 2, 100, "precision_span"),
                ("a2", "safe", 1, 50, None),
                ("a2", "safe", 2, 50, None),
            ):
                rows.append(
                    (
                        "session", model, method, "a" * 64, "module", sys.executable, "now", "machine", 2, "", matrix_id,
                        "ok", 500, 1_000, 1, 100, elapsed_ns, 1, fallback, None,
                    )
                )
            connection.executemany(model_timing.RESULT_UPSERT, rows)
            connection.commit()
            connection.close()

            output = StringIO()
            with redirect_stdout(output):
                model_timing._report(
                    SimpleNamespace(session="session", results_database=results_path, baseline="production:fast")
                )

        self.assertIn("a2\tsafe\t1\t2.000000", output.getvalue())

    def test_worker_keeps_one_model_loaded_for_repeated_native_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            module_dir = root / "module"
            module_dir.mkdir()
            binary = module_dir / "fracessa_core.py"
            binary.write_text(
                "def compute_matrix(*args, **kwargs):\n"
                "    return {'status': 0, 'error_message': '', 'ess_count': 3, 'elapsed_ns': 100, 'safe_fallback': None}\n"
                "compute_matrix.__doc__ = 'method:'\n"
            )
            model = model_timing.Model("a1", "safe", module_dir, Path(sys.executable), binary, "a" * 64)
            messages = queue.Queue()
            cpu_id = min(os.sched_getaffinity(0))
            with mock.patch.object(model_timing, "ROOT", root):
                worker = model_timing._spawn_worker(model, cpu_id, "test", 1, messages)
                try:
                    self.assertIsNotNone(messages.get(timeout=10)[2])
                    assert worker.process.stdin is not None
                    worker.process.stdin.write(
                        '{"type":"run","matrix_id":1,"matrix":"2#0,1,0","method":"safe",'
                        '"target_ns":200,"calibration_ns":100}\n'
                    )
                    worker.process.stdin.flush()
                    started = messages.get(timeout=10)
                    result = messages.get(timeout=10)
                    self.assertEqual(json.loads(started[2])["type"], "started")
                    self.assertEqual(json.loads(result[2])["iterations"], 2)
                finally:
                    model_timing._close_worker(worker, graceful=True)

    def test_model_run_records_and_stops_a_timed_out_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            module_dir = root / "module"
            module_dir.mkdir()
            binary = module_dir / "fracessa_core.py"
            binary.write_text(
                "import time\n"
                "def compute_matrix(*args, **kwargs):\n"
                "    time.sleep(0.1)\n"
                "    return {'status': 0, 'error_message': '', 'ess_count': 1, 'elapsed_ns': 100, 'safe_fallback': None}\n"
                "compute_matrix.__doc__ = 'method:'\n"
            )
            model = model_timing.Model("a1", "safe", module_dir, Path(sys.executable), binary, "a" * 64)
            connection = model_timing._initialize_results_database(root / "results.sqlite3")
            cpu_id = min(os.sched_getaffinity(0))
            with mock.patch.object(model_timing, "ROOT", root):
                with redirect_stdout(StringIO()):
                    model_timing._run_model(
                        connection,
                        "session",
                        model,
                        [(1, 3, "0,0,0,0,0,0", 1, 100, 100)],
                        [cpu_id],
                        100,
                        10_000_000,
                        "",
                        False,
                    )
            row = connection.execute("SELECT status, elapsed_ns FROM model_timings").fetchone()
            connection.close()

        self.assertEqual(row, ("timeout", None))


if __name__ == "__main__":
    unittest.main()
