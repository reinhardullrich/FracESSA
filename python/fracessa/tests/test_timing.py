from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest import mock

from fracessa import timing


class TimingTests(unittest.TestCase):
    def test_matrix_selection_skips_rows_without_candidate_baselines(self):
        with sqlite3.connect(":memory:") as connection:
            connection.executescript(
                timing.DEFAULT_DATABASE.with_name("schema.sql").read_text()
            )
            connection.executemany(
                """INSERT INTO matrices (
                       matrix_id, dimension, size_class, is_cs, matrix,
                       candidate_count, ess_count, candidate_structure,
                       ess_structure, origin
                   ) VALUES (?, 2, 'small', 0, '0,1,0', ?, ?, ?, ?, 'unit')""",
                [
                    (1, 1, 1, "{}", "{}"),
                    (2, None, None, None, None),
                ],
            )

            self.assertEqual(
                timing._load_matrices(connection, None, "all"),
                [(1, 2, "0,1,0", 1, None, None)],
            )
            with self.assertRaisesRegex(ValueError, "no candidate baseline"):
                timing._load_matrices(connection, [2], "all")

    def test_cli_seconds_are_normalized_and_calibrated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "timing.sqlite3"
            schema = timing.DEFAULT_DATABASE.with_name("schema.sql").read_text()
            connection = sqlite3.connect(database)
            connection.executescript(schema)
            connection.execute(
                """INSERT INTO matrices (
                       matrix_id, dimension, size_class, is_cs, matrix,
                       fast_calibration_ns,
                       candidate_count, ess_count, candidate_structure,
                       ess_structure, origin
                   ) VALUES (1, 2, 'small', 0, '0,1,0', 123000, 4, 4, '{}', '{}', 'unit')"""
            )
            connection.commit()
            connection.close()

            executable = root / "fake-fracessa"
            executable.write_text("#!/bin/sh\nprintf '4\\n0.000123\\n'\n")
            executable.chmod(0o755)

            output = StringIO()
            with (
                mock.patch("fracessa.timing._pin_cpu", return_value={0}),
                mock.patch("fracessa.timing._restore_affinity"),
                mock.patch(
                    "fracessa.timing.perf_counter_ns",
                    side_effect=[0, 1_200_000],
                ),
                redirect_stdout(output),
            ):
                status = timing.main(
                    [
                        "run",
                        "--database",
                        str(database),
                        "--backend",
                        "cli",
                        "--executable",
                        str(executable),
                        "--build-label",
                        "legacy",
                        "--source-ref",
                        "main",
                        "--revision",
                        "abc123",
                        "--method",
                        "fast",
                        "--fast-default",
                        "--cli-unit",
                        "s",
                        "--cpu",
                        "0",
                        "--target-seconds",
                        "0.001",
                    ]
                )

            self.assertEqual(status, 0)
            with sqlite3.connect(database) as connection:
                row = connection.execute(
                    """SELECT backend, mode, target_ns, iterations,
                              measured_wall_ns, elapsed_ns, ess_count,
                              source_ref, revision
                       FROM timings"""
                ).fetchone()
            self.assertEqual(
                row,
                (
                    "cli",
                    "fast",
                    1_000_000,
                    9,
                    1_200_000,
                    123_000,
                    4,
                    "main",
                    "abc123",
                ),
            )
            self.assertIn("iterations=9 calibration_us=123.000 median_ns=123000", output.getvalue())

    def test_native_measurement_at_target_uses_one_run(self):
        calls = []

        def runner(matrix_id, matrix, method):
            calls.append((matrix_id, matrix, method))
            return 7, 1_100_000_000

        with mock.patch(
            "fracessa.timing.perf_counter_ns", side_effect=[0, 5]
        ):
            result = timing._measure_target(runner, 3, "2#0,1,0", "safe", 1_000_000_000, 1_100_000_000)

        self.assertEqual(result, (7, 1_100_000_000, 1, 5))
        self.assertEqual(len(calls), 1)

    def test_timeout_calibration_uses_one_run(self):
        calls = []

        def runner(matrix_id, matrix, method):
            calls.append((matrix_id, matrix, method))
            return 7, 12_000_000_000

        with mock.patch("fracessa.timing.perf_counter_ns", side_effect=[0, 12_000_000_001]):
            result = timing._measure_target(runner, 3, "2#0,1,0", "fast", 500_000_000, -1)

        self.assertEqual(result, (7, 12_000_000_000, 1, 12_000_000_001))
        self.assertEqual(len(calls), 1)

    def test_native_duration_sizes_sample_and_result_is_median(self):
        elapsed = iter([100, 1_000, 90, 100])

        def runner(matrix_id, matrix, method):
            return 7, next(elapsed)

        with mock.patch(
            "fracessa.timing.perf_counter_ns", side_effect=[0, 10_000]
        ):
            result = timing._measure_target(runner, 3, "2#0,1,0", "safe", 400, 100)

        self.assertEqual(result, (7, 100, 4, 10_000))

    def test_default_target_is_half_a_second(self):
        arguments = timing._parser().parse_args(
            [
                "run", "--backend", "pybind", "--build-label", "current",
                "--source-ref", "main", "--revision", "abc123", "--method", "fast", "--cpu", "2",
            ]
        )
        self.assertEqual(arguments.target_seconds, 0.5)

    def test_run_rejects_a_missing_calibration_before_loading_the_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "timing.sqlite3"
            with sqlite3.connect(database) as connection:
                connection.executescript(timing.DEFAULT_DATABASE.with_name("schema.sql").read_text())
                connection.execute(
                    """INSERT INTO matrices (
                           matrix_id, dimension, size_class, is_cs, matrix,
                           candidate_count, ess_count, candidate_structure,
                           ess_structure, origin
                       ) VALUES (1, 3, 'small', 0, '0,1,0,0,0,0', 1, 1, '{}', '{}', 'unit')"""
                )

            arguments = timing._parser().parse_args(
                [
                    "run", "--database", str(database), "--backend", "pybind", "--build-label", "current",
                    "--source-ref", "main", "--revision", "abc123", "--method", "fast", "--cpu", "2",
                ]
            )
            with self.assertRaisesRegex(ValueError, "matrix IDs have no fast calibration: \\[1\\]"):
                timing._run(arguments)

    def test_pybind_method_and_historical_interfaces_are_mapped(self):
        for method in ("fast", "safe"):
            arguments = timing._pybind_arguments("2#0,1,0", 1, method, "method")
            self.assertEqual(arguments["method"], method)

        self.assertEqual(
            timing._pybind_arguments("2#0,1,0", 1, "fast", "mode")["mode"],
            "unsafe",
        )
        self.assertEqual(
            timing._pybind_arguments("2#0,1,0", 1, "safe", "mode")["mode"],
            "exact",
        )

        old = timing._pybind_arguments("2#0,1,0", 1, "fast", "booleans")
        self.assertFalse(old["exact"])
        self.assertTrue(old["unsafe"])

    def test_report_includes_matrix_growth_base(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "timing.sqlite3"
            schema = timing.DEFAULT_DATABASE.with_name("schema.sql").read_text()
            with sqlite3.connect(database) as connection:
                connection.executescript(schema)
                connection.execute(
                    """INSERT INTO matrices (
                           matrix_id, dimension, size_class, is_cs, matrix,
                           candidate_count, ess_count, candidate_structure,
                           ess_structure, origin
                       ) VALUES (1, 3, 'small', 1, '1', 8, 8, '{}', '{}', 'unit')"""
                )
                connection.execute(
                    """INSERT INTO timings (
                           session, recorded_at, machine, cpu_id, comment,
                           build_label, source_ref, revision, binary_sha256,
                           backend, mode, matrix_id, target_ns, iterations,
                           measured_wall_ns, elapsed_ns, ess_count
                       ) VALUES (
                           'test', 'now', 'machine', 0, '', 'current', 'main',
                           'abc123', ?, 'pybind', 'safe', 1, 1000000000, 1,
                           123, 100, 8
                       )""",
                    ("0" * 64,),
                )

            output = StringIO()
            with redirect_stdout(output):
                status = timing.main(
                    ["report", "test", "--database", str(database)]
                )

            self.assertEqual(status, 0)
            self.assertIn(
                "matrix_id\tis_cs\tdimension\ttarget_s", output.getvalue()
            )
            self.assertIn(
                "current\tsafe\t1\t1\t3\t1\t1\t0.000000\t100\t8\t8\t"
                "2.000000\tok",
                output.getvalue(),
            )
