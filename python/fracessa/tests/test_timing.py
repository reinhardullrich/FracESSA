from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest import mock

from fracessa import timing


class TimingTests(unittest.TestCase):
    def test_cli_seconds_are_normalized_and_adaptively_averaged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "timing.sqlite3"
            schema = timing.DEFAULT_DATABASE.with_name("schema.sql").read_text()
            connection = sqlite3.connect(database)
            connection.executescript(schema)
            connection.execute(
                """INSERT INTO matrices (
                       matrix_id, dimension, size_class, is_cs, matrix,
                       candidate_count, ess_count, candidate_structure,
                       ess_structure, origin
                   ) VALUES (1, 2, 'small', 0, '0,1,0', 4, 4, '{}', '{}', 'unit')"""
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
                    side_effect=[
                        0,
                        100_000,
                        200_000,
                        300_000,
                        800_000,
                        1_200_000,
                    ],
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
                        "--mode",
                        "unsafe",
                        "--unsafe-default",
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
                    "unsafe",
                    1_000_000,
                    17,
                    1_000_000,
                    123_000,
                    4,
                    "main",
                    "abc123",
                ),
            )
            self.assertIn("iterations=17 average_ns=123000", output.getvalue())

    def test_measurement_at_target_uses_only_the_pilot(self):
        calls = []

        def runner(matrix_id, matrix, mode):
            calls.append((matrix_id, matrix, mode))
            return 7, 1_100_000_000

        with mock.patch(
            "fracessa.timing.perf_counter_ns", side_effect=[0, 1_100_000_000]
        ):
            result = timing._measure_target(
                runner, 3, "2#0,1,0", "exact", 1_000_000_000
            )

        self.assertEqual(result, (7, 1_100_000_000, 1, 1_100_000_000))
        self.assertEqual(len(calls), 1)

    def test_pybind_modes_use_explicit_safe_and_unsafe_flags(self):
        safe = timing._pybind_arguments("2#0,1,0", 1, "safe", True)
        unsafe = timing._pybind_arguments("2#0,1,0", 1, "unsafe", True)
        exact = timing._pybind_arguments("2#0,1,0", 1, "exact", True)

        self.assertFalse(safe["exact"])
        self.assertFalse(safe["unsafe"])
        self.assertTrue(unsafe["unsafe"])
        self.assertTrue(exact["exact"])
        with self.assertRaisesRegex(ValueError, "does not expose safe mode"):
            timing._pybind_arguments("2#0,1,0", 1, "safe", False)

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
                           'abc123', ?, 'pybind', 'exact', 1, 1000000000, 1,
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
                "current\texact\t1\t1\t3\t1\t1\t0.000000\t100\t8\t8\t"
                "2.000000\tok",
                output.getvalue(),
            )
