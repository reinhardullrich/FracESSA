from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest

from pyfracessa import (
    MPConfig,
    Matrix,
    RunConfig,
    run,
    run_multiprocessing,
)
from pyfracessa.core import load_native_module


class NativeIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_native_module()

    def test_run_native(self):
        matrix_id = (1 << 63) - 1
        result = run(
            "safe",
            Matrix(matrix_id=matrix_id, matrix="2#0,1,0"),
            RunConfig(include_candidates=True),
            run_id="native_single",
        )

        self.assertNotIn("success", result)
        self.assertEqual(result["status"], 0)
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(result["ess_count"], 1)
        self.assertEqual(result["candidate_structure"], {2: 1})
        self.assertEqual(result["ess_structure"], {2: 1})
        self.assertIsNone(result["safe_fallback"])
        self.assertEqual(result["matrix_id"], matrix_id)
        self.assertEqual(result["candidates"][0]["matrix_id"], matrix_id)
        self.assertIsNone(result["candidates"][0]["multiplier"])

    def test_fast_and_safe_route_native(self):
        database = Path(__file__).resolve().parents[2] / "testdata/fracessa_testdata.sqlite3"
        with closing(sqlite3.connect(database)) as connection:
            rows = connection.execute(
                "SELECT matrix_id, dimension, matrix FROM matrices WHERE matrix_id IN (46, 2208) ORDER BY matrix_id"
            ).fetchall()

        for matrix_id, dimension, values in rows:
            with self.subTest(matrix_id=matrix_id):
                matrix = Matrix(matrix_id=matrix_id, matrix=f"{dimension}#{values}")
                safe = run("safe", matrix, RunConfig(include_candidates=True), run_id=f"native_safe_{matrix_id}")
                fast = run("fast", matrix, RunConfig(include_candidates=True), run_id=f"native_fast_{matrix_id}")

                self.assertEqual(safe["status"], 0)
                self.assertEqual(fast["status"], 0)
                self.assertEqual(safe["ess_count"], 1)
                self.assertEqual(fast["ess_count"], 1)

    def test_whole_matrix_safe_fallback_is_exposed(self):
        native = load_native_module()
        database = Path(__file__).resolve().parents[2] / "testdata/fracessa_testdata.sqlite3"
        with closing(sqlite3.connect(database)) as connection:
            dimension, values = connection.execute(
                "SELECT dimension, matrix FROM matrices WHERE matrix_id = 2109"
            ).fetchone()
        cases = {
            "2#0,1,0": None,
            "2#1,1000000000,1": "precision_span",
            "2#0,0,0": None,
            f"{dimension}#{values}": "equilibration_invalid",
            "5#0,0,0,0,1,0,0,0,1,0,1,0,0,1,0": "equilibration_non_convergence",
        }

        for matrix, expected in cases.items():
            with self.subTest(matrix=matrix):
                self.assertEqual(native.classify_safe_fallback(matrix), expected)
                result = native.compute_matrix("fast", matrix)
                self.assertEqual(result["safe_fallback"], expected)

    def test_removed_or_unknown_native_method_is_rejected(self):
        native = load_native_module()
        for method in ("test", "unknown"):
            with self.subTest(method=method):
                result = native.compute_matrix(method, "2#0,1,0")
                self.assertEqual(result["status"], 4)

    def test_removed_native_mode_keyword_is_rejected(self):
        native = load_native_module()
        with self.assertRaises(TypeError):
            native.compute_matrix(matrix="2#0,1,0", mode="exact")

    def test_native_candidate_contract(self):
        native = load_native_module()
        result = native.compute_matrix("safe", "2#0,1,0", include_candidates=True)
        candidate = result["candidates"][0]

        self.assertEqual(
            candidate,
            {
                "candidate_id": 1,
                "vector": "1/2,1/2",
                "support": 3,
                "support_size": 2,
                "extended_support": 3,
                "extended_support_size": 2,
                "multiplier": None,
                "is_ess": True,
                "stability": "T_reduced_hessian_nd",
                "payoff": "1/2",
                "payoff_dbl": 0.5,
            },
        )
        self.assertEqual(
            {key: type(value) for key, value in candidate.items()},
            {
                "candidate_id": int,
                "vector": str,
                "support": int,
                "support_size": int,
                "extended_support": int,
                "extended_support_size": int,
                "multiplier": type(None),
                "is_ess": bool,
                "stability": str,
                "payoff": str,
                "payoff_dbl": float,
            },
        )

    def test_circular_native_returns_one_weighted_representative(self):
        result = run(
            "safe",
            Matrix(matrix_id=2, matrix="5#1,3"),
            RunConfig(include_candidates=True),
            run_id="native_circular",
        )

        self.assertEqual(result["ess_count"], 5)
        self.assertEqual(result["candidate_count"], 5)
        self.assertEqual(result["candidate_structure"], {3: 5})
        self.assertEqual(result["ess_structure"], {3: 5})
        self.assertEqual(len(result["candidates"]), 1)
        self.assertEqual(result["candidates"][0]["multiplier"], 5)

        without_rows = run("safe", Matrix(matrix_id=2, matrix="5#1,3"), RunConfig(include_candidates=False))
        self.assertEqual(without_rows["candidate_count"], 5)
        self.assertEqual(without_rows["candidate_structure"], {3: 5})
        self.assertEqual(without_rows["ess_structure"], {3: 5})
        self.assertEqual(without_rows["candidates"], [])

    def test_multiword_native_returns_python_integers(self):
        dimension = 65
        matrix = f"{dimension}#" + ",".join(["1"] * (dimension // 2))
        expected_support = (1 << dimension) - 1

        for method in ("safe", "fast"):
            with self.subTest(method=method):
                result = run(
                    method,
                    Matrix(matrix_id=65, matrix=matrix),
                    RunConfig(full_support=True, include_candidates=True),
                    run_id=f"native_multiword_{method}",
                )
                self.assertEqual(result["status"], 0)
                self.assertEqual(result["candidate_structure"], {dimension: 1})
                self.assertEqual(result["ess_structure"], {dimension: 1})
                self.assertEqual(result["candidates"][0]["support"], expected_support)
                self.assertEqual(result["candidates"][0]["extended_support"], expected_support)
                self.assertIs(type(result["candidates"][0]["support"]), int)

    def test_invalid_matrix_strings_return_parser_error(self):
        invalid_matrices = {
            "2##0": "Multiple '#' characters found in matrix string",
            "2#0,1": "Expected 1 (CS) or 3 (Sym) values, got 2",
            "2#0,1/0,0": "Rational denominator cannot be zero: 1/0",
        }

        for matrix, error_message in invalid_matrices.items():
            with self.subTest(matrix=matrix):
                result = run("safe", Matrix(matrix_id=1, matrix=matrix), run_id="invalid")
                self.assertEqual(result["status"], 1)
                self.assertEqual(result["error_message"], error_message)

    def test_run_multiprocessing_native(self):
        matrices = [
            Matrix(matrix_id=1, matrix="2#0,1,0"),
            Matrix(matrix_id=2, matrix="2#3,3/2,4"),
            Matrix(matrix_id=3, matrix="3#4,13/2,1/2,5,11/2,3"),
        ]

        results = list(
            run_multiprocessing(
                "safe",
                matrices,
                config=RunConfig(include_candidates=False),
                mp_config=MPConfig(workers=2),
                run_id="native_mp",
            )
        )

        self.assertEqual(len(results), 3)
        self.assertCountEqual([result["matrix_id"] for result in results], [1, 2, 3])
        self.assertTrue(all(result["status"] == 0 for result in results))

    def test_logging_calls_from_threads_complete_safely(self):
        native = load_native_module()
        original_directory = Path.cwd()

        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "log").mkdir()
            try:
                os.chdir(tmpdir)
                with ThreadPoolExecutor(max_workers=4) as executor:
                    results = list(
                        executor.map(
                            lambda matrix_id: native.compute_matrix(
                                "safe",
                                "2#0,1,0",
                                enable_logging=True,
                                matrix_id=matrix_id,
                            ),
                            range(8),
                        )
                    )
            finally:
                os.chdir(original_directory)

            self.assertTrue((Path(tmpdir) / "log" / "fracessa.log").exists())

        self.assertTrue(all(result["status"] == 0 for result in results))
