from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
import os
from pathlib import Path
import tempfile
import unittest

from wrapper_v1 import (
    MPConfig,
    MatrixJob,
    RunConfig,
    load_native_module,
    run_jobs_mp,
    run_one,
)


class NativeIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_native_module()

    def test_run_one_native(self):
        matrix_id = (1 << 63) - 1
        result = run_one(
            MatrixJob(matrix_id=matrix_id, matrix="2#0,1,0"),
            RunConfig(include_candidates=True),
            run_id="native_single",
        )

        self.assertTrue(result["success"])
        self.assertEqual(result["status"], 0)
        self.assertEqual(result["ess_count"], 1)
        self.assertEqual(result["matrix_id"], matrix_id)
        self.assertEqual(result["candidates"][0]["matrix_id"], matrix_id)
        self.assertIsNone(result["candidates"][0]["multiplier"])

    def test_native_candidate_contract(self):
        native = load_native_module()
        result = native.compute_matrix("2#0,1,0", include_candidates=True)
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
                "stability": "T_pd_frc",
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
        result = run_one(
            MatrixJob(matrix_id=2, matrix="5#1,3"),
            RunConfig(include_candidates=True, exact=True),
            run_id="native_circular",
        )

        self.assertEqual(result["ess_count"], 5)
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(len(result["candidates"]), 1)
        self.assertEqual(result["candidates"][0]["multiplier"], 5)

    def test_invalid_matrix_strings_return_parser_error(self):
        invalid_matrices = {
            "2##0": "Multiple '#' characters found in matrix string",
            "64#0": "Parser supports dimensions in [1, 63], got 64",
            "2#0,1": "Expected 1 (CS) or 3 (Sym) values, got 2",
            "2#0,1/0,0": "Rational denominator cannot be zero: 1/0",
        }

        for matrix, error_message in invalid_matrices.items():
            with self.subTest(matrix=matrix):
                result = run_one(MatrixJob(matrix_id=1, matrix=matrix), run_id="invalid")
                self.assertFalse(result["success"])
                self.assertEqual(result["status"], 1)
                self.assertEqual(result["error_message"], error_message)

    def test_run_jobs_mp_native(self):
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
        jobs = [
            MatrixJob(matrix_id=1, matrix="2#0,1,0"),
            MatrixJob(matrix_id=2, matrix="2#3,3/2,4"),
            MatrixJob(matrix_id=3, matrix="3#4,13/2,1/2,5,11/2,3"),
        ]

        results = list(
            run_jobs_mp(
                jobs,
                run_config=RunConfig(include_candidates=False),
                mp_config=MPConfig(workers=2, start_method=start_method),
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
