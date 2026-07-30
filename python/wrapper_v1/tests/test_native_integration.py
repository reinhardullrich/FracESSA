import json
import multiprocessing as mp
import unittest
from pathlib import Path

from wrapper_v1 import (
    MPConfig,
    MPQueueRunner,
    MatrixJob,
    RunConfig,
    load_native_module,
    run_jobs_mp,
    run_one,
)

class NativeIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            load_native_module()
            cls.native_available = True
        except ModuleNotFoundError:
            cls.native_available = False

    def setUp(self):
        if not self.native_available:
            self.skipTest("fracessa_core not available; run ./build.sh first")

    def test_bounded_and_unsafe_routes_native(self):
        verification_file = Path(__file__).resolve().parents[2] / "verification/verification_matrices.json"
        with verification_file.open("r", encoding="utf-8") as fh:
            matrix = next(matrix for matrix in json.load(fh)["matrices"] if matrix["id"] == 46)
        job = MatrixJob(matrix_id=46, matrix=f"{matrix['dimension']}#{matrix['matrix']}")

        bounded = run_one(
            job,
            RunConfig(include_candidates=True),
            run_id="native_bounded",
        )
        unsafe = run_one(
            job,
            RunConfig(include_candidates=True, unsafe=True),
            run_id="native_unsafe",
        )

        self.assertTrue(bounded.summary.success)
        self.assertTrue(unsafe.summary.success)
        self.assertEqual(bounded.summary.ess_count, 1)
        self.assertEqual(unsafe.summary.ess_count, 0)
        self.assertIsNone(bounded.candidates[0].multiplier)
        self.assertEqual(unsafe.candidates, [])

    def test_exact_and_unsafe_are_rejected(self):
        result = run_one(
            MatrixJob(matrix_id=1, matrix="2#0,1,0"),
            RunConfig(exact=True, unsafe=True),
            run_id="native_mode_conflict",
        )
        self.assertFalse(result.summary.success)
        self.assertEqual(result.summary.status, 4)

    def test_circular_native_returns_one_weighted_representative(self):
        result = run_one(
            MatrixJob(matrix_id=2, matrix="5#1,3"),
            RunConfig(include_candidates=True, exact=True),
            run_id="native_circular",
        )

        self.assertEqual(result.summary.ess_count, 5)
        self.assertEqual(result.summary.candidate_count, 1)
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(result.candidates[0].multiplier, 5)

    def test_run_jobs_mp_native(self):
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
        jobs = [
            MatrixJob(matrix_id=1, matrix="2#0,1,0"),
            MatrixJob(matrix_id=2, matrix="2#3,3/2,4"),
            MatrixJob(matrix_id=3, matrix="3#4,13/2,1/2,5,11/2,3"),
        ]

        try:
            results = list(
                run_jobs_mp(
                    jobs,
                    run_config=RunConfig(include_candidates=False),
                    mp_config=MPConfig(workers=2, ordered_results=True, start_method=start_method),
                    run_id="native_mp",
                )
            )
        except PermissionError as exc:
            self.skipTest(f"multiprocessing primitives blocked in this environment: {exc}")

        self.assertEqual(len(results), 3)
        self.assertTrue(all(r.summary.status == 0 for r in results))

    def test_queue_runner_native(self):
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"

        try:
            runner = MPQueueRunner(
                run_config=RunConfig(include_candidates=False),
                mp_config=MPConfig(workers=2, ordered_results=False, start_method=start_method),
                run_id="native_queue",
            )
        except PermissionError as exc:
            self.skipTest(f"multiprocessing primitives blocked in this environment: {exc}")

        try:
            for i in range(3):
                runner.submit(MatrixJob(matrix_id=100 + i, matrix="2#0,1,0"))

            runner.close_input()
            results = list(runner.iter_results(expected_results=runner.submitted))
        finally:
            runner.shutdown()

        self.assertEqual(len(results), 3)
        self.assertTrue(all(r.summary.status == 0 for r in results))
