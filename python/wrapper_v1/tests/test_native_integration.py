import multiprocessing as mp
import unittest

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
            self.skipTest("fracessa_core not available; run ./fracessa/build.sh first")

    def test_run_one_native(self):
        result = run_one(
            MatrixJob(matrix_id=1, matrix="2#0,1,0"),
            RunConfig(include_candidates=True),
            run_id="native_single",
        )
        self.assertTrue(result.summary.success)
        self.assertEqual(result.summary.status, 0)
        self.assertEqual(result.summary.ess_count, 1)

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
