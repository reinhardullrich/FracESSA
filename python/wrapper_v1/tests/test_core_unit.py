import unittest
from unittest import mock

from wrapper_v1 import core
from wrapper_v1.types import MatrixJob, RunConfig


class _FakeNative:
    STATUS_OK = 0
    STATUS_PARSE_ERROR = 1
    STATUS_EXEC_ERROR = 4
    STATUS_INTERNAL_ERROR = 255

    def __init__(self):
        self.last_kwargs = None

    def compute_matrix(self, **kwargs):
        self.last_kwargs = kwargs
        candidates = []
        if kwargs.get("include_candidates", False):
            candidates = [
                {
                    "candidate_id": 7,
                    "vector": "1/2,1/2",
                    "support": 3,
                    "support_size": 2,
                    "extended_support": 3,
                    "extended_support_size": 2,
                    "multiplier": None,
                    "is_ess": True,
                    "stability": "pure",
                    "payoff": "1",
                    "payoff_dbl": 1.0,
                }
            ]
        return {
            "status": 0,
            "success": True,
            "error_message": "",
            "ess_count": 2,
            "elapsed_ns": 1234,
            "candidates": candidates,
        }


class CoreUnitTests(unittest.TestCase):
    def test_native_status_map(self):
        fake = _FakeNative()
        with mock.patch("wrapper_v1.core.load_native_module", return_value=fake):
            status_map = core.native_status_map()

        self.assertEqual(
            status_map,
            {
                "OK": 0,
                "PARSE_ERROR": 1,
                "EXEC_ERROR": 4,
                "INTERNAL_ERROR": 255,
            },
        )

    def test_compute_job_uses_cli_string_if_already_prefixed(self):
        fake = _FakeNative()
        job = MatrixJob(matrix_id=11, matrix="2#0,1,0")
        cfg = RunConfig(include_candidates=True)

        with mock.patch("wrapper_v1.core.load_native_module", return_value=fake):
            result = core.compute_job(job=job, config=cfg, run_id="unit")

        self.assertEqual(fake.last_kwargs["matrix"], "2#0,1,0")
        self.assertEqual(result["matrix_id"], 11)
        self.assertEqual(result["ess_count"], 2)
        self.assertEqual(result["elapsed_ns"], 1234)
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(result["candidates"][0]["candidate_id"], 7)

    def test_compute_job_adds_dimension_prefix_from_metadata(self):
        fake = _FakeNative()
        job = MatrixJob(matrix_id=12, matrix="0,1,0", metadata={"dimension": 2})
        cfg = RunConfig(include_candidates=False)

        with mock.patch("wrapper_v1.core.load_native_module", return_value=fake):
            result = core.compute_job(job=job, config=cfg, run_id="unit")

        self.assertEqual(fake.last_kwargs["matrix"], "2#0,1,0")
        self.assertEqual(result["candidate_count"], 0)

    def test_compute_job_values_only_without_dimension_fails(self):
        fake = _FakeNative()
        job = MatrixJob(matrix_id=14, matrix="0,1,0")
        cfg = RunConfig()

        with mock.patch("wrapper_v1.core.load_native_module", return_value=fake):
            with self.assertRaises(ValueError):
                core.compute_job(job=job, config=cfg, run_id="unit")

    def test_compute_job_rejects_non_integer_metadata_dimension(self):
        fake = _FakeNative()
        for dimension in (True, 2.9, "2"):
            with self.subTest(dimension=dimension):
                job = MatrixJob(
                    matrix_id=14,
                    matrix="0,1,0",
                    metadata={"dimension": dimension},
                )
                with mock.patch("wrapper_v1.core.load_native_module", return_value=fake):
                    with self.assertRaisesRegex(TypeError, "dimension.*must be an int"):
                        core.compute_job(job=job, config=RunConfig(), run_id="unit")
