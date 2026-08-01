import unittest
from unittest import mock

from fracessa import core
from fracessa.types import Matrix, RunConfig


class _FakeNative:
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
    def test_compute_matrix_uses_cli_string_if_already_prefixed(self):
        fake = _FakeNative()
        matrix = Matrix(matrix_id=11, matrix="2#0,1,0")
        cfg = RunConfig(include_candidates=True)

        with mock.patch("fracessa.core.load_native_module", return_value=fake):
            result = core.compute_matrix(method="fast", matrix=matrix, config=cfg, run_id="unit")

        self.assertEqual(fake.last_kwargs["matrix"], "2#0,1,0")
        self.assertEqual(fake.last_kwargs["method"], "fast")
        self.assertEqual(result["matrix_id"], 11)
        self.assertEqual(result["ess_count"], 2)
        self.assertEqual(result["elapsed_ns"], 1234)
        self.assertEqual(result["candidate_count"], 1)
        self.assertEqual(result["candidates"][0]["candidate_id"], 7)

    def test_compute_matrix_adds_dimension_prefix_from_metadata(self):
        fake = _FakeNative()
        matrix = Matrix(matrix_id=12, matrix="0,1,0", metadata={"dimension": 2})
        cfg = RunConfig(include_candidates=False)

        with mock.patch("fracessa.core.load_native_module", return_value=fake):
            result = core.compute_matrix(method="safe", matrix=matrix, config=cfg, run_id="unit")

        self.assertEqual(fake.last_kwargs["matrix"], "2#0,1,0")
        self.assertEqual(result["candidate_count"], 0)

    def test_compute_matrix_values_only_without_dimension_fails(self):
        fake = _FakeNative()
        matrix = Matrix(matrix_id=14, matrix="0,1,0")
        cfg = RunConfig()

        with mock.patch("fracessa.core.load_native_module", return_value=fake):
            with self.assertRaises(ValueError):
                core.compute_matrix(method="safe", matrix=matrix, config=cfg, run_id="unit")

    def test_compute_matrix_rejects_non_integer_metadata_dimension(self):
        fake = _FakeNative()
        for dimension in (True, 2.9, "2"):
            with self.subTest(dimension=dimension):
                matrix = Matrix(
                    matrix_id=14,
                    matrix="0,1,0",
                    metadata={"dimension": dimension},
                )
                with mock.patch("fracessa.core.load_native_module", return_value=fake):
                    with self.assertRaisesRegex(TypeError, "dimension.*must be an int"):
                        core.compute_matrix(method="safe", matrix=matrix, config=RunConfig(), run_id="unit")
