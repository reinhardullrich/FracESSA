import inspect
import unittest

import pyfracessa
from pyfracessa.types import MPConfig, Matrix, RunConfig, StatusCode, _validate_search_method


class TypesTests(unittest.TestCase):
    def test_run_config_has_no_timeout_field(self):
        cfg = RunConfig()
        self.assertFalse(hasattr(cfg, "timeout_s"))
        self.assertFalse(hasattr(cfg, "mode"))
        self.assertFalse(hasattr(cfg, "cyclic_symmetry_filter"))

    def test_search_method_is_required_and_validated(self):
        for method in ("fast", "safe"):
            _validate_search_method(method)
        with self.assertRaisesRegex(TypeError, "method must be a str"):
            _validate_search_method(1)
        for method in ("test", "verified", "exact", "unsafe", "unknown"):
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, "fast or safe"):
                _validate_search_method(method)

    def test_status_codes(self):
        self.assertEqual(int(StatusCode.OK), 0)
        self.assertEqual(int(StatusCode.PARSE_ERROR), 1)
        self.assertEqual(int(StatusCode.EXEC_ERROR), 4)
        self.assertEqual(int(StatusCode.INTERNAL_ERROR), 255)

    def test_mp_config_validates_values(self):
        self.assertGreaterEqual(MPConfig().workers, 1)
        with self.assertRaises(ValueError):
            MPConfig(workers=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, prefetch_per_worker=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, queue_maxsize=0)

    def test_matrix_validates_its_public_contract(self):
        for matrix_id in (-(1 << 63), (1 << 63) - 1):
            self.assertEqual(Matrix(matrix_id, "1#0").matrix_id, matrix_id)

        for matrix_id in (True, 1.0, "1"):
            with self.subTest(matrix_id=matrix_id):
                with self.assertRaisesRegex(TypeError, "matrix_id must be an int"):
                    Matrix(matrix_id, "1#0")

        for matrix_id in (-(1 << 63) - 1, 1 << 63):
            with self.subTest(matrix_id=matrix_id):
                with self.assertRaisesRegex(ValueError, "signed 64-bit"):
                    Matrix(matrix_id, "1#0")

        with self.assertRaisesRegex(TypeError, "matrix must be a str"):
            Matrix(1, 0)
        with self.assertRaisesRegex(TypeError, "metadata must be a dict or None"):
            Matrix(1, "1#0", [])

    def test_public_execution_api(self):
        self.assertEqual(
            list(inspect.signature(pyfracessa.compute_matrix).parameters),
            ["method", "matrix", "config", "run_id"],
        )
        self.assertEqual(
            list(inspect.signature(pyfracessa.run).parameters),
            ["method", "matrices", "config", "run_id", "sink"],
        )
        self.assertEqual(
            list(inspect.signature(pyfracessa.run_multiprocessing).parameters),
            ["method", "matrices", "config", "run_id", "sink", "mp_config"],
        )
        self.assertIsNone(
            inspect.signature(pyfracessa.run_multiprocessing).parameters["mp_config"].default
        )
        for function in (pyfracessa.compute_matrix, pyfracessa.run, pyfracessa.run_multiprocessing):
            self.assertIs(inspect.signature(function).parameters["method"].default, inspect.Parameter.empty)
        self.assertEqual(
            pyfracessa.__all__,
            [
                "StatusCode",
                "Matrix",
                "RunConfig",
                "MPConfig",
                "new_run_id",
                "compute_matrix",
                "run",
                "run_multiprocessing",
                "load_matrices_from_json",
                "create_sink",
                "CsvSink",
                "JsonSink",
                "ParquetSink",
            ],
        )
