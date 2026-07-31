import inspect
import unittest

import fracessa
from fracessa.types import MPConfig, Matrix, RunConfig, StatusCode


class TypesTests(unittest.TestCase):
    def test_run_config_has_no_timeout_field(self):
        cfg = RunConfig()
        self.assertFalse(hasattr(cfg, "timeout_s"))
        self.assertEqual(cfg.mode, "verified")

    def test_run_config_validates_mode(self):
        for mode in ("verified", "exact", "unsafe", "very_unsafe"):
            self.assertEqual(RunConfig(mode=mode).mode, mode)
        with self.assertRaisesRegex(TypeError, "mode must be a str"):
            RunConfig(mode=1)
        with self.assertRaisesRegex(ValueError, "verified, exact, unsafe, or very_unsafe"):
            RunConfig(mode="unknown")

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
            list(inspect.signature(fracessa.run).parameters),
            ["matrices", "config", "run_id", "sink"],
        )
        self.assertEqual(
            list(inspect.signature(fracessa.run_multiprocessing).parameters),
            ["matrices", "config", "run_id", "sink", "mp_config"],
        )
        self.assertIsNone(
            inspect.signature(fracessa.run_multiprocessing).parameters["mp_config"].default
        )
        self.assertEqual(
            fracessa.__all__,
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
