import unittest

from wrapper_v1.types import MPConfig, MatrixJob, RunConfig, StatusCode


class TypesTests(unittest.TestCase):
    def test_run_config_has_no_timeout_field(self):
        cfg = RunConfig()
        self.assertFalse(hasattr(cfg, "timeout_s"))

    def test_status_codes(self):
        self.assertEqual(int(StatusCode.OK), 0)
        self.assertEqual(int(StatusCode.PARSE_ERROR), 1)
        self.assertEqual(int(StatusCode.EXEC_ERROR), 4)
        self.assertEqual(int(StatusCode.INTERNAL_ERROR), 255)

    def test_mp_config_validates_values(self):
        with self.assertRaises(ValueError):
            MPConfig(workers=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, prefetch_per_worker=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, queue_maxsize=0)

    def test_matrix_job_validates_its_public_contract(self):
        for matrix_id in (-(1 << 63), (1 << 63) - 1):
            self.assertEqual(MatrixJob(matrix_id, "1#0").matrix_id, matrix_id)

        for matrix_id in (True, 1.0, "1"):
            with self.subTest(matrix_id=matrix_id):
                with self.assertRaisesRegex(TypeError, "matrix_id must be an int"):
                    MatrixJob(matrix_id, "1#0")

        for matrix_id in (-(1 << 63) - 1, 1 << 63):
            with self.subTest(matrix_id=matrix_id):
                with self.assertRaisesRegex(ValueError, "signed 64-bit"):
                    MatrixJob(matrix_id, "1#0")

        with self.assertRaisesRegex(TypeError, "matrix must be a str"):
            MatrixJob(1, 0)
        with self.assertRaisesRegex(TypeError, "metadata must be a dict or None"):
            MatrixJob(1, "1#0", [])
