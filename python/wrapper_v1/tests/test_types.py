import unittest

from wrapper_v1.types import MPConfig, RunConfig, StatusCode


class TypesTests(unittest.TestCase):
    def test_run_config_has_no_timeout_field(self):
        cfg = RunConfig()
        self.assertFalse(hasattr(cfg, "timeout_s"))

    def test_status_codes(self):
        self.assertEqual(int(StatusCode.OK), 0)
        self.assertEqual(int(StatusCode.PARSE_ERROR), 1)
        self.assertEqual(int(StatusCode.DIMENSION_OUT_OF_RANGE), 2)
        self.assertEqual(int(StatusCode.INVALID_VALUE_COUNT), 3)
        self.assertEqual(int(StatusCode.EXEC_ERROR), 4)
        self.assertEqual(int(StatusCode.INTERNAL_ERROR), 255)

    def test_mp_config_validates_values(self):
        with self.assertRaises(ValueError):
            MPConfig(workers=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, chunk_size=0)
        with self.assertRaises(ValueError):
            MPConfig(workers=1, queue_maxsize=0)
