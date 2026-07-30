import unittest
from unittest import mock

from wrapper_v1 import single
from wrapper_v1.types import MatrixJob


class _TestSink:
    def __init__(self):
        self.written = 0
        self.closed = False

    def write_result(self, result):
        self.written += 1

    def close(self):
        self.closed = True


class _FailingSink:
    def __init__(self):
        self.close_called = False

    def write_result(self, result):
        raise RuntimeError("write failure")

    def close(self):
        self.close_called = True
        raise RuntimeError("close failure")


def _fake_result(matrix_id: int) -> dict:
    return {
        "run_id": "unit",
        "matrix_id": matrix_id,
        "status": 0,
        "success": True,
        "ess_count": 1,
        "elapsed_ns": 1,
        "candidate_count": 0,
        "error_message": "",
        "candidates": [],
        "metadata": None,
    }


class SingleUnitTests(unittest.TestCase):
    def test_run_one(self):
        job = MatrixJob(matrix_id=5, matrix="2#0,1,0")
        with mock.patch("wrapper_v1.single.compute_job", return_value=_fake_result(5)):
            result = single.run_one(job)

        self.assertEqual(result["matrix_id"], 5)

    def test_run_many_and_run_many_to_sink(self):
        jobs = [MatrixJob(matrix_id=1, matrix="2#0,1,0"), MatrixJob(matrix_id=2, matrix="2#3,3/2,4")]

        def _compute(job, config, run_id):
            return _fake_result(job.matrix_id)

        with mock.patch("wrapper_v1.single.compute_job", side_effect=_compute):
            results = list(single.run_many(jobs))

        self.assertEqual([result["matrix_id"] for result in results], [1, 2])

        sink = _TestSink()
        with mock.patch("wrapper_v1.single.compute_job", side_effect=_compute):
            written = single.run_many_to_sink(jobs, sink)

        self.assertEqual(written, 2)
        self.assertEqual(sink.written, 2)
        self.assertTrue(sink.closed)

    def test_run_many_to_sink_preserves_write_error_when_cleanup_fails(self):
        sink = _FailingSink()
        job = MatrixJob(matrix_id=1, matrix="2#0,1,0")

        with mock.patch("wrapper_v1.single.compute_job", return_value=_fake_result(1)):
            with self.assertRaisesRegex(RuntimeError, "write failure"):
                single.run_many_to_sink([job], sink)

        self.assertTrue(sink.close_called)
