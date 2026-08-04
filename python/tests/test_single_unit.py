import unittest
from unittest import mock

from pyfracessa import single
from pyfracessa.types import Matrix


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
        "ess_count": 1,
        "elapsed_ns": 1,
        "safe_fallback": None,
        "candidate_count": 0,
        "error_message": "",
        "candidates": [],
        "metadata": None,
    }


class SingleUnitTests(unittest.TestCase):
    def test_run_accepts_one_matrix(self):
        matrix = Matrix(matrix_id=5, matrix="2#0,1,0")
        with mock.patch("pyfracessa.single.compute_matrix", return_value=_fake_result(5)):
            result = single.run("safe", matrix)

        self.assertEqual(result["matrix_id"], 5)

    def test_run_accepts_many_matrices_and_a_sink(self):
        matrices = [Matrix(matrix_id=1, matrix="2#0,1,0"), Matrix(matrix_id=2, matrix="2#3,3/2,4")]
        run_ids = []

        def _compute(method, matrix, config, run_id):
            self.assertEqual(method, "safe")
            run_ids.append(run_id)
            return _fake_result(matrix.matrix_id)

        with mock.patch("pyfracessa.single.compute_matrix", side_effect=_compute):
            results = list(single.run("safe", matrices, run_id="shared"))

        self.assertEqual([result["matrix_id"] for result in results], [1, 2])
        self.assertEqual(run_ids, ["shared", "shared"])

        sink = _TestSink()
        run_ids.clear()
        with mock.patch("pyfracessa.single.compute_matrix", side_effect=_compute):
            written = single.run("safe", matrices, sink=sink)

        self.assertEqual(written, 2)
        self.assertEqual(sink.written, 2)
        self.assertTrue(sink.closed)

    def test_run_with_sink_preserves_write_error_when_cleanup_fails(self):
        sink = _FailingSink()
        matrix = Matrix(matrix_id=1, matrix="2#0,1,0")

        with mock.patch("pyfracessa.single.compute_matrix", return_value=_fake_result(1)):
            with self.assertRaisesRegex(RuntimeError, "write failure"):
                single.run("safe", matrix, sink=sink)

        self.assertTrue(sink.close_called)
