from queue import Empty
from types import SimpleNamespace
import unittest
from unittest import mock

from fracessa import mp as mp_mod
from fracessa.types import MPConfig, Matrix, RunConfig


def _fake_result(matrix_id: int) -> dict:
    return {
        "run_id": "unit",
        "matrix_id": matrix_id,
        "status": 0,
        "success": True,
        "ess_count": 1,
        "elapsed_ns": 1,
        "safe_fallback": None,
        "candidate_count": 0,
        "error_message": "",
        "candidates": [],
        "metadata": None,
    }


class _FakeRunner:
    instances = []

    def __init__(self, method, run_config, mp_config, run_id=None):
        self.method = method
        self.run_config = run_config
        self.mp_config = mp_config
        self.run_id = run_id
        self._submitted = 0
        self._received = 0
        self._ready = []
        self._closed = False
        self._shutdown = False
        self.cancelled = None
        self.max_submitted_not_received = 0
        self.__class__.instances.append(self)

    def submit(self, matrix):
        if self._closed:
            raise AssertionError("submit called after close_input")

        self._submitted += 1
        self._ready.append(matrix)
        outstanding = self._submitted - self._received
        self.max_submitted_not_received = max(self.max_submitted_not_received, outstanding)
        if outstanding > self.mp_config.queue_maxsize:
            raise AssertionError("batch runner submitted too far ahead of result draining")

    def close_input(self):
        self._closed = True

    def get_result(self):
        if not self._ready:
            raise AssertionError("result requested before a matrix was submitted")

        matrix = self._ready.pop()
        self._received += 1
        return _fake_result(matrix.matrix_id)

    def shutdown(self, *, cancel):
        self._shutdown = True
        self.cancelled = cancel


class MPUnitTests(unittest.TestCase):
    def test_run_multiprocessing_accepts_one_matrix(self):
        matrix = Matrix(matrix_id=7, matrix="2#0,1,0")

        _FakeRunner.instances.clear()
        with mock.patch("fracessa.mp._QueueRunner", _FakeRunner):
            result = mp_mod.run_multiprocessing(
                "safe",
                matrix,
                run_id="unit_mp",
            )

        self.assertEqual(result["matrix_id"], 7)
        self.assertEqual(_FakeRunner.instances[-1].mp_config, MPConfig())

    def test_run_multiprocessing_accepts_a_sink(self):
        matrices = [Matrix(matrix_id=i, matrix="2#0,1,0") for i in range(2)]
        sink = mock.Mock()

        _FakeRunner.instances.clear()
        with mock.patch("fracessa.mp._QueueRunner", _FakeRunner):
            count = mp_mod.run_multiprocessing(
                "safe",
                matrices,
                sink=sink,
                mp_config=MPConfig(workers=1),
            )

        self.assertEqual(count, 2)
        self.assertEqual(sink.write_result.call_count, 2)
        sink.close.assert_called_once_with()

    def test_run_multiprocessing_yields_as_completed_and_bounds_submission(self):
        matrices = [Matrix(matrix_id=i, matrix="2#0,1,0") for i in range(10)]
        cfg = MPConfig(workers=1, prefetch_per_worker=10, queue_maxsize=2)

        _FakeRunner.instances.clear()
        with mock.patch("fracessa.mp._QueueRunner", _FakeRunner):
            results = list(
                mp_mod.run_multiprocessing(
                    "safe",
                    matrices=matrices,
                    config=RunConfig(),
                    mp_config=cfg,
                    run_id="unit_mp",
                )
            )

        runner = _FakeRunner.instances[-1]
        self.assertEqual([result["matrix_id"] for result in results], [*range(1, 10), 0])
        self.assertLessEqual(runner.max_submitted_not_received, cfg.queue_maxsize)
        self.assertTrue(runner._closed)
        self.assertTrue(runner._shutdown)
        self.assertFalse(runner.cancelled)

    def test_run_multiprocessing_cancels_when_consumer_stops_early(self):
        matrices = [Matrix(matrix_id=i, matrix="2#0,1,0") for i in range(3)]
        cfg = MPConfig(workers=1, prefetch_per_worker=1, queue_maxsize=1)

        _FakeRunner.instances.clear()
        with mock.patch("fracessa.mp._QueueRunner", _FakeRunner):
            results = mp_mod.run_multiprocessing("safe", matrices, mp_config=cfg)
            next(results)
            results.close()

        runner = _FakeRunner.instances[-1]
        self.assertTrue(runner._shutdown)
        self.assertTrue(runner.cancelled)

    def test_run_multiprocessing_rejects_logging_before_starting_workers(self):
        matrices = [Matrix(matrix_id=1, matrix="2#0,1,0")]

        with mock.patch("fracessa.mp._QueueRunner") as runner:
            with self.assertRaisesRegex(ValueError, "enable_logging"):
                list(
                    mp_mod.run_multiprocessing(
                        "safe",
                        matrices,
                        config=RunConfig(enable_logging=True),
                        mp_config=MPConfig(workers=1),
                    )
                )

        runner.assert_not_called()

    def test_submit_serializes_before_queueing(self):
        runner = object.__new__(mp_mod._QueueRunner)
        runner._input_closed = False
        runner._input_queue = mock.Mock()
        matrix = Matrix(matrix_id=1, matrix="2#0,1,0", metadata={"bad": lambda: None})

        with self.assertRaises(Exception):
            runner.submit(matrix)

        runner._input_queue.put.assert_not_called()

    def test_worker_serializes_results_before_queueing(self):
        input_queue = mock.Mock()
        input_queue.get.side_effect = [
            bytes(mp_mod.ForkingPickler.dumps(Matrix(matrix_id=1, matrix="2#0,1,0"))),
            None,
        ]
        output_queue = mock.Mock()

        with mock.patch("fracessa.mp._safe_compute", return_value={"bad": lambda: None}):
            with self.assertRaises(Exception):
                mp_mod._queue_worker(input_queue, output_queue, "safe", RunConfig(), "unit")

        output_queue.put.assert_not_called()

    def test_safe_compute_does_not_reconvert_an_invalid_mutated_id(self):
        matrix = Matrix(matrix_id=1, matrix="1#0")
        matrix.matrix_id = "invalid"

        with mock.patch("fracessa.mp.compute_matrix", side_effect=RuntimeError("forced")):
            result = mp_mod._safe_compute("safe", matrix, RunConfig(), "unit")

        self.assertEqual(result["matrix_id"], -1)
        self.assertEqual(result["status"], 255)
        self.assertIn("forced", result["error_message"])

    def test_get_result_reports_dead_worker(self):
        runner = object.__new__(mp_mod._QueueRunner)
        runner._output_queue = mock.Mock()
        runner._output_queue.get.side_effect = Empty
        runner._workers = [SimpleNamespace(name="fracessa-worker-0", exitcode=7)]

        with self.assertRaisesRegex(RuntimeError, "fracessa-worker-0=7"):
            runner.get_result()

    def test_max_pending_matrices_uses_worker_prefetch_and_queue_cap(self):
        self.assertEqual(
            mp_mod._max_pending_matrices(
                MPConfig(workers=4, prefetch_per_worker=3, queue_maxsize=20)
            ),
            12,
        )
        self.assertEqual(
            mp_mod._max_pending_matrices(
                MPConfig(workers=4, prefetch_per_worker=3, queue_maxsize=5)
            ),
            5,
        )
