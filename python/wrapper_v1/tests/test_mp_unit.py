import unittest
from unittest import mock

from wrapper_v1 import mp as mp_mod
from wrapper_v1.types import MPConfig, MatrixJob, MatrixResult, RunConfig, SummaryRow


def _fake_result(matrix_id: int) -> MatrixResult:
    return MatrixResult(
        summary=SummaryRow(
            run_id="unit",
            matrix_id=matrix_id,
            status=0,
            success=True,
            ess_count=1,
            elapsed_us=1,
            candidate_count=0,
            error_message="",
        ),
        candidates=[],
        metadata=None,
    )


class _FakeRunner:
    instances = []

    def __init__(self, run_config, mp_config, run_id=None):
        self.run_config = run_config
        self.mp_config = mp_config
        self.run_id = run_id
        self._submitted = 0
        self._received = 0
        self._ready = []
        self._closed = False
        self._shutdown = False
        self.max_submitted_not_received = 0
        self.__class__.instances.append(self)

    def submit(self, job):
        if self._closed:
            raise AssertionError("submit called after close_input")

        seq = self._submitted
        self._submitted += 1
        self._ready.append((seq, job))
        outstanding = self._submitted - self._received
        self.max_submitted_not_received = max(self.max_submitted_not_received, outstanding)
        if outstanding > self.mp_config.queue_maxsize:
            raise AssertionError("batch runner submitted too far ahead of result draining")
        return seq

    def close_input(self):
        self._closed = True

    def _get_result_item(self):
        if not self._ready:
            raise AssertionError("result requested before a job was submitted")

        seq, job = self._ready.pop(0)
        self._received += 1
        return seq, _fake_result(job.matrix_id)

    def shutdown(self):
        self._shutdown = True


class MPUnitTests(unittest.TestCase):
    def test_run_jobs_mp_drains_while_submitting(self):
        jobs = [MatrixJob(matrix_id=i, matrix="2#0,1,0") for i in range(10)]
        cfg = MPConfig(workers=1, chunk_size=10, queue_maxsize=2, ordered_results=True)

        _FakeRunner.instances.clear()
        with mock.patch("wrapper_v1.mp.MPQueueRunner", _FakeRunner):
            results = list(
                mp_mod.run_jobs_mp(
                    jobs=jobs,
                    run_config=RunConfig(),
                    mp_config=cfg,
                    run_id="unit_mp",
                )
            )

        runner = _FakeRunner.instances[-1]
        self.assertEqual([result.summary.matrix_id for result in results], list(range(10)))
        self.assertLessEqual(runner.max_submitted_not_received, cfg.queue_maxsize)
        self.assertTrue(runner._closed)
        self.assertTrue(runner._shutdown)

    def test_max_buffered_results_uses_worker_chunk_and_queue_cap(self):
        self.assertEqual(
            mp_mod._max_buffered_results(MPConfig(workers=4, chunk_size=3, queue_maxsize=20)),
            12,
        )
        self.assertEqual(
            mp_mod._max_buffered_results(MPConfig(workers=4, chunk_size=3, queue_maxsize=5)),
            5,
        )
