from __future__ import annotations

from collections.abc import Iterable, Iterator
import multiprocessing as mp
import time

from .core import compute_job, new_run_id
from .types import MPConfig, MatrixJob, MatrixResult, RunConfig, StatusCode, SummaryRow

_SENTINEL = None


def _safe_compute(job: MatrixJob, config: RunConfig, run_id: str) -> MatrixResult:
    try:
        return compute_job(job=job, config=config, run_id=run_id)
    except Exception as exc:  # defensive: worker must never crash the pool protocol
        summary = SummaryRow(
            run_id=run_id,
            matrix_id=int(job.matrix_id),
            status=int(StatusCode.INTERNAL_ERROR),
            success=False,
            ess_count=0,
            elapsed_us=0,
            candidate_count=0,
            error_message=f"Worker exception: {exc}",
        )
        return MatrixResult(summary=summary, candidates=[], metadata=job.metadata)


def _queue_worker(
    input_queue,
    output_queue,
    config: RunConfig,
    run_id: str,
):
    while True:
        item = input_queue.get()
        if item is _SENTINEL:
            return

        seq, job = item
        result = _safe_compute(job=job, config=config, run_id=run_id)
        output_queue.put((seq, result))


def _max_buffered_results(mp_config: MPConfig) -> int:
    """
    Bound submitted-but-not-yielded jobs.

    The batch API may process millions of matrices. Submitting all jobs before
    draining results can deadlock bounded queues, so keep only a small window
    ahead of the consumer. `chunk_size` acts as per-worker prefetch.
    """
    worker_window = mp_config.workers * mp_config.chunk_size
    return max(1, min(mp_config.queue_maxsize, worker_window))


class MPQueueRunner:
    """
    Multiprocessing queue runner for on-the-fly matrix generation.

    Workflow:
    1) `submit(job)` zero or many times,
    2) `close_input()` when done submitting,
    3) iterate `iter_results(expected_results=...)`,
    4) `shutdown()`.
    """

    def __init__(self, run_config: RunConfig, mp_config: MPConfig, run_id: str | None = None):
        self.run_config = run_config
        self.mp_config = mp_config
        self.run_id = run_id or new_run_id("mp")

        self._ctx = mp.get_context(mp_config.start_method)
        self._input_queue = self._ctx.Queue(maxsize=mp_config.queue_maxsize)
        self._output_queue = self._ctx.Queue(maxsize=mp_config.queue_maxsize)
        self._workers: list[mp.Process] = []

        self._submitted = 0
        self._input_closed = False

        for worker_idx in range(mp_config.workers):
            proc = self._ctx.Process(
                target=_queue_worker,
                name=f"fracessa-worker-{worker_idx}",
                args=(self._input_queue, self._output_queue, self.run_config, self.run_id),
                daemon=True,
            )
            proc.start()
            self._workers.append(proc)

    @property
    def submitted(self) -> int:
        return self._submitted

    def submit(self, job: MatrixJob) -> int:
        if self._input_closed:
            raise RuntimeError("Cannot submit after close_input()")

        seq = self._submitted
        self._input_queue.put((seq, job))
        self._submitted += 1
        return seq

    def close_input(self) -> None:
        if self._input_closed:
            return

        for _ in self._workers:
            self._input_queue.put(_SENTINEL)
        self._input_closed = True

    def iter_results(self, expected_results: int | None = None) -> Iterator[MatrixResult]:
        expected = self._submitted if expected_results is None else expected_results
        if expected < 0:
            raise ValueError("expected_results must be >= 0")

        if self.mp_config.ordered_results:
            pending: dict[int, MatrixResult] = {}
            next_seq = 0

            for _ in range(expected):
                seq, result = self._get_result_item()
                pending[int(seq)] = result

                while next_seq in pending:
                    yield pending.pop(next_seq)
                    next_seq += 1
        else:
            for _ in range(expected):
                _, result = self._get_result_item()
                yield result

    def _get_result_item(self):
        return self._output_queue.get()

    def shutdown(self, join_timeout_s: float = 5.0) -> None:
        if not self._input_closed:
            self.close_input()

        deadline = time.time() + join_timeout_s
        for proc in self._workers:
            timeout = max(0.0, deadline - time.time())
            proc.join(timeout=timeout)
            if proc.is_alive():
                proc.terminate()
                proc.join(timeout=1.0)

        self._input_queue.close()
        self._output_queue.close()


def run_jobs_mp(
    jobs: Iterable[MatrixJob],
    run_config: RunConfig,
    mp_config: MPConfig,
    run_id: str | None = None,
) -> Iterator[MatrixResult]:
    runner = MPQueueRunner(run_config=run_config, mp_config=mp_config, run_id=run_id)
    try:
        jobs_iter = iter(jobs)
        max_buffered = _max_buffered_results(mp_config)

        submitted = 0
        yielded = 0
        next_seq = 0
        pending: dict[int, MatrixResult] = {}
        input_exhausted = False
        input_closed = False

        def submit_until_window() -> None:
            nonlocal submitted, input_exhausted, input_closed
            while not input_exhausted and submitted - yielded < max_buffered:
                try:
                    job = next(jobs_iter)
                except StopIteration:
                    input_exhausted = True
                    if not input_closed:
                        runner.close_input()
                        input_closed = True
                    return

                runner.submit(job)
                submitted += 1

        submit_until_window()

        while yielded < submitted or not input_exhausted:
            if yielded >= submitted:
                submit_until_window()
                if yielded >= submitted and input_exhausted:
                    break

            seq, result = runner._get_result_item()
            seq = int(seq)

            if mp_config.ordered_results:
                pending[seq] = result
                while next_seq in pending:
                    yield pending.pop(next_seq)
                    yielded += 1
                    next_seq += 1
                    submit_until_window()
            else:
                yield result
                yielded += 1
                submit_until_window()

        if not input_closed:
            runner.close_input()
    finally:
        runner.shutdown()


def run_jobs_mp_to_sink(
    jobs: Iterable[MatrixJob],
    sink,
    run_config: RunConfig,
    mp_config: MPConfig,
    run_id: str | None = None,
) -> int:
    count = 0
    try:
        for result in run_jobs_mp(jobs=jobs, run_config=run_config, mp_config=mp_config, run_id=run_id):
            sink.write_result(result)
            count += 1
    finally:
        sink.close()
    return count
