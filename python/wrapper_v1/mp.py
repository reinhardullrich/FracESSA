from __future__ import annotations

from collections.abc import Iterable, Iterator
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler
import pickle
from queue import Empty
import time

from .core import compute_job, new_run_id
from .sinks import _consume_to_sink
from .types import MPConfig, MatrixJob, RunConfig, StatusCode

_SENTINEL = None


def _safe_compute(job: MatrixJob, config: RunConfig, run_id: str) -> dict:
    try:
        return compute_job(job=job, config=config, run_id=run_id)
    except Exception as exc:  # defensive: worker must never crash the pool protocol
        return {
            "run_id": run_id,
            "matrix_id": job.matrix_id if type(job.matrix_id) is int else -1,
            "status": int(StatusCode.INTERNAL_ERROR),
            "success": False,
            "ess_count": 0,
            "elapsed_ns": 0,
            "candidate_count": 0,
            "error_message": f"Worker exception: {exc}",
            "candidates": [],
            "metadata": job.metadata,
        }


def _queue_worker(
    input_queue,
    output_queue,
    config: RunConfig,
    run_id: str,
):
    while True:
        payload = input_queue.get()
        if payload is _SENTINEL:
            return

        job = pickle.loads(payload)
        result = _safe_compute(job=job, config=config, run_id=run_id)
        output_queue.put(bytes(ForkingPickler.dumps(result)))


def _max_pending_jobs(mp_config: MPConfig) -> int:
    """
    Bound submitted-but-not-yielded jobs.

    The batch API may process millions of matrices. Submitting all jobs before
    draining results can deadlock bounded queues, so keep only a small window
    ahead of the consumer.
    """
    worker_window = mp_config.workers * mp_config.prefetch_per_worker
    return max(1, min(mp_config.queue_maxsize, worker_window))


class _QueueRunner:
    def __init__(self, run_config: RunConfig, mp_config: MPConfig, run_id: str | None = None):
        self.run_config = run_config
        self.mp_config = mp_config
        self.run_id = run_id or new_run_id("mp")

        self._ctx = mp.get_context(mp_config.start_method)
        self._input_queue = self._ctx.Queue()
        self._output_queue = self._ctx.Queue(maxsize=mp_config.queue_maxsize)
        self._workers: list[mp.Process] = []

        self._input_closed = False

        try:
            for worker_idx in range(mp_config.workers):
                proc = self._ctx.Process(
                    target=_queue_worker,
                    name=f"fracessa-worker-{worker_idx}",
                    args=(self._input_queue, self._output_queue, self.run_config, self.run_id),
                    daemon=True,
                )
                proc.start()
                self._workers.append(proc)
        except BaseException:
            self.shutdown(cancel=True)
            raise

    def submit(self, job: MatrixJob) -> None:
        if self._input_closed:
            raise RuntimeError("Cannot submit after close_input()")

        payload = bytes(ForkingPickler.dumps(job))
        self._input_queue.put(payload)

    def close_input(self) -> None:
        if self._input_closed:
            return

        for _ in self._workers:
            self._input_queue.put(_SENTINEL)
        self._input_closed = True

    def get_result(self):
        while True:
            try:
                return pickle.loads(self._output_queue.get(timeout=0.1))
            except Empty:
                failed = [proc for proc in self._workers if proc.exitcode not in (None, 0)]
                if failed:
                    details = ", ".join(f"{proc.name}={proc.exitcode}" for proc in failed)
                    raise RuntimeError(f"FracESSA worker exited without a result: {details}")
                if all(proc.exitcode is not None for proc in self._workers):
                    raise RuntimeError("All FracESSA workers exited before producing the expected result")

    def shutdown(self, *, cancel: bool, join_timeout_s: float = 5.0) -> None:
        if cancel:
            self._input_queue.cancel_join_thread()
            for proc in self._workers:
                if proc.is_alive():
                    proc.terminate()
        elif not self._input_closed:
            self.close_input()

        deadline = time.monotonic() + join_timeout_s
        for proc in self._workers:
            timeout = max(0.0, deadline - time.monotonic())
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
) -> Iterator[dict]:
    if run_config.enable_logging:
        raise ValueError("RunConfig.enable_logging is not supported with multiprocessing")

    runner = _QueueRunner(run_config=run_config, mp_config=mp_config, run_id=run_id)
    completed = False
    try:
        jobs_iter = iter(jobs)
        max_pending = _max_pending_jobs(mp_config)

        submitted = 0
        yielded = 0
        input_exhausted = False

        def submit_until_window() -> None:
            nonlocal submitted, input_exhausted
            while not input_exhausted and submitted - yielded < max_pending:
                try:
                    job = next(jobs_iter)
                except StopIteration:
                    input_exhausted = True
                    runner.close_input()
                    return

                runner.submit(job)
                submitted += 1

        submit_until_window()

        while yielded < submitted or not input_exhausted:
            if yielded >= submitted:
                submit_until_window()
                if yielded >= submitted and input_exhausted:
                    break

            yield runner.get_result()
            yielded += 1
            submit_until_window()

        completed = True
    finally:
        runner.shutdown(cancel=not completed)


def run_jobs_mp_to_sink(
    jobs: Iterable[MatrixJob],
    sink,
    run_config: RunConfig,
    mp_config: MPConfig,
    run_id: str | None = None,
) -> int:
    return _consume_to_sink(
        run_jobs_mp(
            jobs=jobs,
            run_config=run_config,
            mp_config=mp_config,
            run_id=run_id,
        ),
        sink,
    )
