"""Run PyFracESSA matrices in processes with bounded shared queues."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler
import pickle
from queue import Empty
import time

from .core import compute_matrix, new_run_id
from .sinks import _consume_to_sink
from .types import MPConfig, Matrix, RunConfig, SearchMethod, StatusCode, _validate_search_method

_SENTINEL = None


def _safe_compute(method: SearchMethod, matrix: Matrix, config: RunConfig, run_id: str) -> dict:
    """Compute a matrix or convert an unexpected worker error to a result row."""

    try:
        return compute_matrix(method=method, matrix=matrix, config=config, run_id=run_id)
    except Exception as exc:  # defensive: worker must never crash the pool protocol
        return {
            "run_id": run_id,
            "matrix_id": matrix.matrix_id if type(matrix.matrix_id) is int else -1,
            "status": int(StatusCode.INTERNAL_ERROR),
            "success": False,
            "ess_count": 0,
            "elapsed_ns": 0,
            "candidate_count": 0,
            "error_message": f"Worker exception: {exc}",
            "candidates": [],
            "metadata": matrix.metadata,
        }


def _queue_worker(
    input_queue,
    output_queue,
    method: SearchMethod,
    config: RunConfig,
    run_id: str,
):
    """Consume serialized matrices and publish serialized result dictionaries."""

    while True:
        payload = input_queue.get()
        if payload is _SENTINEL:
            return

        matrix = pickle.loads(payload)
        result = _safe_compute(method=method, matrix=matrix, config=config, run_id=run_id)
        output_queue.put(bytes(ForkingPickler.dumps(result)))


def _max_pending_matrices(mp_config: MPConfig) -> int:
    """Return the bound for submitted-but-not-yielded matrices.

    The batch API may process millions of matrices. Submitting all matrices before
    draining results can deadlock bounded queues, so keep only a small window
    ahead of the consumer.
    """
    worker_window = mp_config.workers * mp_config.prefetch_per_worker
    return max(1, min(mp_config.queue_maxsize, worker_window))


class _QueueRunner:
    """Own the shared queues and worker processes for one wrapper run."""

    def __init__(self, method: SearchMethod, run_config: RunConfig, mp_config: MPConfig, run_id: str | None = None):
        """Start configured workers and create their shared queues."""

        self.method = method
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
                    args=(self._input_queue, self._output_queue, self.method, self.run_config, self.run_id),
                    daemon=True,
                )
                proc.start()
                self._workers.append(proc)
        except BaseException:
            self.shutdown(cancel=True)
            raise

    def submit(self, matrix: Matrix) -> None:
        """Serialize and enqueue one matrix.

        Raises:
            RuntimeError: If the input queue has already been closed.
            Exception: If the matrix or its metadata cannot be serialized.
        """

        if self._input_closed:
            raise RuntimeError("Cannot submit after close_input()")

        payload = bytes(ForkingPickler.dumps(matrix))
        self._input_queue.put(payload)

    def close_input(self) -> None:
        """Send one stop sentinel per worker; repeated calls do nothing."""

        if self._input_closed:
            return

        for _ in self._workers:
            self._input_queue.put(_SENTINEL)
        self._input_closed = True

    def get_result(self):
        """Wait for and deserialize the next completed result.

        Raises:
            RuntimeError: If workers exit before producing the expected result.
        """

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
        """Stop and join workers, then close both queues.

        Args:
            cancel: Terminate live workers instead of finishing queued input.
            join_timeout_s: Total graceful join budget before termination.
        """

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


def _run_matrices_multiprocessing(
    method: SearchMethod,
    matrices: Iterable[Matrix],
    config: RunConfig,
    mp_config: MPConfig,
    run_id: str,
) -> Iterator[dict]:
    """Yield multiprocessing results in completion order.

    Submission remains bounded by :func:`_max_pending_matrices`. Closing the
    iterator before completion terminates its workers.

    Raises:
        ValueError: If native logging is enabled for multiprocessing.
    """

    if config.enable_logging:
        raise ValueError("RunConfig.enable_logging is not supported with multiprocessing")

    runner = _QueueRunner(method=method, run_config=config, mp_config=mp_config, run_id=run_id)
    completed = False
    try:
        matrices_iter = iter(matrices)
        max_pending = _max_pending_matrices(mp_config)

        submitted = 0
        yielded = 0
        input_exhausted = False

        def submit_until_window() -> None:
            """Fill the bounded pending-work window or close exhausted input."""

            nonlocal submitted, input_exhausted
            while not input_exhausted and submitted - yielded < max_pending:
                try:
                    matrix = next(matrices_iter)
                except StopIteration:
                    input_exhausted = True
                    runner.close_input()
                    return

                runner.submit(matrix)
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


def run_multiprocessing(
    method: SearchMethod,
    matrices: Matrix | Iterable[Matrix],
    config: RunConfig | None = None,
    run_id: str | None = None,
    sink=None,
    mp_config: MPConfig | None = None,
) -> dict | Iterator[dict] | int:
    """Run matrix analysis across worker processes.

    A single :class:`Matrix` blocks and returns one result dictionary. An
    iterable returns a lazy completion-order iterator unless ``sink`` is
    provided; with a sink, all results are written eagerly and the number
    written is returned.

    Args:
        method: Required candidate-search method, ``"fast"`` or ``"safe"``.
        matrices: One matrix or an iterable of matrices.
        config: Analysis options; defaults to :class:`RunConfig`.
        run_id: Output identifier; a timestamp-based ID is generated when omitted.
        sink: Optional object providing ``write_result()``, ``close()``, and
            optionally ``abort()``.
        mp_config: Process scheduling options; defaults to :class:`MPConfig`.

    Returns:
        One result dictionary, a lazy result iterator, or a written-result count.

    Raises:
        ValueError: If ``config.enable_logging`` is true.
    """

    _validate_search_method(method)
    cfg = config if config is not None else RunConfig()
    mp_cfg = mp_config if mp_config is not None else MPConfig()
    rid = run_id or new_run_id("mp")
    is_single = isinstance(matrices, Matrix)
    source = (matrices,) if is_single else matrices
    results = _run_matrices_multiprocessing(
        method=method,
        matrices=source,
        config=cfg,
        mp_config=mp_cfg,
        run_id=rid,
    )

    if sink is not None:
        return _consume_to_sink(results, sink)
    if is_single:
        return list(results)[0]
    return results
