"""Define the public PyFracESSA data types and configuration."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import os
from typing import Any, Literal

SearchMethod = Literal["fast", "safe", "test"]


def _validate_search_method(method: str) -> None:
    """Reject a missing or unknown candidate-search method."""

    if not isinstance(method, str):
        raise TypeError("method must be a str")
    if method not in {"fast", "safe", "test"}:
        raise ValueError("method must be fast, safe, or test")


class StatusCode(IntEnum):
    """Status values returned in every result dictionary.

    Attributes:
        OK: Analysis completed successfully.
        PARSE_ERROR: The native parser rejected the matrix text.
        EXEC_ERROR: Analysis failed after parsing.
        INTERNAL_ERROR: The Python worker or native boundary caught an unexpected failure.
    """

    OK = 0
    PARSE_ERROR = 1
    EXEC_ERROR = 4
    INTERNAL_ERROR = 255


@dataclass(slots=True)
class Matrix:
    """Describe one matrix submitted for analysis.

    Attributes:
        matrix_id: Signed 64-bit identifier preserved in every output row.
        matrix: Exact rational matrix text in ``dimension#values`` form. Values-only text is also accepted when ``metadata``
            contains an integer ``dimension``.
        metadata: Optional user dictionary copied unchanged to the result and, when used, the sink's metadata output.

    ``matrix`` accepts the full upper triangle of a symmetric matrix or the compact circular-symmetric layout described in
    :doc:`getting-started`. Entries must be integers or integer fractions; decimal notation is not accepted.
    """

    matrix_id: int
    matrix: str
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        """Validate the matrix input contract at construction time."""

        if type(self.matrix_id) is not int:
            raise TypeError("Matrix.matrix_id must be an int")
        if not -(1 << 63) <= self.matrix_id < (1 << 63):
            raise ValueError("Matrix.matrix_id must fit in a signed 64-bit integer")
        if not isinstance(self.matrix, str):
            raise TypeError("Matrix.matrix must be a str")
        if self.metadata is not None and not isinstance(self.metadata, dict):
            raise TypeError("Matrix.metadata must be a dict or None")


@dataclass(slots=True)
class RunConfig:
    """Configure native analysis for sequential or multiprocessing execution.

    Attributes:
        full_support: Check the full support first. If the selected candidate route finds it and exact stability accepts it as an
            ESS, stop; otherwise continue the ordinary search without checking that support twice.
        include_candidates: Include individual representative dictionaries in ``result["candidates"]``. Counts and support-size
            structures for the selected method are returned regardless of this setting.
        enable_logging: Write the native diagnostic trace for sequential execution. Multiprocessing rejects this option because
            workers cannot safely share the rotating log file.
    """

    full_support: bool = False
    include_candidates: bool = False
    enable_logging: bool = False

@dataclass(slots=True)
class MPConfig:
    """Configure process scheduling for :func:`run_multiprocessing`.

    Attributes:
        workers: Number of worker processes. The default is the CPUs available to the Python process; each worker analyzes one
            matrix at a time.
        prefetch_per_worker: Pending-matrix allowance per worker. Together with ``queue_maxsize``, this bounds submitted work.
        queue_maxsize: Maximum serialized-result queue size and upper bound on the pending-matrix window.
        start_method: Python multiprocessing start method. The portable default, ``"spawn"``, requires the calling script to use
            an ``if __name__ == "__main__":`` guard.
    """

    workers: int = max(1, getattr(os, "process_cpu_count", os.cpu_count)() or 1)
    prefetch_per_worker: int = 128
    queue_maxsize: int = 4096
    start_method: Literal["spawn", "forkserver", "fork"] = "spawn"

    def __post_init__(self) -> None:
        """Reject non-positive worker and queue settings."""

        if self.workers < 1:
            raise ValueError("MPConfig.workers must be >= 1")
        if self.prefetch_per_worker < 1:
            raise ValueError("MPConfig.prefetch_per_worker must be >= 1")
        if self.queue_maxsize < 1:
            raise ValueError("MPConfig.queue_maxsize must be >= 1")
