"""Define the public PyFracESSA data types and configuration."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import os
from typing import Any, Literal


class StatusCode(IntEnum):
    """Status values returned by native matrix computation."""

    OK = 0
    PARSE_ERROR = 1
    EXEC_ERROR = 4
    INTERNAL_ERROR = 255


@dataclass(slots=True)
class Matrix:
    """Describe one matrix submitted for analysis.

    Attributes:
        matrix_id: Signed 64-bit identifier preserved in every output row.
        matrix: Native matrix text, normally in ``dimension#values`` format.
        metadata: Optional user data copied to the result and metadata sidecar.
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
        mode: Candidate-search mode; verified is the default.
        full_support: Request full-support analysis from the native engine.
        include_candidates: Include candidate dictionaries in each result.
        enable_logging: Enable native logging for sequential execution only.
    """

    mode: Literal["verified", "exact", "unsafe", "very_unsafe"] = "verified"
    full_support: bool = False
    include_candidates: bool = False
    enable_logging: bool = False

    def __post_init__(self) -> None:
        """Reject unknown native analysis modes."""

        if not isinstance(self.mode, str):
            raise TypeError("RunConfig.mode must be a str")
        if self.mode not in {"verified", "exact", "unsafe", "very_unsafe"}:
            raise ValueError(
                "RunConfig.mode must be verified, exact, unsafe, or very_unsafe"
            )


@dataclass(slots=True)
class MPConfig:
    """Configure process scheduling for :func:`run_multiprocessing`.

    Attributes:
        workers: Number of worker processes; defaults to available CPUs.
        prefetch_per_worker: Maximum queued work window per worker.
        queue_maxsize: Maximum number of serialized results in the result queue.
        start_method: Python multiprocessing start method.
    """

    workers: int = max(1, os.process_cpu_count() or 1)
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
