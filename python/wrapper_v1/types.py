from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Any, Literal


class StatusCode(IntEnum):
    OK = 0
    PARSE_ERROR = 1
    EXEC_ERROR = 4
    INTERNAL_ERROR = 255


@dataclass(slots=True)
class MatrixJob:
    matrix_id: int
    matrix: str
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if type(self.matrix_id) is not int:
            raise TypeError("MatrixJob.matrix_id must be an int")
        if not -(1 << 63) <= self.matrix_id < (1 << 63):
            raise ValueError("MatrixJob.matrix_id must fit in a signed 64-bit integer")
        if not isinstance(self.matrix, str):
            raise TypeError("MatrixJob.matrix must be a str")
        if self.metadata is not None and not isinstance(self.metadata, dict):
            raise TypeError("MatrixJob.metadata must be a dict or None")


@dataclass(slots=True)
class RunConfig:
    exact: bool = False
    full_support: bool = False
    include_candidates: bool = False
    enable_logging: bool = False


@dataclass(slots=True)
class MPConfig:
    workers: int
    prefetch_per_worker: int = 128
    queue_maxsize: int = 4096
    start_method: Literal["spawn", "forkserver", "fork"] = "spawn"

    def __post_init__(self) -> None:
        if self.workers < 1:
            raise ValueError("MPConfig.workers must be >= 1")
        if self.prefetch_per_worker < 1:
            raise ValueError("MPConfig.prefetch_per_worker must be >= 1")
        if self.queue_maxsize < 1:
            raise ValueError("MPConfig.queue_maxsize must be >= 1")
