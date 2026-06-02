from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any, Literal


class StatusCode(IntEnum):
    OK = 0
    PARSE_ERROR = 1
    DIMENSION_OUT_OF_RANGE = 2
    INVALID_VALUE_COUNT = 3
    EXEC_ERROR = 4
    INTERNAL_ERROR = 255


@dataclass(slots=True)
class MatrixJob:
    matrix_id: int
    matrix: str
    metadata: dict[str, Any] | None = None


@dataclass(slots=True)
class RunConfig:
    exact: bool = False
    full_support: bool = False
    include_candidates: bool = False
    include_timing: bool = True
    unsafe: bool = False
    enable_logging: bool = False


@dataclass(slots=True)
class MPConfig:
    workers: int
    chunk_size: int = 128
    queue_maxsize: int = 4096
    ordered_results: bool = True
    start_method: Literal["spawn", "forkserver", "fork"] = "spawn"

    def __post_init__(self) -> None:
        if self.workers < 1:
            raise ValueError("MPConfig.workers must be >= 1")
        if self.chunk_size < 1:
            raise ValueError("MPConfig.chunk_size must be >= 1")
        if self.queue_maxsize < 1:
            raise ValueError("MPConfig.queue_maxsize must be >= 1")


@dataclass(slots=True)
class SummaryRow:
    run_id: str
    matrix_id: int
    status: int
    success: bool
    ess_count: int
    elapsed_us: int
    candidate_count: int
    error_message: str


@dataclass(slots=True)
class CandidateRow:
    run_id: str
    matrix_id: int
    candidate_id: int
    vector: str
    support: int
    support_size: int
    extended_support: int
    extended_support_size: int
    shift_reference: int
    is_ess: bool
    stability: str
    payoff: str
    payoff_dbl: float


@dataclass(slots=True)
class MatrixResult:
    summary: SummaryRow
    candidates: list[CandidateRow] = field(default_factory=list)
    metadata: dict[str, Any] | None = None


@dataclass(slots=True)
class RunStats:
    run_id: str
    submitted: int
    completed: int
    ok: int
    failed: int
    elapsed_s: float
