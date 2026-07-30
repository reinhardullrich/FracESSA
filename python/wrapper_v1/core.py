from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
import threading

from .types import CandidateRow, MatrixJob, MatrixResult, RunConfig, StatusCode, SummaryRow

_native_module = None
_native_lock = threading.Lock()


def new_run_id(prefix: str = "run") -> str:
    return f"{prefix}_{datetime.now(tz=timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"


def _repo_root() -> Path:
    # python/wrapper_v1/core.py -> repo root is two levels up.
    return Path(__file__).resolve().parents[2]


def _native_search_paths() -> list[Path]:
    root = _repo_root()
    build = root / "cpp" / "build"
    return [
        build,
        build / "Release",
        build / "RelWithDebInfo",
        build / "Debug",
    ]


def load_native_module():
    global _native_module
    if _native_module is not None:
        return _native_module

    with _native_lock:
        if _native_module is not None:
            return _native_module

        try:
            _native_module = import_module("fracessa_core")
            return _native_module
        except ModuleNotFoundError:
            pass

        import sys

        for path in _native_search_paths():
            if not path.exists():
                continue
            path_str = str(path)
            if path_str not in sys.path:
                sys.path.insert(0, path_str)
            try:
                _native_module = import_module("fracessa_core")
                return _native_module
            except ModuleNotFoundError:
                continue

        raise ModuleNotFoundError(
            "Could not import native module 'fracessa_core'. Build it first with ./build.sh"
        )


def native_status_map() -> dict[str, int]:
    native = load_native_module()
    return {
        "OK": int(native.STATUS_OK),
        "PARSE_ERROR": int(native.STATUS_PARSE_ERROR),
        "DIMENSION_OUT_OF_RANGE": int(native.STATUS_DIMENSION_OUT_OF_RANGE),
        "INVALID_VALUE_COUNT": int(native.STATUS_INVALID_VALUE_COUNT),
        "EXEC_ERROR": int(native.STATUS_EXEC_ERROR),
        "INTERNAL_ERROR": int(native.STATUS_INTERNAL_ERROR),
    }


def _matrix_cli_string(job: MatrixJob) -> str:
    text = job.matrix.strip()
    if "#" in text:
        return text

    if not job.metadata:
        raise ValueError(
            "MatrixJob.matrix must be in CLI format 'dimension#values' when metadata is absent"
        )

    dimension = job.metadata.get("dimension")
    if dimension is None:
        raise ValueError(
            "MatrixJob.metadata['dimension'] is required when MatrixJob.matrix has no 'dimension#' prefix"
        )

    return f"{int(dimension)}#{text}"


def _candidate_row(run_id: str, matrix_id: int, row: dict) -> CandidateRow:
    return CandidateRow(
        run_id=run_id,
        matrix_id=matrix_id,
        candidate_id=int(row.get("candidate_id", 0)),
        vector=str(row.get("vector", "")),
        support=int(row.get("support", 0)),
        support_size=int(row.get("support_size", 0)),
        extended_support=int(row.get("extended_support", 0)),
        extended_support_size=int(row.get("extended_support_size", 0)),
        multiplier=(
            int(row["multiplier"]) if row.get("multiplier") is not None else None
        ),
        is_ess=bool(row.get("is_ess", False)),
        stability=str(row.get("stability", "")),
        payoff=str(row.get("payoff", "0")),
        payoff_dbl=float(row.get("payoff_dbl", 0.0)),
    )


def compute_job(job: MatrixJob, config: RunConfig, run_id: str) -> MatrixResult:
    native = load_native_module()
    matrix_cli = _matrix_cli_string(job)

    native_out = native.compute_matrix(
        matrix=matrix_cli,
        include_candidates=config.include_candidates,
        exact=config.exact,
        full_support=config.full_support,
        enable_logging=config.enable_logging,
        matrix_id=int(job.matrix_id),
    )

    status = int(native_out.get("status", StatusCode.INTERNAL_ERROR))
    success = bool(native_out.get("success", False))

    candidates: list[CandidateRow] = []
    for row in native_out.get("candidates", []):
        candidates.append(_candidate_row(run_id=run_id, matrix_id=job.matrix_id, row=row))

    elapsed_us = int(native_out.get("elapsed_us", 0)) if config.include_timing else 0

    summary = SummaryRow(
        run_id=run_id,
        matrix_id=int(job.matrix_id),
        status=status,
        success=success,
        ess_count=int(native_out.get("ess_count", 0)),
        elapsed_us=elapsed_us,
        candidate_count=len(candidates),
        error_message=str(native_out.get("error_message", "")),
    )

    return MatrixResult(summary=summary, candidates=candidates, metadata=job.metadata)


def result_to_dict(result: MatrixResult) -> dict:
    return {
        "summary": asdict(result.summary),
        "candidates": [asdict(c) for c in result.candidates],
        "metadata": result.metadata,
    }
