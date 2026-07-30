from __future__ import annotations

from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
import threading

from .types import Matrix, RunConfig

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


def _matrix_cli_string(matrix: Matrix) -> str:
    text = matrix.matrix.strip()
    if "#" in text:
        return text

    if not matrix.metadata:
        raise ValueError(
            "Matrix.matrix must be in CLI format 'dimension#values' when metadata is absent"
        )

    dimension = matrix.metadata.get("dimension")
    if dimension is None:
        raise ValueError(
            "Matrix.metadata['dimension'] is required when Matrix.matrix has no 'dimension#' prefix"
        )
    if type(dimension) is not int:
        raise TypeError("Matrix.metadata['dimension'] must be an int")

    return f"{dimension}#{text}"


def compute_matrix(matrix: Matrix, config: RunConfig, run_id: str) -> dict:
    native = load_native_module()
    matrix_cli = _matrix_cli_string(matrix)
    matrix_id = matrix.matrix_id

    native_out = native.compute_matrix(
        matrix=matrix_cli,
        include_candidates=config.include_candidates,
        exact=config.exact,
        full_support=config.full_support,
        enable_logging=config.enable_logging,
        matrix_id=matrix_id,
    )

    candidates = [
        {"run_id": run_id, "matrix_id": matrix_id, **candidate}
        for candidate in native_out["candidates"]
    ]
    return {
        "run_id": run_id,
        "matrix_id": matrix_id,
        "status": native_out["status"],
        "success": native_out["success"],
        "ess_count": native_out["ess_count"],
        "elapsed_ns": native_out["elapsed_ns"],
        "candidate_count": len(candidates),
        "error_message": native_out["error_message"],
        "candidates": candidates,
        "metadata": matrix.metadata,
    }
