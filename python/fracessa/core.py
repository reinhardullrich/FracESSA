"""Adapt PyFracESSA objects to the native Pybind module."""

from __future__ import annotations

from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
import threading

from .types import Matrix, RunConfig, SearchMethod, _validate_search_method

_native_module = None
_native_lock = threading.Lock()


def new_run_id(prefix: str = "run") -> str:
    """Return a UTC timestamp-based identifier using ``prefix``."""

    return f"{prefix}_{datetime.now(tz=timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"


def _repo_root() -> Path:
    """Return the repository root containing the Python and C++ sources."""

    return Path(__file__).resolve().parents[2]


def _native_search_paths() -> list[Path]:
    """Return local CMake build directories that may contain ``fracessa_core``."""

    root = _repo_root()
    build = root / "cpp" / "build"
    return [
        build,
        build / "Release",
        build / "RelWithDebInfo",
        build / "Debug",
    ]


def load_native_module():
    """Load and cache the ``fracessa_core`` extension module.

    The normal Python import path is tried first, followed by the repository's
    common CMake build directories.

    Returns:
        The imported native extension module.

    Raises:
        ModuleNotFoundError: If ``fracessa_core`` cannot be found.
    """

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
    """Return ``matrix`` in the native ``dimension#values`` text format.

    A matrix that already contains ``#`` is returned after trimming surrounding
    whitespace. Otherwise its dimension is read from ``matrix.metadata``.

    Raises:
        ValueError: If a values-only matrix has no dimension metadata.
        TypeError: If the metadata dimension is not a built-in integer.
    """

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


def compute_matrix(method: SearchMethod, matrix: Matrix, config: RunConfig, run_id: str) -> dict:
    """Compute one matrix with the native extension and normalize its result.

    Candidate rows are augmented with ``run_id`` and ``matrix_id``. The returned
    dictionary is the canonical result shape consumed by every file sink.

    Args:
        method: Required candidate-search method: ``"fast"``, ``"safe"``, or experimental ``"test"``.
        matrix: Validated matrix input.
        config: Native analysis options.
        run_id: Identifier attached to the result and candidate rows.

    Returns:
        A flat result dictionary containing status, timing, counts, candidates,
        and the input metadata.
    """

    _validate_search_method(method)
    native = load_native_module()
    matrix_cli = _matrix_cli_string(matrix)
    matrix_id = matrix.matrix_id

    native_out = native.compute_matrix(
        method=method,
        matrix=matrix_cli,
        include_candidates=config.include_candidates,
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
        "safe_fallback": native_out["safe_fallback"],
        "candidate_count": len(candidates),
        "error_message": native_out["error_message"],
        "candidates": candidates,
        "metadata": matrix.metadata,
    }
