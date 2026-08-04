"""Provide shared PyFracESSA sink schemas and lifecycle helpers."""

from __future__ import annotations

from contextlib import ExitStack, contextmanager
import json
from itertools import islice
import math
from pathlib import Path


_SUMMARY_FIELDS = (
    "run_id",
    "matrix_id",
    "status",
    "candidate_count",
    "ess_count",
    "candidate_structure",
    "ess_structure",
    "elapsed_ns",
    "safe_fallback",
    "error_message",
)

_CANDIDATE_FIELDS = (
    "run_id",
    "matrix_id",
    "candidate_id",
    "vector",
    "support",
    "support_size",
    "extended_support",
    "extended_support_size",
    "multiplier",
    "is_ess",
    "stability",
    "payoff",
    "payoff_dbl",
)

_ROW_BATCH_SIZE = 1024


def _summary_row(result: dict, encode_structures: bool = False) -> dict:
    """Return the common summary fields, optionally flattening structures."""

    row = {name: result[name] for name in _SUMMARY_FIELDS}
    if encode_structures:
        for name in ("candidate_structure", "ess_structure"):
            row[name] = json.dumps(row[name], separators=(",", ":"), sort_keys=True)
    return row


def _json_value(value):
    """Return a JSON-safe copy with every non-finite float replaced by ``None``."""

    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _call_all(*actions) -> None:
    """Call every action and re-raise the first failure after all have run."""

    first_error = None
    for action in actions:
        try:
            action()
        except BaseException as exc:
            if first_error is None:
                first_error = exc
    if first_error is not None:
        raise first_error


def _close_all(*resources) -> None:
    """Close every resource and re-raise the first close failure."""

    _call_all(*(resource.close for resource in resources))


def _rollback_call(action) -> None:
    """Run one best-effort rollback action and suppress any failure."""

    try:
        action()
    except BaseException:
        pass


@contextmanager
def _abort_on_error(sink):
    """Abort ``sink`` if an exception escapes the managed operation."""

    try:
        yield
    except BaseException:
        _rollback_call(sink.abort)
        raise


def _consume_to_sink(results, sink) -> int:
    """Write an entire result iterable and finalize its sink.

    On failure, a closable result iterator is closed before the sink is aborted.
    If the sink has no ``abort`` method, ``close`` is used as cleanup instead.

    Args:
        results: Iterable of canonical result dictionaries.
        sink: Destination implementing ``write_result`` and ``close``.

    Returns:
        Number of results successfully written before finalization.
    """

    results = iter(results)
    count = 0
    try:
        for result in results:
            sink.write_result(result)
            count += 1
        sink.close()
    except BaseException:
        close_results = getattr(results, "close", None)
        if close_results is not None:
            _rollback_call(close_results)
        cleanup = getattr(sink, "abort", None) or sink.close
        _rollback_call(cleanup)
        raise
    return count


@contextmanager
def _open_output_triplet(
    summary_path,
    candidates_path,
    metadata_path,
    mode: str,
    **kwargs,
):
    """Open exclusive summary, candidate, and metadata outputs transactionally.

    The context yields ``((summary, candidates, metadata), rollback_stack)``.
    Unless the caller detaches the rollback stack after successful finalization,
    all opened files are closed and removed.
    """

    summary_path = Path(summary_path)
    candidates_path = Path(candidates_path)
    metadata_path = (
        Path(metadata_path)
        if metadata_path is not None
        else Path(f"{summary_path}.metadata.json")
    )
    rollback = ExitStack()
    try:
        opened = []
        for path, open_mode, open_kwargs in (
            (summary_path, mode, kwargs),
            (candidates_path, mode, kwargs),
            (metadata_path, "x", {"encoding": "utf-8"}),
        ):
            output_file = path.open(open_mode, **open_kwargs)
            rollback.callback(_rollback_call, path.unlink)
            rollback.callback(_rollback_call, output_file.close)
            opened.append(output_file)

        yield tuple(opened), rollback
    except BaseException:
        rollback.close()
        raise


class _JsonArrayWriter:
    """Stream dictionaries into one valid JSON array without retaining them."""

    def __init__(self, output_file):
        """Start a JSON array in an already-open text file."""

        self._file = output_file
        self._written = 0
        self._file.write("[\n")

    def write(self, row: dict) -> None:
        """Append one JSON-safe dictionary to the array."""

        payload = json.dumps(_json_value(row), ensure_ascii=True, allow_nan=False)
        if self._written:
            self._file.write(",\n")
        self._file.write(payload)
        self._written += 1

    def close(self) -> None:
        """Finish the JSON array and close its file."""

        _call_all(lambda: self._file.write("\n]\n"), self._file.close)


class _MetadataWriter:
    """Write non-null per-matrix metadata rows to a JSON array."""

    def __init__(self, output_file):
        """Wrap an already-open metadata text file."""

        self._writer = _JsonArrayWriter(output_file)

    def write_result(self, result: dict) -> None:
        """Write metadata from one result when metadata is present."""

        if result["metadata"] is not None:
            self._writer.write(
                {
                    "run_id": result["run_id"],
                    "matrix_id": result["matrix_id"],
                    "metadata": result["metadata"],
                }
            )

    def close(self) -> None:
        """Finish and close the metadata array."""

        self._writer.close()


class _RowBuffer:
    """Flush rows to a callback in fixed-size batches."""

    def __init__(self, write_rows):
        """Create a buffer whose batches are passed to ``write_rows``."""

        self._write_rows = write_rows
        self._rows = []

    def append(self, row: dict) -> None:
        """Append one row and flush when the batch is full."""

        self._rows.append(row)
        if len(self._rows) >= _ROW_BATCH_SIZE:
            self.flush()

    def extend(self, rows) -> None:
        """Append an iterable of rows while preserving the batch-size limit."""

        rows = iter(rows)
        while chunk := list(islice(rows, _ROW_BATCH_SIZE - len(self._rows))):
            self._rows.extend(chunk)
            if len(self._rows) == _ROW_BATCH_SIZE:
                self.flush()

    def flush(self) -> None:
        """Write and clear the current partial batch, if any."""

        if self._rows:
            self._write_rows(self._rows)
            self._rows.clear()


def create_sink(kind: str, output_dir: str | Path, run_id: str):
    """Create a CSV, JSON, or Parquet sink with run-specific file names.

    The output directory is created when needed. Each sink writes separate
    summary, candidate, and metadata files and refuses to overwrite any of them.

    Args:
        kind: Case-insensitive ``"csv"``, ``"json"``, or ``"parquet"``.
        output_dir: Directory that will contain the output triplet.
        run_id: Identifier used in all three file names.

    Returns:
        A new :class:`CsvSink`, :class:`JsonSink`, or :class:`ParquetSink`.

    Raises:
        ValueError: If ``kind`` is unsupported.
        FileExistsError: If any target output already exists.
    """

    sink_kind = kind.lower()
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if sink_kind == "csv":
        from .sinks_csv import CsvSink

        return CsvSink(
            summary_path=out_dir / f"{run_id}_summary.csv",
            candidates_path=out_dir / f"{run_id}_candidates.csv",
            metadata_path=out_dir / f"{run_id}_csv_metadata.json",
        )

    if sink_kind == "json":
        from .sinks_json import JsonSink

        return JsonSink(
            summary_path=out_dir / f"{run_id}_summary.json",
            candidates_path=out_dir / f"{run_id}_candidates.json",
            metadata_path=out_dir / f"{run_id}_json_metadata.json",
        )

    if sink_kind == "parquet":
        from .sinks_parquet import ParquetSink

        return ParquetSink(
            summary_path=out_dir / f"{run_id}_summary.parquet",
            candidates_path=out_dir / f"{run_id}_candidates.parquet",
            metadata_path=out_dir / f"{run_id}_parquet_metadata.json",
        )

    raise ValueError(f"Unsupported sink kind: {kind}")
