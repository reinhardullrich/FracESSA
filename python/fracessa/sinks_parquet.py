"""Write batched PyFracESSA Parquet output with JSON metadata."""

from __future__ import annotations

from pathlib import Path

from .sinks import (
    _MetadataWriter,
    _RowBuffer,
    _SUMMARY_FIELDS,
    _abort_on_error,
    _call_all,
    _open_output_triplet,
    _rollback_call,
)


class ParquetSink:
    """Write batched summary and candidate Parquet files plus metadata JSON.

    Rows are buffered in fixed-size batches before being passed to PyArrow.
    Files are opened exclusively and removed by :meth:`abort` after any failure.
    """

    def __init__(
        self,
        summary_path: str | Path,
        candidates_path: str | Path,
        metadata_path: str | Path | None = None,
    ):
        """Open a new transactional Parquet output triplet.

        Args:
            summary_path: Destination for the summary Parquet table.
            candidates_path: Destination for the candidate Parquet table.
            metadata_path: Optional JSON sidecar path. By default it is derived
                from ``summary_path``.

        Raises:
            RuntimeError: If PyArrow is unavailable.
            FileExistsError: If any destination already exists.
        """

        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise RuntimeError("ParquetSink requires pyarrow") from exc

        self._pa = pa

        with _open_output_triplet(
            summary_path,
            candidates_path,
            metadata_path,
            "xb",
        ) as output:
            files, self._rollback = output
            self._summary_file, self._candidate_file, metadata_file = files
            self._state = "active"
            self._metadata_writer = _MetadataWriter(metadata_file)

            self._summary_schema = pa.schema([
                ("run_id", pa.string()),
                ("matrix_id", pa.int64()),
                ("status", pa.int32()),
                ("ess_count", pa.int64()),
                ("elapsed_ns", pa.int64()),
                ("safe_fallback", pa.string()),
                ("candidate_count", pa.int64()),
                ("error_message", pa.string()),
            ])

            self._candidate_schema = pa.schema([
                ("run_id", pa.string()),
                ("matrix_id", pa.int64()),
                ("candidate_id", pa.int64()),
                ("vector", pa.string()),
                ("support", pa.uint64()),
                ("support_size", pa.int64()),
                ("extended_support", pa.uint64()),
                ("extended_support_size", pa.int64()),
                ("multiplier", pa.int64()),
                ("is_ess", pa.bool_()),
                ("stability", pa.string()),
                ("payoff", pa.string()),
                ("payoff_dbl", pa.float64()),
            ])

            self._summary_writer = pq.ParquetWriter(
                self._summary_file,
                self._summary_schema,
            )
            self._rollback.callback(_rollback_call, self._summary_writer.close)
            self._candidate_writer = pq.ParquetWriter(
                self._candidate_file,
                self._candidate_schema,
            )
            self._rollback.callback(_rollback_call, self._candidate_writer.close)
            self._summary_buffer = _RowBuffer(self._write_summary_rows)
            self._candidate_buffer = _RowBuffer(self._write_candidate_rows)

    def _write_summary_rows(self, rows) -> None:
        """Convert and append one batch of summary dictionaries."""

        table = self._pa.Table.from_pylist(rows, schema=self._summary_schema)
        self._summary_writer.write_table(table)

    def _write_candidate_rows(self, rows) -> None:
        """Convert and append one batch of candidate dictionaries."""

        table = self._pa.Table.from_pylist(rows, schema=self._candidate_schema)
        self._candidate_writer.write_table(table)

    def write_result(self, result: dict) -> None:
        """Buffer one canonical result and its candidate and metadata rows."""

        if self._state != "active":
            raise RuntimeError("Cannot write to a closed or aborted Parquet sink")
        with _abort_on_error(self):
            self._summary_buffer.append({name: result[name] for name in _SUMMARY_FIELDS})
            self._candidate_buffer.extend(result["candidates"])
            self._metadata_writer.write_result(result)

    def close(self) -> None:
        """Flush, finalize, and preserve outputs; repeated calls do nothing."""

        if self._state != "active":
            return
        with _abort_on_error(self):
            _call_all(
                self._summary_buffer.flush,
                self._candidate_buffer.flush,
                self._summary_writer.close,
                self._candidate_writer.close,
                self._summary_file.close,
                self._candidate_file.close,
                self._metadata_writer.close,
            )
        self._state = "closed"
        self._rollback.pop_all()

    def abort(self) -> None:
        """Close and remove incomplete outputs; repeated calls are safe."""

        if self._state == "closed":
            return
        self._state = "aborted"
        self._rollback.close()
