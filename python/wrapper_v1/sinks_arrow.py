from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from .types import MatrixResult


class ArrowSink:
    """
    Writes Arrow IPC files (stream format) for summary and candidate rows.
    Requires `pyarrow`.
    """

    def __init__(self, summary_path: str | Path, candidates_path: str | Path):
        try:
            import pyarrow as pa
            import pyarrow.ipc as ipc
        except ImportError as exc:
            raise RuntimeError("ArrowSink requires pyarrow") from exc

        self._pa = pa
        self._ipc = ipc

        self._summary_file = Path(summary_path).open("wb")
        self._candidates_file = Path(candidates_path).open("wb")

        self._summary_schema = pa.schema([
            ("run_id", pa.string()),
            ("matrix_id", pa.int64()),
            ("status", pa.int32()),
            ("success", pa.bool_()),
            ("ess_count", pa.int64()),
            ("elapsed_us", pa.int64()),
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
            ("shift_reference", pa.int64()),
            ("is_ess", pa.bool_()),
            ("stability", pa.string()),
            ("payoff", pa.string()),
            ("payoff_dbl", pa.float64()),
        ])

        self._summary_writer = ipc.new_stream(self._summary_file, self._summary_schema)
        self._candidate_writer = ipc.new_stream(self._candidates_file, self._candidate_schema)

    def write_result(self, result: MatrixResult) -> None:
        summary_dict = asdict(result.summary)
        summary_batch = self._pa.RecordBatch.from_pylist([summary_dict], schema=self._summary_schema)
        self._summary_writer.write_batch(summary_batch)

        if result.candidates:
            candidate_rows = [asdict(c) for c in result.candidates]
            candidate_batch = self._pa.RecordBatch.from_pylist(candidate_rows, schema=self._candidate_schema)
            self._candidate_writer.write_batch(candidate_batch)

    def close(self) -> None:
        self._summary_writer.close()
        self._candidate_writer.close()
        self._summary_file.close()
        self._candidates_file.close()
