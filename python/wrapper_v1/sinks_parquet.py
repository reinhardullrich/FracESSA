from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from .types import MatrixResult


class ParquetSink:
    """
    Writes Parquet files for summary and candidate rows.
    Requires `pyarrow`.
    """

    def __init__(self, summary_path: str | Path, candidates_path: str | Path):
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise RuntimeError("ParquetSink requires pyarrow") from exc

        self._pa = pa
        self._pq = pq

        self._summary_path = Path(summary_path)
        self._candidate_path = Path(candidates_path)

        self._summary_writer = None
        self._candidate_writer = None

    def write_result(self, result: MatrixResult) -> None:
        summary_dict = asdict(result.summary)
        summary_table = self._pa.Table.from_pylist([summary_dict])

        if self._summary_writer is None:
            self._summary_writer = self._pq.ParquetWriter(
                self._summary_path,
                summary_table.schema,
            )
        self._summary_writer.write_table(summary_table)

        if result.candidates:
            candidate_rows = [asdict(c) for c in result.candidates]
            candidate_table = self._pa.Table.from_pylist(candidate_rows)
            if self._candidate_writer is None:
                self._candidate_writer = self._pq.ParquetWriter(
                    self._candidate_path,
                    candidate_table.schema,
                )
            self._candidate_writer.write_table(candidate_table)

    def close(self) -> None:
        if self._summary_writer is not None:
            self._summary_writer.close()
        if self._candidate_writer is not None:
            self._candidate_writer.close()
