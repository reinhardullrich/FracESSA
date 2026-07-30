from __future__ import annotations

import csv
from pathlib import Path

from .sinks import (
    _CANDIDATE_FIELDS,
    _MetadataWriter,
    _SUMMARY_FIELDS,
    _abort_on_error,
    _close_all,
    _open_output_triplet,
)


class CsvSink:
    """
    Writes two CSV files and one metadata JSON file:
    - summary rows
    - candidate rows
    - per-matrix metadata
    """

    def __init__(
        self,
        summary_path: str | Path,
        candidates_path: str | Path,
        metadata_path: str | Path | None = None,
    ):
        with _open_output_triplet(
            summary_path,
            candidates_path,
            metadata_path,
            "x",
            encoding="utf-8",
            newline="",
        ) as output:
            files, self._rollback = output
            self._summary_file, self._candidates_file, metadata_file = files
            self._state = "active"
            self._metadata_writer = _MetadataWriter(metadata_file)

            self._summary_writer = csv.DictWriter(
                self._summary_file,
                fieldnames=_SUMMARY_FIELDS,
            )
            self._candidates_writer = csv.DictWriter(
                self._candidates_file,
                fieldnames=_CANDIDATE_FIELDS,
            )
            self._summary_writer.writeheader()
            self._candidates_writer.writeheader()

    def write_result(self, result: dict) -> None:
        if self._state != "active":
            raise RuntimeError("Cannot write to a closed or aborted CSV sink")
        with _abort_on_error(self):
            summary_dict = {
                key: value
                for key, value in result.items()
                if key not in {"candidates", "metadata"}
            }
            self._summary_writer.writerow(summary_dict)

            for candidate_dict in result["candidates"]:
                self._candidates_writer.writerow(candidate_dict)
            self._metadata_writer.write_result(result)

    def close(self) -> None:
        if self._state != "active":
            return
        with _abort_on_error(self):
            _close_all(self._summary_file, self._candidates_file, self._metadata_writer)
        self._state = "closed"
        self._rollback.pop_all()

    def abort(self) -> None:
        if self._state == "closed":
            return
        self._state = "aborted"
        self._rollback.close()
