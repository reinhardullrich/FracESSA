"""Stream PyFracESSA results, candidates, and metadata to CSV."""

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
    """Write summary and candidate CSV files plus a metadata JSON sidecar.

    Files are opened exclusively. A successful :meth:`close` preserves the
    output triplet; :meth:`abort` closes and removes all three files.
    """

    def __init__(
        self,
        summary_path: str | Path,
        candidates_path: str | Path,
        metadata_path: str | Path | None = None,
    ):
        """Open a new transactional CSV output triplet.

        Args:
            summary_path: Destination for one summary row per matrix.
            candidates_path: Destination for zero or more candidate rows per matrix.
            metadata_path: Optional JSON sidecar path. By default it is derived
                from ``summary_path``.

        Raises:
            FileExistsError: If any destination already exists.
        """

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
        """Write one canonical result and all of its candidate and metadata rows."""

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
        """Finalize and preserve every output; repeated calls do nothing."""

        if self._state != "active":
            return
        with _abort_on_error(self):
            _close_all(self._summary_file, self._candidates_file, self._metadata_writer)
        self._state = "closed"
        self._rollback.pop_all()

    def abort(self) -> None:
        """Close and remove incomplete outputs; repeated calls are safe."""

        if self._state == "closed":
            return
        self._state = "aborted"
        self._rollback.close()
