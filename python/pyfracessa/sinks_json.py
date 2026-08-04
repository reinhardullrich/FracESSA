"""Stream PyFracESSA results, candidates, and metadata to JSON."""

from __future__ import annotations

from pathlib import Path

from .sinks import (
    _JsonArrayWriter,
    _MetadataWriter,
    _abort_on_error,
    _close_all,
    _open_output_triplet,
    _summary_row,
)


class JsonSink:
    """Write summary, candidate, and metadata arrays without retaining rows.

    Files are opened exclusively. A successful :meth:`close` preserves the
    output triplet; :meth:`abort` closes and removes all three files.
    """

    def __init__(
        self,
        summary_path: str | Path,
        candidates_path: str | Path,
        metadata_path: str | Path | None = None,
    ):
        """Open a new transactional JSON output triplet.

        Args:
            summary_path: Destination for the summary array.
            candidates_path: Destination for the candidate array.
            metadata_path: Optional metadata-array path. By default it is
                derived from ``summary_path``.

        Raises:
            FileExistsError: If any destination already exists.
        """

        with _open_output_triplet(
            summary_path,
            candidates_path,
            metadata_path,
            "x",
            encoding="utf-8",
        ) as output:
            files, self._rollback = output
            summary_file, candidates_file, metadata_file = files
            self._state = "active"
            self._summary_writer = _JsonArrayWriter(summary_file)
            self._candidates_writer = _JsonArrayWriter(candidates_file)
            self._metadata_writer = _MetadataWriter(metadata_file)

    def write_result(self, result: dict) -> None:
        """Write one canonical result and all of its candidate and metadata rows."""

        if self._state != "active":
            raise RuntimeError("Cannot write to a closed or aborted JSON sink")
        with _abort_on_error(self):
            self._summary_writer.write(_summary_row(result))

            for candidate_dict in result["candidates"]:
                self._candidates_writer.write(candidate_dict)
            self._metadata_writer.write_result(result)

    def close(self) -> None:
        """Finalize and preserve every JSON array; repeated calls do nothing."""

        if self._state != "active":
            return
        with _abort_on_error(self):
            _close_all(
                self._summary_writer,
                self._candidates_writer,
                self._metadata_writer,
            )
        self._state = "closed"
        self._rollback.pop_all()

    def abort(self) -> None:
        """Close and remove incomplete outputs; repeated calls are safe."""

        if self._state == "closed":
            return
        self._state = "aborted"
        self._rollback.close()
