from __future__ import annotations

from pathlib import Path

from .sinks import (
    _JsonArrayWriter,
    _MetadataWriter,
    _abort_on_error,
    _close_all,
    _open_output_triplet,
)


class JsonSink:
    """
    Writes summary, candidate, and metadata JSON arrays without retaining rows.
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
        ) as output:
            files, self._rollback = output
            summary_file, candidates_file, metadata_file = files
            self._state = "active"
            self._summary_writer = _JsonArrayWriter(summary_file)
            self._candidates_writer = _JsonArrayWriter(candidates_file)
            self._metadata_writer = _MetadataWriter(metadata_file)

    def write_result(self, result: dict) -> None:
        if self._state != "active":
            raise RuntimeError("Cannot write to a closed or aborted JSON sink")
        with _abort_on_error(self):
            summary_dict = {
                key: value
                for key, value in result.items()
                if key not in {"candidates", "metadata"}
            }
            self._summary_writer.write(summary_dict)

            for candidate_dict in result["candidates"]:
                self._candidates_writer.write(candidate_dict)
            self._metadata_writer.write_result(result)

    def close(self) -> None:
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
        if self._state == "closed":
            return
        self._state = "aborted"
        self._rollback.close()
