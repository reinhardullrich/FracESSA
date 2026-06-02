from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .types import MatrixResult


class JsonSink:
    """
    Writes two JSON array files without keeping all rows in memory.
    """

    def __init__(self, summary_path: str | Path, candidates_path: str | Path):
        self._summary_file = Path(summary_path).open("w", encoding="utf-8")
        self._candidates_file = Path(candidates_path).open("w", encoding="utf-8")

        self._summary_file.write("[\n")
        self._candidates_file.write("[\n")

        self._summary_written = 0
        self._candidate_written = 0

    def write_result(self, result: MatrixResult) -> None:
        summary_dict = asdict(result.summary)
        if self._summary_written > 0:
            self._summary_file.write(",\n")
        self._summary_file.write(json.dumps(summary_dict, ensure_ascii=True))
        self._summary_written += 1

        for candidate in result.candidates:
            candidate_dict = asdict(candidate)
            if self._candidate_written > 0:
                self._candidates_file.write(",\n")
            self._candidates_file.write(json.dumps(candidate_dict, ensure_ascii=True))
            self._candidate_written += 1

    def close(self) -> None:
        self._summary_file.write("\n]\n")
        self._candidates_file.write("\n]\n")
        self._summary_file.close()
        self._candidates_file.close()
