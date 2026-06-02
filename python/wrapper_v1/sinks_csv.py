from __future__ import annotations

from dataclasses import asdict
import csv
from pathlib import Path

from .types import MatrixResult


class CsvSink:
    """
    Writes two CSV files:
    - summary rows
    - candidate rows
    """

    def __init__(self, summary_path: str | Path, candidates_path: str | Path):
        self._summary_file = Path(summary_path).open("w", encoding="utf-8", newline="")
        self._candidates_file = Path(candidates_path).open("w", encoding="utf-8", newline="")

        self._summary_writer = None
        self._candidates_writer = None

    def write_result(self, result: MatrixResult) -> None:
        summary_dict = asdict(result.summary)
        if self._summary_writer is None:
            self._summary_writer = csv.DictWriter(
                self._summary_file,
                fieldnames=list(summary_dict.keys()),
            )
            self._summary_writer.writeheader()
        self._summary_writer.writerow(summary_dict)

        for candidate in result.candidates:
            candidate_dict = asdict(candidate)
            if self._candidates_writer is None:
                self._candidates_writer = csv.DictWriter(
                    self._candidates_file,
                    fieldnames=list(candidate_dict.keys()),
                )
                self._candidates_writer.writeheader()
            self._candidates_writer.writerow(candidate_dict)

    def close(self) -> None:
        self._summary_file.close()
        self._candidates_file.close()
