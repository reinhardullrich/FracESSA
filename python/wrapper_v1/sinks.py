from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Protocol

from .types import MatrixResult


class ResultSink(Protocol):
    def write_result(self, result: MatrixResult) -> None:
        ...

    def close(self) -> None:
        ...


class StreamSink:
    """
    JSON-lines stream sink for live consumption.
    """

    def __init__(self, stream=None, flush_each: bool = False):
        self._stream = stream if stream is not None else sys.stdout
        self._flush_each = flush_each

    def write_result(self, result: MatrixResult) -> None:
        payload = {
            "summary": asdict(result.summary),
            "candidates": [asdict(c) for c in result.candidates],
            "metadata": result.metadata,
        }
        self._stream.write(json.dumps(payload, separators=(",", ":")) + "\n")
        if self._flush_each:
            self._stream.flush()

    def close(self) -> None:
        try:
            self._stream.flush()
        except Exception:
            pass


class MultiSink:
    def __init__(self, sinks: list[ResultSink]):
        self._sinks = sinks

    def write_result(self, result: MatrixResult) -> None:
        for sink in self._sinks:
            sink.write_result(result)

    def close(self) -> None:
        for sink in self._sinks:
            sink.close()


def create_sink(kind: str, output_dir: str | Path, run_id: str, stream=None):
    """
    Build a sink backend by name.

    Supported kinds:
    - csv
    - json
    - arrow
    - parquet
    - stream
    """
    sink_kind = kind.lower()
    if sink_kind == "stream":
        return StreamSink(stream=stream)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if sink_kind == "csv":
        from .sinks_csv import CsvSink

        return CsvSink(
            summary_path=out_dir / f"{run_id}_summary.csv",
            candidates_path=out_dir / f"{run_id}_candidates.csv",
        )

    if sink_kind == "json":
        from .sinks_json import JsonSink

        return JsonSink(
            summary_path=out_dir / f"{run_id}_summary.json",
            candidates_path=out_dir / f"{run_id}_candidates.json",
        )

    if sink_kind == "arrow":
        from .sinks_arrow import ArrowSink

        return ArrowSink(
            summary_path=out_dir / f"{run_id}_summary.arrow",
            candidates_path=out_dir / f"{run_id}_candidates.arrow",
        )

    if sink_kind == "parquet":
        from .sinks_parquet import ParquetSink

        return ParquetSink(
            summary_path=out_dir / f"{run_id}_summary.parquet",
            candidates_path=out_dir / f"{run_id}_candidates.parquet",
        )

    raise ValueError(f"Unsupported sink kind: {kind}")
