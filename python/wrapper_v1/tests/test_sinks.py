import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from wrapper_v1.sinks import MultiSink, StreamSink, create_sink
from wrapper_v1.sinks_csv import CsvSink
from wrapper_v1.sinks_json import JsonSink
from wrapper_v1.types import CandidateRow, MatrixResult, SummaryRow


def _sample_result() -> MatrixResult:
    summary = SummaryRow(
        run_id="run1",
        matrix_id=3,
        status=0,
        success=True,
        ess_count=1,
        elapsed_us=12,
        candidate_count=1,
        error_message="",
    )
    candidate = CandidateRow(
        run_id="run1",
        matrix_id=3,
        candidate_id=1,
        vector="1/2,1/2",
        support=3,
        support_size=2,
        extended_support=3,
        extended_support_size=2,
        multiplier=None,
        is_ess=True,
        stability="pure",
        payoff="1",
        payoff_dbl=1.0,
    )
    return MatrixResult(summary=summary, candidates=[candidate], metadata={"source": "unit"})


class _CountingSink:
    def __init__(self):
        self.count = 0
        self.closed = False

    def write_result(self, result):
        self.count += 1

    def close(self):
        self.closed = True


class SinkTests(unittest.TestCase):
    def test_stream_sink_writes_json_lines(self):
        buffer = io.StringIO()
        sink = StreamSink(stream=buffer, flush_each=True)
        sink.write_result(_sample_result())
        sink.close()

        lines = buffer.getvalue().strip().splitlines()
        self.assertEqual(len(lines), 1)
        payload = json.loads(lines[0])
        self.assertEqual(payload["summary"]["matrix_id"], 3)

    def test_csv_sink_writes_summary_and_candidates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "summary.csv"
            candidates_path = Path(tmpdir) / "candidates.csv"
            sink = CsvSink(summary_path, candidates_path)
            sink.write_result(_sample_result())
            sink.close()

            summary_lines = summary_path.read_text(encoding="utf-8").strip().splitlines()
            candidate_lines = candidates_path.read_text(encoding="utf-8").strip().splitlines()
            with candidates_path.open(encoding="utf-8", newline="") as fh:
                candidate_rows = list(csv.DictReader(fh))

        self.assertGreaterEqual(len(summary_lines), 2)
        self.assertGreaterEqual(len(candidate_lines), 2)
        self.assertIn("multiplier", candidate_lines[0])
        self.assertEqual(candidate_rows[0]["multiplier"], "")

    def test_json_sink_writes_arrays(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "summary.json"
            candidates_path = Path(tmpdir) / "candidates.json"
            sink = JsonSink(summary_path, candidates_path)
            sink.write_result(_sample_result())
            sink.close()

            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            candidates = json.loads(candidates_path.read_text(encoding="utf-8"))

        self.assertEqual(len(summary), 1)
        self.assertEqual(len(candidates), 1)
        self.assertIsNone(candidates[0]["multiplier"])

    def test_parquet_multiplier_stays_nullable_across_batches(self):
        try:
            import pyarrow.parquet as pq
            from wrapper_v1.sinks_parquet import ParquetSink
        except ImportError:
            self.skipTest("pyarrow not installed")

        ordinary = _sample_result()
        circular = _sample_result()
        circular.candidates[0].multiplier = 5

        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "summary.parquet"
            candidates_path = Path(tmpdir) / "candidates.parquet"
            sink = ParquetSink(summary_path, candidates_path)
            sink.write_result(ordinary)
            sink.write_result(circular)
            sink.close()
            multipliers = pq.read_table(candidates_path)["multiplier"].to_pylist()

        self.assertEqual(multipliers, [None, 5])

    def test_create_sink_and_multi_sink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_sink = create_sink("csv", tmpdir, "rid")
            csv_sink.write_result(_sample_result())
            csv_sink.close()

            self.assertTrue((Path(tmpdir) / "rid_summary.csv").exists())
            self.assertTrue((Path(tmpdir) / "rid_candidates.csv").exists())

        sink_a = _CountingSink()
        sink_b = _CountingSink()
        multi = MultiSink([sink_a, sink_b])
        multi.write_result(_sample_result())
        multi.close()

        self.assertEqual(sink_a.count, 1)
        self.assertEqual(sink_b.count, 1)
        self.assertTrue(sink_a.closed)
        self.assertTrue(sink_b.closed)

    def test_create_sink_invalid_kind_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                create_sink("bad-kind", tmpdir, "rid")
