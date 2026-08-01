import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fracessa import Matrix, run
from fracessa.sinks import _JsonArrayWriter, _RowBuffer, create_sink
from fracessa.sinks_csv import CsvSink
from fracessa.sinks_json import JsonSink


def _sample_result(matrix_id: int = 3) -> dict:
    return {
        "run_id": "run1",
        "matrix_id": matrix_id,
        "status": 0,
        "success": True,
        "ess_count": 1,
        "elapsed_ns": 12,
        "candidate_count": 1,
        "error_message": "",
        "candidates": [
            {
                "run_id": "run1",
                "matrix_id": matrix_id,
                "candidate_id": 1,
                "vector": "1/2,1/2",
                "support": 3,
                "support_size": 2,
                "extended_support": 3,
                "extended_support_size": 2,
                "multiplier": None,
                "is_ess": True,
                "stability": "pure",
                "payoff": "1",
                "payoff_dbl": 1.0,
            }
        ],
        "metadata": {"source": "unit"},
    }


def _available_file_kinds() -> dict[str, str]:
    kinds = {"csv": "csv", "json": "json"}
    try:
        __import__("pyarrow")
    except ImportError:
        return kinds
    kinds["parquet"] = "parquet"
    return kinds


class SinkTests(unittest.TestCase):
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
            metadata = json.loads(
                Path(f"{summary_path}.metadata.json").read_text(encoding="utf-8")
            )

        self.assertGreaterEqual(len(summary_lines), 2)
        self.assertGreaterEqual(len(candidate_lines), 2)
        self.assertIn("multiplier", candidate_lines[0])
        self.assertEqual(candidate_rows[0]["multiplier"], "")
        self.assertEqual(metadata[0]["metadata"], {"source": "unit"})

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
            from fracessa.sinks_parquet import ParquetSink
        except ImportError:
            self.skipTest("pyarrow not installed")

        ordinary = _sample_result()
        circular = _sample_result()
        circular["candidates"][0]["multiplier"] = 5

        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = Path(tmpdir) / "summary.parquet"
            candidates_path = Path(tmpdir) / "candidates.parquet"
            sink = ParquetSink(summary_path, candidates_path)
            sink.write_result(ordinary)
            sink.write_result(circular)
            sink.close()
            multipliers = pq.read_table(candidates_path)["multiplier"].to_pylist()

        self.assertEqual(multipliers, [None, 5])

    def test_json_sink_replaces_non_finite_floats_with_null(self):
        result = _sample_result()
        result["candidates"][0]["payoff_dbl"] = float("inf")
        result["metadata"] = {"values": [float("nan"), float("-inf")]}

        def reject_non_standard_number(value):
            raise AssertionError(f"non-standard JSON number: {value}")

        with tempfile.TemporaryDirectory() as tmpdir:
            sink = create_sink("json", tmpdir, "finite_json")
            sink.write_result(result)
            sink.close()
            candidates = json.loads(
                (Path(tmpdir) / "finite_json_candidates.json").read_text(encoding="utf-8"),
                parse_constant=reject_non_standard_number,
            )
            metadata = json.loads(
                (Path(tmpdir) / "finite_json_json_metadata.json").read_text(encoding="utf-8"),
                parse_constant=reject_non_standard_number,
            )

        self.assertIsNone(candidates[0]["payoff_dbl"])
        self.assertEqual(metadata[0]["metadata"]["values"], [None, None])

    def test_sink_initialization_failure_rolls_back_and_allows_retry(self):
        calls = 0

        def fail_second_writer(output_file):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise RuntimeError("forced initialization failure")
            return _JsonArrayWriter(output_file)

        with tempfile.TemporaryDirectory() as tmpdir:
            paths = [
                Path(tmpdir) / "summary.json",
                Path(tmpdir) / "candidates.json",
                Path(tmpdir) / "metadata.json",
            ]
            with mock.patch(
                "fracessa.sinks_json._JsonArrayWriter",
                side_effect=fail_second_writer,
            ):
                with self.assertRaisesRegex(RuntimeError, "forced initialization failure"):
                    JsonSink(*paths)

            self.assertTrue(all(not path.exists() for path in paths))

            retry = JsonSink(*paths)
            retry.close()
            self.assertTrue(all(path.exists() for path in paths))

    def test_file_sink_write_failure_rolls_back_and_allows_retry(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind, extension in _available_file_kinds().items():
                with self.subTest(kind=kind):
                    run_id = f"write_failure_{kind}"
                    paths = [
                        Path(tmpdir) / f"{run_id}_summary.{extension}",
                        Path(tmpdir) / f"{run_id}_candidates.{extension}",
                        Path(tmpdir) / f"{run_id}_{kind}_metadata.json",
                    ]
                    result = _sample_result()
                    result["metadata"] = {"not_json": {1}}

                    sink = create_sink(kind, tmpdir, run_id)
                    with self.assertRaises(TypeError):
                        sink.write_result(result)

                    self.assertTrue(all(not path.exists() for path in paths))
                    retry = create_sink(kind, tmpdir, run_id)
                    retry.close()
                    self.assertTrue(all(path.exists() for path in paths))

    def test_file_sink_pipeline_failure_rolls_back_and_allows_retry(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind, extension in _available_file_kinds().items():
                with self.subTest(kind=kind):
                    run_id = f"pipeline_failure_{kind}"
                    paths = [
                        Path(tmpdir) / f"{run_id}_summary.{extension}",
                        Path(tmpdir) / f"{run_id}_candidates.{extension}",
                        Path(tmpdir) / f"{run_id}_{kind}_metadata.json",
                    ]
                    sink = create_sink(kind, tmpdir, run_id)

                    with mock.patch(
                        "fracessa.single.compute_matrix",
                        side_effect=RuntimeError("forced computation failure"),
                    ):
                        with self.assertRaisesRegex(RuntimeError, "forced computation failure"):
                            run(
                                "safe",
                                [Matrix(1, "2#0,1,0")],
                                sink=sink,
                                run_id=run_id,
                            )

                    self.assertTrue(all(not path.exists() for path in paths))
                    retry = create_sink(kind, tmpdir, run_id)
                    retry.close()
                    self.assertTrue(all(path.exists() for path in paths))

    def test_file_sink_close_failure_rolls_back_and_allows_retry(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_id = "close_failure_json"
            paths = [
                Path(tmpdir) / f"{run_id}_summary.json",
                Path(tmpdir) / f"{run_id}_candidates.json",
                Path(tmpdir) / f"{run_id}_json_metadata.json",
            ]
            sink = create_sink("json", tmpdir, run_id)
            original_close = sink._summary_writer.close

            def fail_after_close():
                original_close()
                raise RuntimeError("forced close failure")

            sink._summary_writer.close = fail_after_close
            with self.assertRaisesRegex(RuntimeError, "forced close failure"):
                sink.close()

            self.assertTrue(all(not path.exists() for path in paths))
            retry = create_sink("json", tmpdir, run_id)
            retry.close()
            self.assertTrue(all(path.exists() for path in paths))

    def test_create_sink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_sink = create_sink("csv", tmpdir, "rid")
            json_sink = create_sink("json", tmpdir, "rid")
            for sink in (csv_sink, json_sink):
                sink.write_result(_sample_result())
                sink.close()

            self.assertTrue((Path(tmpdir) / "rid_summary.csv").exists())
            self.assertTrue((Path(tmpdir) / "rid_candidates.csv").exists())
            metadata = json.loads(
                (Path(tmpdir) / "rid_csv_metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata[0]["metadata"], {"source": "unit"})
            self.assertTrue((Path(tmpdir) / "rid_json_metadata.json").exists())

    def test_create_sink_invalid_kind_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind in ("stream", "arrow", "bad-kind"):
                with self.subTest(kind=kind):
                    with self.assertRaises(ValueError):
                        create_sink(kind, tmpdir, "rid")

    def test_file_sinks_refuse_to_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind, extension in _available_file_kinds().items():
                with self.subTest(kind=kind):
                    run_id = f"existing_{kind}"
                    summary_path = Path(tmpdir) / f"{run_id}_summary.{extension}"
                    candidates_path = Path(tmpdir) / f"{run_id}_candidates.{extension}"
                    metadata_path = Path(tmpdir) / f"{run_id}_{kind}_metadata.json"
                    summary_path.write_bytes(b"existing summary")
                    candidates_path.write_bytes(b"existing candidates")
                    metadata_path.write_bytes(b"existing metadata")

                    with self.assertRaises(FileExistsError):
                        create_sink(kind, tmpdir, run_id)

                    self.assertEqual(summary_path.read_bytes(), b"existing summary")
                    self.assertEqual(candidates_path.read_bytes(), b"existing candidates")
                    self.assertEqual(metadata_path.read_bytes(), b"existing metadata")

                    summary_path.unlink()
                    with self.assertRaises(FileExistsError):
                        create_sink(kind, tmpdir, run_id)

                    self.assertFalse(summary_path.exists())
                    self.assertEqual(candidates_path.read_bytes(), b"existing candidates")

                    candidates_path.unlink()
                    with self.assertRaises(FileExistsError):
                        create_sink(kind, tmpdir, run_id)

                    self.assertFalse(summary_path.exists())
                    self.assertFalse(candidates_path.exists())
                    self.assertEqual(metadata_path.read_bytes(), b"existing metadata")

    def test_file_sinks_write_metadata_sidecars(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind in _available_file_kinds():
                with self.subTest(kind=kind):
                    run_id = f"metadata_{kind}"
                    sink = create_sink(kind, tmpdir, run_id)
                    sink.write_result(_sample_result())
                    without_metadata = _sample_result(4)
                    without_metadata["metadata"] = None
                    sink.write_result(without_metadata)
                    sink.close()

                    rows = json.loads(
                        (Path(tmpdir) / f"{run_id}_{kind}_metadata.json").read_text(
                            encoding="utf-8"
                        )
                    )
                    self.assertEqual(
                        rows,
                        [{"run_id": "run1", "matrix_id": 3, "metadata": {"source": "unit"}}],
                    )

    def test_file_sinks_write_readable_empty_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for kind, extension in _available_file_kinds().items():
                with self.subTest(kind=kind):
                    run_id = f"empty_{kind}"
                    create_sink(kind, tmpdir, run_id).close()
                    summary_path = Path(tmpdir) / f"{run_id}_summary.{extension}"
                    candidates_path = Path(tmpdir) / f"{run_id}_candidates.{extension}"

                    if kind == "csv":
                        self.assertTrue(summary_path.read_text(encoding="utf-8").startswith("run_id,"))
                        self.assertTrue(
                            candidates_path.read_text(encoding="utf-8").startswith("run_id,")
                        )
                    elif kind == "json":
                        self.assertEqual(json.loads(summary_path.read_text(encoding="utf-8")), [])
                        self.assertEqual(json.loads(candidates_path.read_text(encoding="utf-8")), [])
                    elif kind == "parquet":
                        import pyarrow.parquet as pq

                        self.assertEqual(pq.read_table(summary_path).num_rows, 0)
                        self.assertEqual(pq.read_table(candidates_path).num_rows, 0)

                    metadata_path = Path(tmpdir) / f"{run_id}_{kind}_metadata.json"
                    self.assertEqual(json.loads(metadata_path.read_text(encoding="utf-8")), [])

    def test_parquet_batches_rows(self):
        try:
            import pyarrow.parquet as pq
        except ImportError:
            self.skipTest("pyarrow is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            sink = create_sink("parquet", tmpdir, "batched")
            for matrix_id in range(1100):
                sink.write_result(_sample_result(matrix_id))
            sink.close()

            summary_file = pq.ParquetFile(Path(tmpdir) / "batched_summary.parquet")
            candidate_file = pq.ParquetFile(Path(tmpdir) / "batched_candidates.parquet")
            self.assertEqual(summary_file.num_row_groups, 2)
            self.assertEqual(candidate_file.num_row_groups, 2)
            self.assertEqual(summary_file.metadata.num_rows, 1100)
            self.assertEqual(candidate_file.metadata.num_rows, 1100)

    def test_row_buffer_splits_one_large_candidate_list(self):
        batch_sizes = []
        buffer = _RowBuffer(lambda rows: batch_sizes.append(len(rows)))

        buffer.extend([{}] * 1100)
        buffer.flush()

        self.assertEqual(batch_sizes, [1024, 76])
