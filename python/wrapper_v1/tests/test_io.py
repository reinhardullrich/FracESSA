import json
import tempfile
import unittest
from pathlib import Path

from wrapper_v1.io import load_jobs_from_json


class IoTests(unittest.TestCase):
    def test_load_jobs_from_json_sorts_and_prefixes(self):
        payload = {
            "matrices": [
                {"id": 5, "dimension": 2, "matrix": "0,1,0", "tag": "b"},
                {"id": 1, "matrix": "3#4,13/2,1/2,5,11/2,3", "tag": "a"},
            ]
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "jobs.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            jobs = load_jobs_from_json(path)

        self.assertEqual([j.matrix_id for j in jobs], [1, 5])
        self.assertEqual(jobs[1].matrix, "2#0,1,0")
        self.assertEqual(jobs[0].metadata.get("tag"), "a")

    def test_load_jobs_from_json_requires_dimension_for_values_only(self):
        payload = [{"id": 9, "matrix": "0,1,0"}]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "jobs.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_jobs_from_json(path)

    def test_load_jobs_from_json_rejects_malformed_row_schema(self):
        cases = (
            ({"matrixes": []}, "must contain 'matrices'"),
            ({"matrices": {}}, "must be a list"),
            ({"matrices": ["not an object"]}, "must be an object"),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "jobs.json"
            for payload, message in cases:
                with self.subTest(payload=payload):
                    path.write_text(json.dumps(payload), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        load_jobs_from_json(path)

    def test_load_jobs_from_json_rejects_lossy_field_coercions(self):
        cases = (
            ({"id": 1.9, "matrix": "1#0"}, "matrix_id must be an int"),
            ({"id": True, "matrix": "1#0"}, "matrix_id must be an int"),
            ({"id": 1, "matrix": 0, "dimension": 1}, "matrix must be a str"),
            ({"id": 1, "matrix": "0", "dimension": 1.9}, "dimension.*must be an int"),
            ({"id": 1, "matrix": "0", "dimension": True}, "dimension.*must be an int"),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "jobs.json"
            for row, message in cases:
                with self.subTest(row=row):
                    path.write_text(json.dumps([row]), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        load_jobs_from_json(path)
