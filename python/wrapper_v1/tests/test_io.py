import json
import tempfile
import unittest
from pathlib import Path

from wrapper_v1.io import load_jobs_from_json, load_jobs_from_verification_json


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

    def test_load_jobs_from_verification_json(self):
        payload = {
            "matrices": [
                {"id": 2, "dimension": 2, "number_ess": 1, "is_cs": False, "matrix": "0,1,0"},
                {"id": 1, "dimension": 3, "number_ess": 2, "is_cs": True, "matrix": "1,3"},
            ]
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "verification_matrices.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            jobs = load_jobs_from_verification_json(path)

        self.assertEqual([j.matrix_id for j in jobs], [1, 2])
        self.assertEqual(jobs[0].matrix, "3#1,3")
        self.assertTrue(jobs[0].metadata["is_cs"])
