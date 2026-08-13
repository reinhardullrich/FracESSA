import json
import tempfile
import unittest
from pathlib import Path

from pyfracessa.io import load_matrices_from_json


class IoTests(unittest.TestCase):
    def test_load_matrices_from_json_sorts_and_prefixes(self):
        payload = {
            "matrices": [
                {"id": 5, "dimension": 2, "matrix": "0,1,0", "tag": "b"},
                {"id": 1, "matrix": "3#4,13/2,1/2,5,11/2,3", "tag": "a"},
            ]
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrices.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            matrices = load_matrices_from_json(path)

        self.assertEqual([matrix.matrix_id for matrix in matrices], [1, 5])
        self.assertEqual(matrices[1].matrix, "2#0,1,0")
        self.assertEqual(matrices[0].metadata.get("tag"), "a")

    def test_load_matrices_from_json_requires_dimension_for_values_only(self):
        payload = [{"id": 9, "matrix": "0,1,0"}]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrices.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_matrices_from_json(path)

    def test_load_matrices_from_json_preserves_matrix_market(self):
        matrix_market = "%%MatrixMarket matrix array integer symmetric\n2 2\n1\n0\n1\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrices.json"
            path.write_text(json.dumps([{"id": 3, "matrix": matrix_market}]), encoding="utf-8")
            matrices = load_matrices_from_json(path)

        self.assertEqual(matrices[0].matrix, matrix_market.strip())

    def test_load_matrices_from_json_rejects_malformed_row_schema(self):
        cases = (
            ({"matrixes": []}, "must contain 'matrices'"),
            ({"matrices": {}}, "must be a list"),
            ({"matrices": ["not an object"]}, "must be an object"),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrices.json"
            for payload, message in cases:
                with self.subTest(payload=payload):
                    path.write_text(json.dumps(payload), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        load_matrices_from_json(path)

    def test_load_matrices_from_json_rejects_lossy_field_coercions(self):
        cases = (
            ({"id": 1.9, "matrix": "1#0"}, "matrix_id must be an int"),
            ({"id": True, "matrix": "1#0"}, "matrix_id must be an int"),
            ({"id": 1, "matrix": 0, "dimension": 1}, "matrix must be a str"),
            ({"id": 1, "matrix": "0", "dimension": 1.9}, "dimension.*must be an int"),
            ({"id": 1, "matrix": "0", "dimension": True}, "dimension.*must be an int"),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "matrices.json"
            for row, message in cases:
                with self.subTest(row=row):
                    path.write_text(json.dumps([row]), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        load_matrices_from_json(path)
