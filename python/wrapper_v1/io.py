from __future__ import annotations

from collections.abc import Iterable, Iterator
import json
from pathlib import Path

from .types import MatrixJob


def load_jobs_from_verification_json(path: str | Path) -> list[MatrixJob]:
    """
    Load matrices from `verification_matrices.json` schema.
    """
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)

    matrices = raw.get("matrices", [])
    jobs: list[MatrixJob] = []
    for m in matrices:
        matrix_id = int(m["id"])
        dimension = int(m["dimension"])
        values = str(m["matrix"])
        matrix_cli = f"{dimension}#{values}"

        metadata = {
            "dimension": dimension,
            "number_ess": int(m.get("number_ess", 0)),
            "is_cs": bool(m.get("is_cs", False)),
        }
        jobs.append(MatrixJob(matrix_id=matrix_id, matrix=matrix_cli, metadata=metadata))

    jobs.sort(key=lambda x: x.matrix_id)
    return jobs


def load_jobs_from_json(
    path: str | Path,
    matrices_key: str = "matrices",
    id_key: str = "id",
    matrix_key: str = "matrix",
    dimension_key: str = "dimension",
) -> list[MatrixJob]:
    """
    Generic JSON loader.

    Supported input shapes:
    - {"matrices": [ ... ]}
    - [ ... ]

    If a matrix string has no `dimension#` prefix, `dimension` must be present.
    """
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)

    if isinstance(raw, dict):
        rows = raw.get(matrices_key, [])
    elif isinstance(raw, list):
        rows = raw
    else:
        raise ValueError(f"Unsupported JSON schema in {p}")

    jobs: list[MatrixJob] = []
    for row in rows:
        matrix_id = int(row[id_key])
        matrix_text = str(row[matrix_key]).strip()

        if "#" not in matrix_text:
            if dimension_key not in row:
                raise ValueError(
                    f"Row {matrix_id} has no 'dimension#' prefix and no '{dimension_key}' field"
                )
            matrix_text = f"{int(row[dimension_key])}#{matrix_text}"

        metadata = {
            k: v
            for k, v in row.items()
            if k not in {id_key, matrix_key}
        }
        jobs.append(MatrixJob(matrix_id=matrix_id, matrix=matrix_text, metadata=metadata))

    jobs.sort(key=lambda x: x.matrix_id)
    return jobs


def iter_jobs(jobs: Iterable[MatrixJob]) -> Iterator[MatrixJob]:
    for job in jobs:
        yield job
