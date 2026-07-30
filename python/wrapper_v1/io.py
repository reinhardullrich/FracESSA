from __future__ import annotations

import json
from pathlib import Path

from .types import MatrixJob


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
        if matrices_key not in raw:
            raise ValueError(f"JSON object in {p} must contain '{matrices_key}'")
        rows = raw[matrices_key]
    elif isinstance(raw, list):
        rows = raw
    else:
        raise ValueError(f"Unsupported JSON schema in {p}")

    if not isinstance(rows, list):
        raise ValueError(f"Matrix rows in {p} must be a list")

    jobs: list[MatrixJob] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"Matrix row {index} in {p} must be an object")

        metadata = {
            k: v
            for k, v in row.items()
            if k not in {id_key, matrix_key}
        }
        try:
            job = MatrixJob(
                matrix_id=row[id_key],
                matrix=row[matrix_key],
                metadata=metadata,
            )
        except KeyError as exc:
            raise ValueError(
                f"Matrix row {index} in {p} is missing field {exc.args[0]!r}"
            ) from exc
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid matrix row {index} in {p}: {exc}") from exc

        matrix_text = job.matrix.strip()

        if "#" not in matrix_text:
            if dimension_key not in row:
                raise ValueError(
                    f"Row {job.matrix_id} has no 'dimension#' prefix and no '{dimension_key}' field"
                )
            dimension = row[dimension_key]
            if type(dimension) is not int:
                raise ValueError(
                    f"Row {job.matrix_id} field '{dimension_key}' must be an int"
                )
            matrix_text = f"{dimension}#{matrix_text}"

        job.matrix = matrix_text
        jobs.append(job)

    jobs.sort(key=lambda x: x.matrix_id)
    return jobs
