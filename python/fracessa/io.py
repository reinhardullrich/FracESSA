"""Load PyFracESSA matrix inputs from JSON files."""

from __future__ import annotations

import json
from pathlib import Path

from .types import Matrix


def load_matrices_from_json(
    path: str | Path,
    matrices_key: str = "matrices",
    id_key: str = "id",
    matrix_key: str = "matrix",
    dimension_key: str = "dimension",
) -> list[Matrix]:
    """Load, validate, normalize, and sort matrices from a JSON file.

    Supported input shapes:

    - ``{"matrices": [...]}``
    - ``[...]``

    Field names can be remapped with the key arguments. Extra row fields are
    preserved as metadata. Values-only matrix strings require an integer
    dimension field and are normalized to ``dimension#values``.

    Args:
        path: JSON file to read.
        matrices_key: Top-level key containing rows when the document is an object.
        id_key: Row key containing the signed 64-bit matrix identifier.
        matrix_key: Row key containing the matrix string.
        dimension_key: Row key used to prefix values-only matrix strings.

    Returns:
        Validated matrices sorted by ``matrix_id``.

    Raises:
        ValueError: If the JSON shape or any matrix row is invalid.
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

    matrices: list[Matrix] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"Matrix row {index} in {p} must be an object")

        metadata = {
            k: v
            for k, v in row.items()
            if k not in {id_key, matrix_key}
        }
        try:
            matrix = Matrix(
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

        matrix_text = matrix.matrix.strip()

        if "#" not in matrix_text:
            if dimension_key not in row:
                raise ValueError(
                    f"Row {matrix.matrix_id} has no 'dimension#' prefix and no '{dimension_key}' field"
                )
            dimension = row[dimension_key]
            if type(dimension) is not int:
                raise ValueError(
                    f"Row {matrix.matrix_id} field '{dimension_key}' must be an int"
                )
            matrix_text = f"{dimension}#{matrix_text}"

        matrix.matrix = matrix_text
        matrices.append(matrix)

    matrices.sort(key=lambda x: x.matrix_id)
    return matrices
