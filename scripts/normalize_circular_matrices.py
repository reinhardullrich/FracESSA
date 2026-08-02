#!/usr/bin/env python3
"""Store eligible circulant matrices compactly and remove newer normalized duplicates."""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import sqlite3


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = REPOSITORY_ROOT / "testdata" / "fracessa_testdata.sqlite3"
NORMALIZATION_PREFIX = "Circular-storage normalization:"
DIMENSION_ONE_NOTE = (
    "Dimension-one circular exception: this matrix is mathematically circulant, but FracESSA's compact format would contain "
    "floor(1/2) = 0 values while the parser requires a value, so it remains in full non-circular storage."
)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _expand_matrix(dimension: int, is_circular: bool, text: str) -> list[list[Fraction]]:
    values = [Fraction(token) for token in text.split(",")]
    if is_circular:
        if len(values) != dimension // 2:
            raise ValueError(f"compact dimension {dimension} requires {dimension // 2} values, got {len(values)}")
        first_row = [Fraction(0) for _ in range(dimension)]
        for distance, value in enumerate(values, start=1):
            first_row[distance] = value
            first_row[dimension - distance] = value
        return [[first_row[(column - row) % dimension] for column in range(dimension)] for row in range(dimension)]

    expected = dimension * (dimension + 1) // 2
    if len(values) != expected:
        raise ValueError(f"upper-triangle dimension {dimension} requires {expected} values, got {len(values)}")
    matrix = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    position = 0
    for row in range(dimension):
        for column in range(row, dimension):
            matrix[row][column] = values[position]
            matrix[column][row] = values[position]
            position += 1
    return matrix


def _is_circulant(matrix: list[list[Fraction]]) -> bool:
    dimension = len(matrix)
    return all(matrix[row][column] == matrix[0][(column - row) % dimension] for row in range(dimension) for column in range(dimension))


def _compact_normalized(matrix: list[list[Fraction]]) -> tuple[Fraction, str]:
    dimension = len(matrix)
    diagonal = matrix[0][0]
    normalized = [[value - diagonal for value in row] for row in matrix]
    compact_values = [normalized[0][distance] for distance in range(1, dimension // 2 + 1)]
    compact_text = ",".join(_fraction_text(value) for value in compact_values)
    if _expand_matrix(dimension, True, compact_text) != normalized:
        raise RuntimeError("compact circular reconstruction does not match the normalized matrix")
    return diagonal, compact_text


def _append_note(description: str | None, note: str) -> str:
    if description and note in description:
        return description
    return f"{description} {note}" if description else note


def _audit(connection: sqlite3.Connection):
    already_compact = []
    conversions = []
    dimension_one = []
    canonical_rows: dict[tuple[int, str], list[tuple[str | None, int]]] = {}
    for matrix_id, dimension, is_circular, text, description, created_at in connection.execute(
        "SELECT matrix_id, dimension, is_cs, matrix, description, created_at FROM matrices ORDER BY matrix_id"
    ):
        matrix = _expand_matrix(dimension, bool(is_circular), text)
        if not _is_circulant(matrix):
            if is_circular:
                raise RuntimeError(f"matrix {matrix_id} is marked circular but is not circulant")
            continue
        if is_circular:
            if any(matrix[index][index] != 0 for index in range(dimension)):
                raise RuntimeError(f"compact matrix {matrix_id} has a nonzero diagonal")
            already_compact.append(matrix_id)
            canonical_text = ",".join(_fraction_text(matrix[0][distance]) for distance in range(1, dimension // 2 + 1))
            canonical_rows.setdefault((dimension, canonical_text), []).append((created_at, matrix_id))
            continue
        if dimension == 1:
            dimension_one.append((matrix_id, _append_note(description, DIMENSION_ONE_NOTE)))
            continue

        diagonal, compact_text = _compact_normalized(matrix)
        diagonal_text = _fraction_text(diagonal)
        note = (
            f"{NORMALIZATION_PREFIX} subtracted the original diagonal value {diagonal_text} from every matrix entry, "
            "producing the strategically equivalent zero-diagonal game stored here in compact circular form. The preceding "
            f"provenance describes the unnormalized matrix; each stored payoff equals its original payoff minus {diagonal_text}."
        )
        conversions.append((matrix_id, dimension, diagonal_text, compact_text, _append_note(description, note)))
        canonical_rows.setdefault((dimension, compact_text), []).append((created_at, matrix_id))

    duplicate_rows = []
    for rows in canonical_rows.values():
        rows.sort(key=lambda row: (row[0] is None, row[0] or "", row[1]))
        duplicate_rows.extend((matrix_id, rows[0][1]) for _, matrix_id in rows[1:])
    return already_compact, conversions, dimension_one, duplicate_rows


def normalize(database: Path, apply_changes: bool) -> None:
    if not database.is_file():
        raise ValueError(f"database does not exist: {database}")
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        already_compact, conversions, dimension_one, duplicate_rows = _audit(connection)
        print(
            f"already_compact={len(already_compact)} conversions={len(conversions)} "
            f"dimension_one_exceptions={len(dimension_one)} duplicates={len(duplicate_rows)}"
        )
        for matrix_id, dimension, diagonal, compact_text, _ in conversions:
            print(f"matrix={matrix_id} dimension={dimension} diagonal={diagonal} compact={compact_text}")
        for matrix_id, _ in dimension_one:
            print(f"matrix={matrix_id} dimension=1 retained_non_circular")
        for matrix_id, retained_id in duplicate_rows:
            print(f"matrix={matrix_id} duplicate_of={retained_id} remove_newer=true")
        if not apply_changes:
            print("dry_run=true; pass --apply to update the database")
            return

        affected_ids = {matrix_id for matrix_id, *_ in conversions} | {matrix_id for matrix_id, _ in duplicate_rows}
        candidate_rows = sum(
            connection.execute("SELECT COUNT(*) FROM candidates WHERE matrix_id = ?", (matrix_id,)).fetchone()[0]
            for matrix_id in affected_ids
        )
        timing_rows = sum(
            connection.execute("SELECT COUNT(*) FROM timings WHERE matrix_id = ?", (matrix_id,)).fetchone()[0]
            for matrix_id in affected_ids
        )
        with connection:
            for matrix_id, description in dimension_one:
                connection.execute("UPDATE matrices SET description = ? WHERE matrix_id = ?", (description, matrix_id))
            for matrix_id, _, _, compact_text, description in conversions:
                connection.execute("DELETE FROM candidates WHERE matrix_id = ?", (matrix_id,))
                connection.execute("DELETE FROM timings WHERE matrix_id = ?", (matrix_id,))
                cursor = connection.execute(
                    """UPDATE matrices
                       SET is_cs = 1, matrix = ?, description = ?,
                           candidate_count = NULL, ess_count = NULL,
                           candidate_structure = NULL, ess_structure = NULL,
                           fast_calibration_ns = NULL, safe_calibration_ns = NULL
                       WHERE matrix_id = ? AND is_cs = 0""",
                    (compact_text, description, matrix_id),
                )
                if cursor.rowcount != 1:
                    raise RuntimeError(f"matrix {matrix_id} changed during circular normalization")
            for matrix_id, _ in duplicate_rows:
                cursor = connection.execute("DELETE FROM matrices WHERE matrix_id = ?", (matrix_id,))
                if cursor.rowcount != 1:
                    raise RuntimeError(f"duplicate matrix {matrix_id} changed during circular normalization")
        print(
            f"applied=true converted={len(conversions)} removed_duplicates={len(duplicate_rows)} "
            f"invalidated_candidate_rows={candidate_rows} "
            f"deleted_timing_rows={timing_rows}"
        )
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--apply", action="store_true")
    arguments = parser.parse_args()
    normalize(arguments.database.resolve(), arguments.apply)


if __name__ == "__main__":
    main()
