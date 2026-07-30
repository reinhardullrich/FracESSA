# FracESSA Test Data

`fracessa_testdata.sqlite3` is the canonical store for test matrices and their
complete expected candidate results.

The current snapshot contains 85 matrices and 49,155 stored candidate
representatives. Their multipliers represent 86,150 candidates and 83,375 ESS.
It contains each distinct matrix from Tables 1 and 2 of the
Bomze-Schachinger-Ullrich ESS-growth paper exactly once. IDs 18 and 26 hold the
exact published Table 1 matrices that replaced same-property alternatives;
IDs 80-90 hold the previously missing Table 2 base and constructed matrices.
Redundant alternatives formerly at IDs 12 and 21 were removed. IDs 56-66 are
staged complete-multipartite many-ESS benchmark matrices. IDs 67-79 are
deterministic random-integer coverage matrices; together with the existing rows,
every dimension from 2 through 25 has at least one circular and one non-circular
matrix. No SQLite matrix suite is currently wired into CTest.

The timing table contains one complete current-build Pybind session: all 85
matrices in both unsafe and exact modes, for 170 adaptive measurements with no
ESS-count mismatch.

## Tables

### `matrices`

Each row stores one exact matrix input and its summary:

- `matrix_id`: stable FracESSA verification ID.
- `dimension`: number of strategies, from 1 through 63.
- `size_class`: `small` for dimensions 1-8, `medium` for 9-16, and
  `large` for 17-63.
- `is_cs`: 1 for compact circular-symmetric input, otherwise 0.
- `matrix`: the exact comma-separated input values, without the `n#` prefix.
- `candidate_count` and `ess_count`: complete weighted baseline counts.
- `candidate_structure`: JSON object mapping support size to weighted candidate
  count.
- `ess_structure`: JSON object mapping support size to weighted ESS count.
- `origin`: where the matrix came from and why it was retained.
- `tags`: JSON array of short qualitative categories such as
  `"numerical_edge"` or `"support_frontier"`.

There is deliberately no matrix-name field. The stable ID identifies a matrix;
`origin` records the useful human context: source or construction, purpose, and
the property the matrix is meant to exercise. Where old provenance was never
recorded, the origin says so rather than inventing a history.

For example, `{"1":8,"4":2}` means eight support-one results and two
support-four results. Empty ESS structure is stored as `{}`.

### `candidates`

Each row mirrors one candidate output row. A circular row stores only the
smallest support in its rotation/reflection orbit, and its non-null `multiplier`
is the number of distinct supports represented. A non-circular row has a null
multiplier. Exact fractions and vectors remain text; `payoff_double` is also
text so database reads cannot alter formatting. `(matrix_id, candidate_id)` is
the primary key, and a support may occur only once for a matrix.

Fixed facts already represented by columns, including size, circular symmetry,
counts, and support-size structures, are not duplicated in `tags`.

### `timings`

Each row is one sequential analyzer timing measurement for one matrix. A
session may contain several builds, but each build is measured by a separate
runner invocation. Rows record the machine and pinned CPU, human build label,
moving source reference such as `main`, immutable revision, binary SHA-256,
Pybind or CLI backend, numerical mode, target and measured wall durations,
iteration count, average native duration in nanoseconds, observed ESS count,
and an optional comment.

The observed count remains separate from the expected count in `matrices`, so a
report can expose unsafe-mode mismatches without hiding or rejecting their
timings. The report derives the Bomze-Schachinger-Ullrich exponential-growth
lower bound `expected_ess ** (1 / dimension)` and prints it with dimension and
circularity; this value is not stored in the database. Old CLI builds simply
have no `safe` rows.

## Scope

`python -m fracessa.timing` reads matrices from this database and writes timing
observations back to `timings`. It is deliberately a sequential, Linux
CPU-affinity runner, not part of the multiprocessing wrapper. One pilot at or
above the target is the complete measurement. Faster cases use the pilot as a
warmup and choose enough measured iterations for about the requested wall
duration, one second by default; the stored result is their average native
duration. The tool supports the current Pybind timer and legacy CLI timers whose
output unit is supplied on the command line. No active matrix-verification
runner is wired into CTest yet.

The schema is defined in `schema.sql`. The C++ runtime does not read this
database; the timing tool uses Python's standard `sqlite3` module.

## Integrity

Basic database checks are:

```bash
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA integrity_check;'
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA foreign_key_check;'
```
