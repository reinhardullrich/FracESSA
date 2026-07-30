# FracESSA Test Data

`fracessa_testdata.sqlite3` stores verification matrices and their complete
candidate baselines. It is an initial migration snapshot; the existing Python
verification files remain active until their consumers are deliberately
switched to this database.

The current snapshot contains 63 matrices and 29,114 stored candidate
representatives. Their multipliers represent 65,962 candidates and 63,369 ESS.
The first 52 matrices mirror the active verification files. IDs 56-66 are
staged complete-multipartite many-ESS benchmark matrices and are not active
CTest fixtures yet.

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

## Scope

Benchmark runs are intentionally not represented yet. Their environment and
result schema will be decided separately, because timing is an observation tied
to a machine and binary rather than a property of a matrix.

The schema is defined in `schema.sql`. The C++ runtime does not read this
database; verification tooling can use Python's standard `sqlite3` module.

## Integrity

Basic database checks are:

```bash
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA integrity_check;'
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA foreign_key_check;'
```
