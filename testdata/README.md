# FracESSA Test Data

`fracessa_testdata.sqlite3` is the canonical store for test matrices and their
complete expected candidate results.

The current snapshot contains 90 matrices and 49,161 stored candidate
representatives. Their multipliers represent 86,156 candidates and 83,381 ESS.
It contains each distinct matrix from Tables 1 and 2 of the
Bomze-Schachinger-Ullrich ESS-growth paper exactly once. IDs 18 and 26 hold the
exact published Table 1 matrices that replaced same-property alternatives;
IDs 80-90 hold the previously missing Table 2 base and constructed matrices.
Redundant alternatives formerly at IDs 12 and 21 were removed. IDs 56-66 are
staged complete-multipartite many-ESS benchmark matrices. IDs 67-79 are
deterministic random-integer coverage matrices; together with the existing rows,
every dimension from 2 through 25 has at least one circular and one non-circular matrix. IDs 45-47 preserve the retired normalized
heuristic and verified-search regression history. IDs 91-93 preserve one exact false rejection for each former fast per-support
candidate condition: the pivot cutoff, probability sign, and outside-payoff margin. Current fast uses its small-pivot and
precision-span fallbacks to recover them. No complete SQLite matrix suite is currently wired into CTest.

The timing table retains four CPU-2 persistent-Pybind build families: Werner fast and safe, the July 27 pre-refactor GitHub build
renamed `classic`, and the paired safe builds immediately before and after the C++ FLINT-wrapper extraction. The paired builds use
one-second native medians for 77 matrices of dimension at least 3; IDs 47, 65-66, and 90 are excluded by the 30-second rule. All
paired ESS counts match. Build label, revision, and binary hash identify every stored build.

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
Pybind or CLI backend, search method in the historically named `mode` column, target and measured wall durations,
iteration count, median native duration in nanoseconds, observed ESS count,
and an optional comment.

The observed count remains separate from the expected count in `matrices`, so a
report can expose fast-method mismatches without hiding or rejecting their
timings. The report derives the Bomze-Schachinger-Ullrich exponential-growth
lower bound `expected_ess ** (1 / dimension)` and prints it with dimension and
circularity; this value is not stored in the database.

## Scope

`python -m fracessa.timing` reads matrices from this database and writes timing
observations back to `timings`. It is deliberately a sequential, Linux
CPU-affinity runner, not part of the multiprocessing wrapper. One Pybind process
stays open for all selected methods and matrices in a build. The first returned
C++ duration chooses `ceil(target / duration)` total samples and remains part
of the sample; a duration at or above the target chooses one run. The stored
result is the median returned `elapsed_ns`. Wall time is recorded as metadata
but does not choose the sample count or result. The tool also supports legacy
CLI timers whose output unit is supplied on the command line, but those
fresh-process rows must not be mixed with persistent-Pybind microbenchmarks. No
active matrix-verification runner is wired into CTest yet.

The schema is defined in `schema.sql`. The C++ runtime does not read this
database; the timing tool uses Python's standard `sqlite3` module.

## Integrity

Basic database checks are:

```bash
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA integrity_check;'
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA foreign_key_check;'
```
