# Python Review

Last verified: 2026-08-03

Scope: maintained Python wrapper, multiprocessing, sinks, generic JSON input,
and wrapper tests. The native `fracessa_core` extension is reviewed separately
in `PYBIND_REVIEW.md`. No active matrix-verification runner exists. Sequential
timing uses the canonical SQLite database under `testdata/`.

Correctness is ranked before speed. This maintained audit record keeps current open findings, if any, plus validation evidence and
completed decisions that prevent rejected wrapper designs from being repeated.

## Open Findings

None.

## Current Validation State

- All 63 PyFracESSA tests passed against the combined Release native module
  with PyArrow available.
- The former JSON/CSV verification, baseline-generation, subprocess benchmark,
  and JSON-fed Callgrind paths have been removed. Their small replacement timing
  tool keeps one build loaded on a pinned CPU, sizes each matrix/method sample from its matrix-owned calibration, and stores the
  median returned nanoseconds in the canonical SQLite database. The default target is 0.5 seconds.
- Generic JSON input requires its configured key, a row list, and object rows;
  malformed schemas fail instead of silently producing no work.
- `Matrix` validates built-in string/dictionary fields and signed 64-bit
  integer IDs at construction; the JSON loader and native adapter no longer
  coerce float, Boolean, numeric-string, or non-string matrix values.
- The unused input pass-through iterator and its collection imports have been
  deleted.
- `testdata/fracessa_testdata.sqlite3` passes SQLite integrity and foreign-key
  checks and contains 1,072 distinct strategically normalized matrices. Its 780 analyzed rows store 67,875 representatives whose
  multipliers recover 106,401 candidates and 91,950 ESS; 762 baselines are exact or exact-fallback, 18 are unverified fast-only,
  and 292 rows remain catalog-only. Dimensions 2-25 each have at least
  one circular and one non-circular matrix, and every distinct matrix from the
  two published Bomze-Schachinger-Ullrich result tables is present once.
- Its retained timing data have 892 persistent-Pybind median rows spanning Werner, `classic`, current quick-suite, paired
  safe-wrapper, four fast-path experiment panels, the fast-timeout retry, and circular normalization. Catalog-only matrices are
  excluded automatically.
  Historical rows for seven matrices whose stored values changed were removed; retained exact results and every current panel row
  match their expected ESS counts. Report rows
  include dimension, circularity, and the derived lower bound
  `gamma_lower_bound = expected_ess ** (1 / dimension)`.
- Fast calibration has 780 measured durations and 292 cutoff sentinels; safe calibration has 762 measured durations and 310
  sentinels. No field is null, all 1,072 rows have audit timestamps, and whole-matrix fallback counts are 955 null, 45
  `precision_span`, four `equilibration_invalid`, and 68 `equilibration_non_convergence`.
- Sequential and multiprocessing paths use one flat result dictionary;
  CSV, JSON, and Parquet are the only output sinks.
- `run` and `run_multiprocessing` are the only public execution functions;
  both accept one `Matrix` or an iterable and optionally consume directly into
  a sink. `compute_matrix` is the sole public low-level native adapter.
- `run_multiprocessing` keeps the same parameter order as `run` and adds only a
  final optional `MPConfig`; `MPConfig()` defaults to the CPUs available to the
  Python process.
- `RunConfig()` contains analysis options only. `run`, `run_multiprocessing`, and `compute_matrix` require `fast`, `safe`, or
  experimental `test` as their first argument; there is no default method.
- Multiprocessing serializes results inside workers before queueing them, so a
  return-trip serialization failure exits the worker instead of hanging.
- Multiprocessing uses one shared matrix queue and one shared result queue, yields
  completion order without a reordering buffer, and bounds pending work with
  `prefetch_per_worker`.
- Multiprocessing rejects `RunConfig.enable_logging=True` before creating any
  workers; native logging remains available for sequential execution.
- File sinks create output triplets exclusively and raise `FileExistsError`
  without changing existing output when a run ID is reused.
- Sink construction and consumption share one retained rollback stack. Any
  initialization, computation, write, or finalization exception closes the
  result iterator and sink resources, removes the attempt's three new paths,
  and re-raises the original error so the same run ID can be retried.
- Every file sink writes per-matrix metadata to a format-specific JSON sidecar;
  rows with `metadata=None` are omitted.
- JSON output replaces non-finite floats with `null` and enables strict
  encoding, so summary, candidate, and metadata files remain standard JSON.
- Empty CSV, JSON, and Parquet summary/candidate datasets are readable
  with stable schemas.
- Parquet buffers 1,024 rows; a 1,100-result regression produces two row groups
  rather than 1,100.
- Release CI installs `pyarrow` before running the wrapper suite, so retained
  Parquet tests cannot be skipped for a missing dependency.
- Pull requests and `main` run Ubuntu and macOS once; Windows is temporarily
  release-tag only until its dependency installation is fast enough for CI.
  Packaging, artifact upload, and publication remain release-tag only.
- Build-backed native integration tests fail when `fracessa_core` is missing;
  they no longer report the binding boundary as a successful skip.
- Sink finalization attempts every close and propagates the first failure;
  cleanup errors never replace an active computation or write error.
- A 64-matrix spawn-mode run completed, early iterator closure left no workers,
  input serialization failed synchronously, and the return-serialization
  regression passed without hanging.
- Eight concurrent JSON sink constructors produced one intact output triplet
  and seven clean `FileExistsError` results.
