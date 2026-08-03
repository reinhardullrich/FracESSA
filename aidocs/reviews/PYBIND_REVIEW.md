# Pybind Review

Last verified: 2026-08-03

Scope: the native `fracessa_core` extension boundary in
`cpp/src/pybind_module.cpp`, including argument/result conversion, native status
codes, GIL handling, timing exposed through the binding, and binding-specific
integration behavior. The analyzer core belongs to `CPP_REVIEW.md`; Python
orchestration belongs to `PYTHON_REVIEW.md`.

Correctness is ranked before speed. This maintained audit record keeps current open findings, if any, plus validation evidence and
completed boundary decisions that prevent contract drift.

## Open Findings

None.

## Current Validation State

- The combined Release build passed all 10 C++/CLI tests and all 63 PyFracESSA
  tests with the native module and PyArrow available.
- Results use the integer `status` plus `error_message` contract without a redundant success Boolean.
- Every safe-parser rejection returns the single `PARSE_ERROR` status plus the
  shared parser's detailed diagnostic without reparsing or writing to `stderr`.
- The native result preserves the analyzer's `size_t` ESS count through Pybind;
  Python receives it without a 32-bit narrowing conversion.
- Matrix IDs use signed 64-bit values from CLI/Python through the analyzer;
  boundary tests cover the maximum value.
- Native candidate conversion has an exact 11-field value-and-type contract
  regression. Ordinary rows use `multiplier=None`; one circular regression
  returns a representative with `multiplier=5` and a weighted ESS count of 5.
- The binding requires `fast`, `safe`, or experimental `test` before the matrix; there is no default or compatibility alias, and
  unknown methods return `EXEC_ERROR`.
- Native single-process and multiprocessing integration tests pass.
- One process-wide native mutex serializes logging-enabled analyzer calls from
  Python threads while non-logging calls remain concurrent.
- `compute_matrix()` always measures with `std::chrono::steady_clock` and
  returns integer nanoseconds as `elapsed_ns`; no timing-suppression option
  remains.
- Shared-parser and full CTest state remains recorded in `CPP_REVIEW.md`.
