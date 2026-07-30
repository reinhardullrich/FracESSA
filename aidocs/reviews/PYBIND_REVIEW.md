# Pybind Review

Last verified: 2026-07-31

Scope: the native `fracessa_core` extension boundary in
`cpp/src/pybind_module.cpp`, including argument/result conversion, native status
codes, GIL handling, timing exposed through the binding, and binding-specific
integration behavior. The analyzer core belongs to `CPP_REVIEW.md`; Python
orchestration belongs to `PYTHON_REVIEW.md`.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Simplicity

### P3: `success` duplicates the status code

`NativeResult::success` at `cpp/src/pybind_module.cpp:45` is always equivalent
to `status == STATUS_OK`; the Python wrapper and sinks carry both fields without
using `success` for control flow.

Required outcome: delete the duplicate field through the native dictionary,
wrapper result, sink schemas, tests, and documentation. Callers can compare the
status when needed.

## Current Validation State

- The combined Release build passed all 11 C++/CLI tests and all 53 PyFracESSA
  tests with the native module and PyArrow available.
- Every safe-parser rejection returns the single `PARSE_ERROR` status plus the
  shared parser's detailed diagnostic without reparsing or writing to `stderr`.
- The native result preserves the analyzer's `size_t` ESS count through Pybind;
  Python receives it without a 32-bit narrowing conversion.
- Matrix IDs use signed 64-bit values from CLI/Python through the analyzer;
  boundary tests cover the maximum value.
- Native candidate conversion has an exact 11-field value-and-type contract
  regression. Ordinary rows use `multiplier=None`; one circular regression
  returns a representative with `multiplier=5` and a weighted ESS count of 5.
- The binding exposes verified search by default, explicit unsafe search, and
  exact search; exact and unsafe together return `EXEC_ERROR`.
- Native single-process and multiprocessing integration tests pass.
- One process-wide native mutex serializes logging-enabled analyzer calls from
  Python threads while non-logging calls remain concurrent.
- `compute_matrix()` always measures with `std::chrono::steady_clock` and
  returns integer nanoseconds as `elapsed_ns`; no timing-suppression option
  remains.
- Shared-parser and full CTest state remains recorded in `CPP_REVIEW.md`.
