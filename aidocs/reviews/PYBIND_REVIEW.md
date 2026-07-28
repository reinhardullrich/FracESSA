# Pybind Review

Last verified: 2026-07-27

Scope: the native `fracessa_core` extension boundary in
`cpp/src/pybind_module.cpp`, including argument/result conversion, native status
codes, GIL handling, timing exposed through the binding, and binding-specific
integration behavior. The analyzer core belongs to `CPP_REVIEW.md`; Python
orchestration belongs to `PYTHON_REVIEW.md`.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Correctness

### P2: Pybind reparses parser failures and can misclassify them

The shared parser returns only `bool`. When it fails, pybind reparses the
same input in `infer_parse_status()` at `cpp/src/pybind_module.cpp:59`.
That second parser covers only structural and dimension errors and otherwise
returns `INVALID_VALUE_COUNT`, so an invalid rational can receive the wrong
public status at `cpp/src/pybind_module.cpp:102`.

Required outcome: return one typed parse status from the shared parser to both
CLI and pybind callers, then remove `infer_safe_parse_status()`.

### P3: Pybind narrows the ESS count

The analyzer stores `ess_count_` as `size_t`, but `NativeResult::ess_count` is
`int` at `cpp/src/pybind_module.cpp:41` and receives an explicit narrowing cast
at `cpp/src/pybind_module.cpp:114`. The search space can exceed `INT_MAX`, while
Python integers do not require this truncation.

Required outcome: preserve `size_t` or `uint64_t` through the native result.

## Speed

### P3: Disabling timing does not remove native clock overhead

`RunConfig.include_timing=False` only replaces the returned value with zero at
`python/wrapper_v1/core.py:139`. The pybind implementation always calls the
clock before and after every matrix at `cpp/src/pybind_module.cpp:109` and
`cpp/src/pybind_module.cpp:111` because the option is not passed across the
binding.

Required outcome: benchmark the complete pybind call with and without the two
clock reads before extending the native API. The public wrapper contract only
promises to suppress the returned value, so this is an optimization opportunity
rather than a behavior bug.

## Current Validation State

- The current Release build produces an importable `fracessa_core` module.
- Wrapper unittests that exercise the native module pass 23/23.
- Native single-process and multiprocessing integration tests pass.
- `compute_matrix()` still always measures elapsed time and returns the narrowed
  native ESS count described above.
- Shared-parser and full CTest state remains recorded in `CPP_REVIEW.md`.
