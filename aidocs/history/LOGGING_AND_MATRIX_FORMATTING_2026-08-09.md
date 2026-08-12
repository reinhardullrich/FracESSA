# Logging And Matrix Formatting

Status: historical implementation record from 2026-08-09. The output contract remains current, while filenames, test counts, and
verification details below describe that implementation snapshot rather than the latest repository inventory.

## Goal

Make the optional native log easy for humans to read and remove logging mechanics from the ESS algorithm without changing the
normal runtime path, numerical decisions, CLI result format, or Python result format.

The log answers five questions:

1. What matrix and search method were run?
2. Which support cardinality is being searched?
3. Which exact candidates were found?
4. Why was each candidate accepted or rejected as an ESS?
5. What were the final candidate and ESS counts?

## Problems Addressed

Before this redesign:

- `fracessa.cpp` directly constructed the sink, printed run separators, formatted the input matrix, tracked support-size progress,
  and wrote candidate CSV.
- `check_stability.cpp` repeated stability assignment and `logger_->info()` calls for every decision path.
- Fraction and integer matrices used a fixed width of 12, which failed as soon as an exact value was longer than 12 characters.
- Candidate supports were printed primarily as bitstrings, which became difficult to read and uselessly long for multiword inputs.
- The diagnostic log duplicated machine-readable candidate CSV and therefore forced logging-only runs to materialize and retain all
  candidate vectors.
- Embedded multiline matrix messages timestamped only the heading, while the rows had no visual structure or row labels.

## Output Contract

Keep each destination responsible for one representation:

- CLI stdout and Python return values remain the machine-readable result contract.
- CSV, JSON, and Parquet sinks remain the persistent candidate-data formats.
- `log/fracessa.log` is a human diagnostic log only.

The diagnostic log no longer contains the semicolon-separated candidate header and rows. `--candidates` remains the way to
request complete candidate data. This separation lets `--log` stop changing candidate-vector materialization and storage.

## Minimal Design

The implementation adds no logger interface, virtual calls, asynchronous logging, event hierarchy, or dependency. `basic_fracessa`
is the only logging owner and has six private semantic helpers:

```cpp
void start_log(search_method method, bool is_cs, std::int64_t matrix_id);
void log_support_size(size_t support_size);
void log_candidate();
void log_reduced_b(const SupportMask& outside_best_replies);
void set_stability_result(bool is_ess, std::string_view reason);
void finish_log();
```

The existing nullable `std::shared_ptr<spdlog::logger>` remains sufficient. Forward-declare `spdlog::logger` in `fracessa.hpp` and
include the concrete spdlog sink only in `fracessa.cpp`. No logger class or new source file is needed.

`log_support_size()` owns the previous support size and replaces the same four-line block in all four generator branches:

```cpp
log_support_size(support_size);
```

`set_stability_result()` performs the state update and diagnostic together:

```cpp
set_stability_result(false, "F_not_copos_nonpositive_diagonal");
return;
```

This replaces the repeated assignment, reason assignment, and conditional log call in every stability path.

The sink is created only when logging is requested, and the run header is emitted after fast/test conversion and whole-matrix fallback
selection. The header can then distinguish `requested=fast` from `effective=safe` and include the exact fallback reason. Use the
nullable logger pointer as the sole enabled flag; a second `conf_with_log_` boolean is unnecessary.

## Human Format

The decorative `*#*#*#` banner was replaced with one searchable run boundary:

```text
==============================================================================
run started: matrix_id=2296 requested=safe effective=safe dimension=9 circular=false
==============================================================================
```

Exact matrices render as bracketed rows. Values are not truncated or padded to the longest arbitrary-precision
fraction:

```text
game matrix (3 x 3)
  0: [0, 1/3, -12/7]
  1: [1/3, 0, 5]
  2: [-12/7, 5, 0]
```

This is shorter and more robust than dynamic padding. Both exact matrix wrappers expose `to_pretty_string()`. Their short
implementations remain separate because rational and FLINT integer entry conversion are different; a generic formatter would save
little code.

Supports are displayed as strategy-index sets rather than bitstrings:

```text
solved candidate representative
  support:              {1, 4}
  extended:             {0, 1, 4, 7}
  reference:            1
  outside best replies: {0, 7}
```

This is the exact support solved before symmetry expansion. For circular games, one logged representative may generate several
stored representatives and each may cover multiple rotations or reflections. Final candidate and ESS counts are therefore labelled
as weighted totals.

The stable reason code remains in the output, with an explicit result:

```text
stability: rejected [F_not_copos_nonpositive_diagonal]
```

The reduced B matrix is printed only on the path that constructs it. Each run ends with:

```text
run finished: weighted candidates=8 weighted ESS=1
```

## Runtime Behavior

The disabled path must remain equivalent to the current code:

- one predictable null-logger check at each existing logging point;
- no formatting, allocation, matrix conversion, or string construction when logging is disabled;
- no new branch inside numerical elimination, bit scanning, or copositivity kernels.

With candidate CSV removed from the diagnostic log:

- reserve `candidates_` only for `with_candidates`;
- ask the exact solver to materialize probability vectors only for `with_candidates`;
- retain and sort circular candidate rows only for `with_candidates`;
- log the exact support and stability decision directly without storing the candidate row.

This makes logging cheaper while leaving the non-logging hot path unchanged.

## Implemented Files

- `cpp/include/fracessa/fracessa.hpp`: declares semantic helpers, retains the nullable logger, and removes logging-specific candidate
  retention state.
- `cpp/src/fracessa.cpp`: centralizes logger creation, run boundaries, support progress, candidate context, and final summary; removes
  candidate CSV from the diagnostic log.
- `cpp/src/check_stability.cpp`: routes every decision through `set_stability_result()` and moves reduced-B diagnostics behind the
  helper.
- `cpp/include/linalg/matrix_fraction.hpp`: replaces fixed-width matrix logging with bracketed exact rows.
- `cpp/include/linalg/matrix_integer.hpp`: applies the same visible format while retaining correct `flint_free()` ownership.
- `cpp/tests/cli_parser_blackbox.py`: verifies the human run boundary, matrix rows, support-set notation, stability reason, multiword
  support formatting, and final summary.
- `aidocs/pyfracessa/README.md`: documents that the native log is diagnostic rather than a candidate-data sink.

No CLI argument, Python configuration field, database schema, candidate schema, or CMake target needs to change.

## Verification Record

- The Release build passed all 10 C++/CLI tests and all 68 Python tests.
- Black-box coverage verifies the run boundary, matrix rows, support-set notation, stability reason, multiword support formatting,
  final summary, absence of candidate CSV in the diagnostic log, and fast-to-safe fallback labels.
- A manual safe run verified exact scaled-reduced-$B$ output.
- No permanent logging benchmark was added because disabled logging performs no formatting or string construction.

## Explicit Non-Goals

- No asynchronous logger.
- No multiprocessing logging support.
- No per-run log-file configuration.
- No JSON log format.
- No generic logging framework.
- No truncation policy for very large matrices in the first patch.

Add one of these features only after a concrete need appears. The implemented patch is strictly a human-readability and separation
cleanup.
