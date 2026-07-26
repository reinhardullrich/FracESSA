# Python Review

Last verified: 2026-07-26

Scope: maintained Python wrapper, multiprocessing, sinks, verification scripts,
speed benchmark, and wrapper tests. The archived reference executables are
treated as data, not reviewed as source.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Correctness

### P1: `MPQueueRunner` can deadlock or wait for impossible results

The exposed low-level API encourages a producer to call `submit()` repeatedly
before consuming results. Bounded input and output queues can then block each
other at `python/wrapper_v1/mp.py:101` and `python/wrapper_v1/mp.py:43`.
`_get_result_item()` blocks forever at `python/wrapper_v1/mp.py:134` if a worker
dies without publishing.

The API also accepts arbitrary `metadata`, but a non-pickleable value can fail
inside the queue feeder after `submit()` has counted the job. Likewise,
`iter_results(expected_results=...)` accepts a count larger than the number
submitted. Both paths leave the consumer waiting for a result that cannot
arrive.

`run_jobs_mp()` interleaves bounded submission and consumption and is safer, but
the raw queue API still fails its on-the-fly generation use case. Early consumer
exit can also enter `close_input()` while workers are blocked publishing output.

Required outcome: make submission and consumption interleave in the documented
low-level workflow, reject impossible expected-result counts synchronously, and
detect dead workers while waiting. The existing five-second join applies during
shutdown/cancellation; it is not a normal per-matrix timeout and should not be
described as one.

### P1: Baseline regeneration can publish failed output

`python/verification/create_baselines.py:338` reports expected-ESS mismatches as
warnings, still counts those runs as successful, writes both baseline files
directly at `python/verification/create_baselines.py:355`, and exits zero even
after execution failures.

Required outcome: require every execution and ESS expectation to pass, write
both temporary files completely, then replace each destination. Any mismatch,
parse failure, or execution failure must return nonzero and leave the existing
baselines untouched.

### P1: The correctness verifier silently accepts malformed field values

`parse_bool()` at `python/verification/ctest_verify_matrix.py:28` maps every
unknown string to `False`. `normalize_candidate_row()` and
`parse_cli_candidates()` map an invalid floating payoff to `0.0` at
`python/verification/ctest_verify_matrix.py:51` and
`python/verification/ctest_verify_matrix.py:83`. A corrupted baseline and a
malformed CLI row can therefore compare equal for false/zero fields.

The candidate parser also silently skips short rows at
`python/verification/ctest_verify_matrix.py:81` rather than failing the test.

Required outcome: parse strict accepted representations and raise on every
unknown, missing, or malformed field. Correctness fixtures must fail closed.

### P2: Generated output IDs can collide and truncate prior runs

`new_run_id()` has one-second resolution at `python/wrapper_v1/core.py:15`.
`python/run_matrices.py:237` uses the same second-resolution timestamp. Two
same-prefix runs in one second produce identical paths, and all file sinks open
with truncating write semantics.

Required outcome: use collision-resistant IDs or exclusive file creation.

### P2: Empty output is inconsistent across sink backends

`CsvSink` creates an empty candidate file without a header when no candidate is
ever written. `ParquetSink` creates no files for an empty run and no candidate
file for a valid run with zero candidates because writers are created lazily at
`python/wrapper_v1/sinks_parquet.py:35` and
`python/wrapper_v1/sinks_parquet.py:45`. JSON and Arrow produce readable empty
datasets.

Reproduced with installed pyarrow 23.0.1: closing an unused Parquet sink leaves
both paths absent; writing one zero-candidate result creates only the summary
path.

Required outcome: every requested sink must produce readable summary and
candidate datasets with stable schemas, including zero-row cases.

### P2: Sink finalization can hide errors or skip remaining sinks

`StreamSink.close()` catches and discards every flush exception at
`python/wrapper_v1/sinks.py:39`. Buffered output can therefore fail only during
the final flush while `run_many_to_sink()` or `run_jobs_mp_to_sink()` still
reports a successful result count.

`MultiSink.close()` stops at the first close exception at
`python/wrapper_v1/sinks.py:54`, so later sinks may never flush or close.

Required outcome: attempt every sink close, then propagate the first failure.
If broken-pipe suppression is wanted for a CLI, handle that explicitly at the
CLI boundary rather than in the reusable sink.

### P3: Legacy Python dimension inference disagrees with the CLI

`python/fracessa_py.py:41` infers a square dimension from value count. A normal
non-circular CLI matrix contains `n*(n+1)/2` upper-triangular values, so valid
values-only input cannot be inferred correctly.

The active speed script always supplies a dimension, so this is a legacy API
contract bug rather than a current benchmark failure. Required outcome: infer
triangular value counts correctly or delete auto-inference. Migrating the speed
script is tracked separately below.

## Speed

### P1: Arrow and Parquet write one summary row per batch

`ArrowSink.write_result()` creates a one-row RecordBatch at
`python/wrapper_v1/sinks_arrow.py:58`; `ParquetSink.write_result()` creates a
one-row Table and row group at `python/wrapper_v1/sinks_parquet.py:31`. This
adds container overhead for every result. A fresh 100-result probe produced 100
Arrow batches and 100 Parquet row groups.

Required outcome: buffer fixed-size row batches and flush the final partial
batch on close.

### P1: The speed benchmark measures cold-process and candidate work

`python/run_matrices.py:22` still uses the legacy subprocess wrapper. Each
repeat launches a fresh executable and requests all candidates at
`python/run_matrices.py:52`. The reported C++ interval starts after input parsing
and ends before output serialization, so it does not include subprocess launch,
parsing, or text output. It does include cold-process analyzer execution and
candidate collection.

Timing is truncated to whole microseconds at `cpp/src/main.cpp:48`, and
`python/run_matrices.py:80` keeps the minimum sample. Both choices distort the
very small matrices most in need of repeated measurement.

Required outcome: warm up the in-process pybind API, disable candidates for the
core benchmark, and report median plus dispersion. Benchmark candidate output
separately.

### P2: Multiprocessing `chunk_size` does not batch IPC

`chunk_size` only changes the submitted-job window in
`_max_buffered_results()` at `python/wrapper_v1/mp.py:46`. Every matrix and
result still incurs one queue operation and one pickle/unpickle cycle. For the
project's very small matrices, IPC can exceed computation time.

Ordered mode also counts only yielded results when reopening the submission
window at `python/wrapper_v1/mp.py:173`; one early long-running matrix can hold
many completed later results and leave workers idle.

Required outcome: use the existing `ordered_results=False` mode when throughput
matters. Rename `chunk_size` to reflect prefetching, or implement real IPC
chunks only after a benchmark shows serialization/queue overhead is material;
do not unbound pending ordered results merely to keep workers busy.

### P2: Benchmark comparisons lack build and machine provenance

Result metadata records mode and wall time but not hostname, CPU, affinity,
compiler, flags, Git commit, or dirty-worktree state at
`python/run_matrices.py:318`. The script nevertheless prints percentage changes
against historical files that may come from another machine or build.

Required outcome: record enough provenance to identify a computer-plus-build
run and refuse or visibly qualify comparisons across incompatible provenance.

### P3: Disabling timing does not remove native clock overhead

`RunConfig.include_timing=False` only replaces the returned value with zero at
`python/wrapper_v1/core.py:139`. The pybind implementation always calls the
clock before and after every matrix at `cpp/src/pybind_module.cpp:109` and
`cpp/src/pybind_module.cpp:111` because the option is not passed across the
binding.

Required outcome: benchmark the complete pybind call with and without the two
clock reads before extending the native API. The current public documentation
only promises to suppress the returned value, so this is an optimization
opportunity rather than a behavior bug.

## Documentation

### P3: The public README labels the speed script as verification

`README.md:34` calls `./python.sh` "Batch Verification", but the script is a
speed benchmark and only warns on ESS mismatches. Correctness verification is
part of `./test.sh`/CTest.

Required outcome: rename that README section and keep correctness claims tied to
CTest.

## Current Validation State

- Python syntax/import path: successful through the normal build/test flow.
- Wrapper unittests: 23/23 passed with the native module present.
- Native single-process and multiprocessing integration tests passed.
- Direct probes reconfirmed the raw queue submission hang, same-second run-ID
  collision, verifier coercion, and suppressed stream flush failures.
- Pyarrow 23.0.1 probes reconfirmed missing empty Parquet datasets and one batch
  or row group per result.
- Python correctness tests do not override the four C++ matrix regressions; the
  full CTest state is recorded in `CPP_REVIEW.md`.
