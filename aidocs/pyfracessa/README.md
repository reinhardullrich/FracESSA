# PyFracESSA

The maintained Python API calls the native `fracessa_core` pybind module
in-process.

## Build And Import

Python 3.11 through 3.14 is supported. Install the published package with:

```bash
python -m pip install pyfracessa
```

The native module is built for the selected interpreter and must be imported by that same Python minor version. To build from the
repository instead:

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build -j"$(nproc)"
PYTHONPATH=python python3 your_script.py
```

The internal native loader searches the normal Python import path and
`cpp/build/{,Release,RelWithDebInfo,Debug}`.

## Generated API Documentation

Every production module, class, function, and method has a standard Python
docstring. The built-in `pydoc`, the lightweight `pdoc` generator, and Sphinx
`autodoc` can all read them. To inspect the installed public API without adding
a documentation dependency:

```bash
PYTHONPATH=python python3 -m pydoc fracessa
```

Use `pdoc` for a small standalone HTML API site or Sphinx when the API reference
must be part of a larger authored manual.

## Core Types

```text
Matrix(matrix_id: int, matrix: str, metadata: dict | None = None)
RunConfig(full_support=False, include_candidates=False, enable_logging=False)
MPConfig(workers=<available CPUs>, prefetch_per_worker=128, queue_maxsize=4096,
         start_method="spawn")
```

`Matrix` validates its public contract immediately: `matrix_id` must be a
built-in Python integer in the signed 64-bit range (booleans are rejected),
`matrix` must be a string, and `metadata` must be a dictionary or `None`.

Every computation returns one plain dictionary with `run_id`, `matrix_id`,
`status`, `candidate_count`, `ess_count`, `candidate_structure`, `ess_structure`, `elapsed_ns`, `safe_fallback`, `error_message`,
`candidates`, and `metadata`. `candidates` is a list of plain dictionaries; there are no result-row classes or conversion step.

`safe_fallback` is `None` unless `fast` bypassed double search for the complete matrix. The possible reasons are
`"precision_span"`, `"equilibration_invalid"`, and `"equilibration_non_convergence"`. An inconclusive pivot for one support does
not set this field.

Each candidate dictionary has a nullable `multiplier`: circular matrices return
one bracelet representative with its orbit count, while ordinary candidates use
`None`. The complete candidate and ESS counts and their support-size structures are always returned, independently of
`include_candidates`. That option controls only whether the individual representative rows are included.

A matrix may use full CLI form (`"3#4,13/2,..."`) or values only when
`metadata["dimension"]` is present.

Every execution call requires `"fast"`, `"safe"`, or experimental `"test"` before the matrix; there is no default. Fast and test
apply the exact precision-span gate, convert and equilibrate the complete game once, and then use independent binary64 candidate
implementations. A matrix-wide preparation failure uses safe exact candidate checks instead. Safe bypasses every floating-point
candidate procedure. Matrix input always uses the validating native parser.

`elapsed_ns` is always the native analyzer duration measured with a monotonic
clock. There is deliberately no computation timeout.

## Timing Database

The repository timing tool runs sequentially on one selected Linux CPU and
stores native analyzer durations in the existing test-data database. It does
not call `run_multiprocessing`.

One invocation measures one build. The generated session is printed; pass that
same session to later invocations to group a moving baseline and temporary
build without loading two `fracessa_core` modules into one interpreter:

```bash
PYTHONPATH=python python3 -m fracessa.timing run \
  --backend pybind --module-dir cpp/build \
  --build-label main --source-ref main --revision "$(git rev-parse HEAD)" \
  --method fast --method safe \
  --cpu 2 --comment "before candidate change"

PYTHONPATH=python python3 -m fracessa.timing report latest --baseline main
```

`source_ref` records a moving name such as `main`; `revision` records the exact
commit measured, and the module or executable SHA-256 is captured automatically.
The default selection is the `small` matrix class; `--size-class` also accepts `medium`, `large`, `super_large`, and `all`. Each
matrix/method starts with
the `fast_calibration_ns` or `safe_calibration_ns` value in its database row. These values are maintained with local-only database
helpers that are intentionally excluded from public clones and releases. The timing runner uses
`ceil(target / calibration)` samples for positive calibrations; a calibration at or above the target or a `-1` timeout
sentinel chooses one run. The default target is 0.5 seconds, and a missing calibration is an error. The stored result is the median
returned native duration. The Pybind module stays loaded for every selected method and matrix in that invocation. Python wall time
is recorded as metadata only and does not select the sample count or result. Use repeated `--matrix-id`, `--size-class`, and
`--target-seconds` to change the selection or target duration.

Current Pybind builds supply nanoseconds directly. For a legacy CLI whose
no-flag method corresponds to today's fast method and whose second output line
is seconds, use `--fast-default --cli-unit s`. A safe-by-default legacy CLI
instead uses `--safe-default`. CLI runs start one child process per sample and
must not be mixed with persistent-Pybind microbenchmarks. The timing adapter
maps `fast` and `safe` onto the old mode or Boolean interfaces when measuring
historical binaries.
Reports include the iteration count and actual measured wall duration, compare
observed ESS counts with the expected database count, and retain mismatching
unsafe timings visibly. Each result row also includes the matrix dimension,
circular-symmetry flag, and the paper-style lower bound
`gamma_lower_bound = expected_ess ** (1 / dimension)`. The bound is derived from
the canonical expected count and matches the generated `matrices.gamma_lower_bound` column.

Status codes:

| Code | Name |
|---:|---|
| 0 | `OK` |
| 1 | `PARSE_ERROR` |
| 4 | `EXEC_ERROR` |
| 255 | `INTERNAL_ERROR` |

Every rejected safe-parser input returns `PARSE_ERROR` with the parser's
detailed diagnostic in `error_message`. The parser throws
`std::invalid_argument`; Pybind catches it and does not write to `stderr`.

## Execution API

```python
from fracessa import Matrix, RunConfig, run

matrix = Matrix(3, "3#4,13/2,1/2,5,11/2,3")
result = run("safe", matrix, RunConfig(include_candidates=True), run_id="example")
print(result["ess_count"])
```

`compute_matrix(method, matrix, config, run_id) -> dict` is the sole public low-level
wrapper call into the native module. Native-module loading is internal.

There are two public execution functions:

- `run(method, matrices, config=None, run_id=None, sink=None)`
- `run_multiprocessing(method, matrices, config=None, run_id=None, sink=None,
  mp_config=None)`

Both accept either one `Matrix` or an iterable of matrices. Without a sink,
one matrix returns one dictionary and an iterable returns an iterator. With a
sink, execution is eager and returns the number of matrices written.

## Multiprocessing Batch API

Use the batch iterator for lists, generators, and large streams. One worker
computes one matrix at a time; parallelism is across matrices.

```python
from fracessa import (
    MPConfig,
    Matrix,
    RunConfig,
    run_multiprocessing,
)

matrices = [
    Matrix(1, "2#0,1,0"),
    Matrix(2, "3#4,13/2,1/2,5,11/2,3"),
]

if __name__ == "__main__":
    for result in run_multiprocessing(
        "fast",
        matrices,
        config=RunConfig(include_candidates=False),
        mp_config=MPConfig(workers=8),
    ):
        print(result["matrix_id"], result["ess_count"])
```

When omitted, `mp_config` uses all CPUs available to the Python process. Its
other fields control the pending-matrix window, result-queue capacity, and
process start method. `RunConfig` controls matrix analysis; `MPConfig` controls
only process scheduling.

Submission is bounded to
`min(queue_maxsize, workers * prefetch_per_worker)` outstanding matrices. The
parent puts one matrix at a time onto one shared matrix queue, every worker takes
from that queue, and all workers return dictionaries through one shared result
queue. Results are yielded as workers finish; there is no ordered mode or IPC
batching. Matrices are serialized before submission, worker exits are detected
while waiting, and closing the iterator early cancels its workers. The queue
runner is internal.

Native logging is sequential-only. Passing `RunConfig(enable_logging=True)` to
`run_multiprocessing()` raises `ValueError` before any worker is created because
independent processes cannot safely rotate the shared native log file. Direct
logging-enabled calls from Python threads are serialized by the native binding;
non-logging calls remain concurrent.

The `if __name__ == "__main__":` guard is required in scripts when using the
default `spawn` start method.

## Input Loader

- `load_matrices_from_json(path, ...)` accepts either a top-level list or
  `{"matrices": [...]}`.

Top-level objects must contain the configured matrix key, its value must be a
list, and every row must be an object. Malformed schemas raise `ValueError`.

Rows without a `dimension#` matrix prefix must provide a built-in integer
dimension field. The loader rejects floats, booleans, and strings for integer
fields and rejects non-string matrix values instead of coercing them.

## Output Sinks

- `CsvSink`: summary and candidate CSV files plus metadata JSON.
- `JsonSink`: summary, candidate, and metadata JSON arrays.
- `ParquetSink`: Parquet summary/candidate files plus metadata JSON; requires
  `pyarrow`.
- `create_sink(kind, output_dir, run_id)`: constructs a named sink.

Every file sink accepts `(summary_path, candidates_path, metadata_path=None)`.
When omitted for direct construction, the metadata path defaults to
`<summary_path>.metadata.json`. `create_sink()` instead uses format-specific
names such as `<run_id>_csv_metadata.json` so each output is unambiguous.

Metadata sidecars contain streamed rows of `run_id`, `matrix_id`, and the
original per-matrix metadata dictionary. Results with `metadata=None` add no
metadata row. All three paths are created exclusively; if any exists,
construction raises `FileExistsError` without overwriting it.

Every JSON file is standards-compliant. Non-finite floating-point approximations
(`NaN` or positive/negative infinity) are written as JSON `null`; unsupported
metadata values still raise and trigger the normal triplet rollback.

Construction and consumption retain one rollback stack for the complete run. If
initialization, computation, writing, or finalization raises, the result
iterator and sink resources are closed, only paths created by that attempt are
removed, and the original error is re-raised. Cleanup errors do not replace the
active error, and the same run ID can then be retried.

Empty outputs retain readable stable schemas. Parquet aggregates up to 1,024
result or candidate rows before writing a row group and flushes the final
partial group during `close()`. Closing a sink attempts every underlying
finalization and propagates the first failure.

For large runs, consume the result iterator continuously or pass `sink=` to one
of the two execution functions so results do not accumulate in memory. Disable
candidates when they are not needed; candidate payloads can dominate output
volume.

## Tests

```bash
PYTHONPATH=python python3 -m unittest discover \
  -s python/tests -p 'test_*.py'
```

The suite includes unit tests plus mandatory native single-process,
multithreaded-logging, and multiprocessing integration tests. Build `fracessa_core` first with the CMake commands above; a missing
module is a test failure.
Ordinary pushes and pull requests run no GitHub Actions. The manually started
release workflow builds and tests all supported platforms before creating its
tag, publishing the standalone archives, and publishing `pyfracessa`.
