# Python Wrapper

The maintained Python API calls the native `fracessa_core` pybind module
in-process. It does not use the legacy subprocess wrapper
`python/fracessa_py.py`.

## Build And Import

From repository root:

```bash
./build.sh
PYTHONPATH=python python3 your_script.py
```

`load_native_module()` searches the normal Python import path and
`cpp/build/{,Release,RelWithDebInfo,Debug}`.

## Core Types

```text
MatrixJob(matrix_id: int, matrix: str, metadata: dict | None = None)
RunConfig(exact=False, full_support=False, include_candidates=False,
          include_timing=True, enable_logging=False, unsafe=False)
MPConfig(workers: int, chunk_size=128, queue_maxsize=4096,
         ordered_results=True, start_method="spawn")
MatrixResult(summary: SummaryRow, candidates: list[CandidateRow], metadata)
```

A matrix may use full CLI form (`"3#4,13/2,..."`) or values only when
`metadata["dimension"]` is present.

`exact=False, unsafe=False` uses candidate-rejector-double: it rejects only
after a one-sided proof and otherwise runs the exact candidate solver.
`unsafe=True` selects the faster heuristic, which can miss exact candidates and
ESS results. `exact=True` bypasses both numerical rejection procedures. Exact
and unsafe together return `EXEC_ERROR`. If candidate-rejector-double is
unavailable on the build or calling thread, default mode also returns
`EXEC_ERROR` before support enumeration and requires explicit exact or unsafe
mode. Matrix input always uses the validating native parser.

`include_timing=False` suppresses the returned timing value; the native call
still measures its execution internally. There is deliberately no computation
timeout.

Status codes:

| Code | Name |
|---:|---|
| 0 | `OK` |
| 1 | `PARSE_ERROR` |
| 2 | `DIMENSION_OUT_OF_RANGE` |
| 3 | `INVALID_VALUE_COUNT` |
| 4 | `EXEC_ERROR` |
| 255 | `INTERNAL_ERROR` |

## Sequential API

```python
from wrapper_v1 import MatrixJob, RunConfig, run_one

job = MatrixJob(3, "3#4,13/2,1/2,5,11/2,3")
result = run_one(job, RunConfig(include_candidates=True), run_id="example")
print(result.summary.ess_count)
```

Native and adapter helpers:

- `load_native_module()` loads `fracessa_core`.
- `compute_job(job, config, run_id) -> MatrixResult` calls the native module.
- `native_status_map()` returns native status constants.
- `result_to_dict(result)` converts dataclass output to plain dictionaries.

Public sequential functions:

- `run_one(job, config=None, run_id=None) -> MatrixResult`
- `run_many(jobs, config=None, run_id=None) -> Iterator[MatrixResult]`
- `run_many_to_sink(jobs, sink, config=None, run_id=None) -> int`

## Multiprocessing Batch API

Use the batch iterator for lists, generators, and large streams. One worker
computes one matrix at a time; parallelism is across matrices.

```python
from wrapper_v1 import (
    MPConfig,
    RunConfig,
    load_jobs_from_verification_json,
    run_jobs_mp,
)

jobs = load_jobs_from_verification_json(
    "python/verification/verification_matrices.json"
)
for result in run_jobs_mp(
    jobs,
    RunConfig(include_candidates=False),
    MPConfig(workers=8, ordered_results=False),
):
    print(result.summary.matrix_id, result.summary.ess_count)
```

Public batch functions:

- `run_jobs_mp(jobs, run_config, mp_config, run_id=None)`
- `run_jobs_mp_to_sink(jobs, sink, run_config, mp_config, run_id=None)`

Submission is bounded to
`min(queue_maxsize, workers * chunk_size)` outstanding jobs. `chunk_size`
controls this window; it does not currently combine jobs into one IPC message.
Set `ordered_results=False` to consume results as workers finish.

## Low-Level Queue API

`MPQueueRunner` exposes `submit()`, `close_input()`, `iter_results()`, and
`shutdown()` for generated input. It is not yet safe for a single thread to
submit an unbounded stream before consuming results: bounded input/output queues
can deadlock, dead workers are not detected while waiting, and shutdown can
terminate legitimate long-running work. See `../reviews/PYTHON_REVIEW.md`.

Until that is fixed, prefer `run_jobs_mp()` for generated iterables. It already
interleaves bounded submission and result consumption.

## Input Loaders

- `load_jobs_from_verification_json(path)` loads the project verification schema.
- `load_jobs_from_json(path, ...)` accepts either a top-level list or
  `{"matrices": [...]}`.

Rows without a `dimension#` matrix prefix must provide a dimension field.

## Output Sinks

- `StreamSink`: one compact JSON object per line.
- `CsvSink`: summary and candidate CSV files.
- `JsonSink`: summary and candidate JSON arrays.
- `ArrowSink`: Arrow IPC files; requires `pyarrow`.
- `ParquetSink`: Parquet files; requires `pyarrow`.
- `MultiSink`: fan-out to several sinks.
- `create_sink(kind, output_dir, run_id)`: constructs a named sink.

For large runs, consume the result iterator continuously or use a sink so
results do not accumulate in memory. Disable candidates when they are not
needed; candidate payloads can dominate output volume.

## Tests

```bash
PYTHONPATH=python python3 -m unittest discover \
  -s python/wrapper_v1/tests -p 'test_*.py'
```

The suite includes unit tests plus native single-process and multiprocessing
integration tests. Build `fracessa_core` first with `./build.sh`.
