# FracESSA Python Wrapper v1 (Pybind)

This wrapper provides a Python API around the native `fracessa_core` pybind module.

It is designed for:
- low overhead single-matrix calls,
- multiprocessing across many matrices,
- queue-driven processing for on-the-fly matrix generation,
- structured outputs to stream/CSV/JSON/Arrow/Parquet.

## 1. Build and Import

From workspace root:

```bash
./fracessa/build.sh
```

Use from Python (workspace root):

```bash
PYTHONPATH=fracessa/python python3 your_script.py
```

In code:

```python
from wrapper_v1 import load_native_module

native = load_native_module()  # loads fracessa_core
```

`load_native_module()` searches:
- current Python import path,
- `fracessa/build`,
- `fracessa/build/Release`,
- `fracessa/build/RelWithDebInfo`,
- `fracessa/build/Debug`.

## 2. Matrix Input Format

`MatrixJob.matrix` accepts:
- full CLI format: `"dimension#values"`, or
- values-only string (e.g. `"4,13/2,..."`) if `metadata["dimension"]` is provided.

Example:

```python
from wrapper_v1 import MatrixJob

job = MatrixJob(matrix_id=3, matrix="3#4,13/2,1/2,5,11/2,3")
```

## 3. Core Types

- `MatrixJob(matrix_id, matrix, metadata=None)`
- `RunConfig(...)`
- `MPConfig(...)`
- `MatrixResult(summary, candidates, metadata)`
- `SummaryRow` and `CandidateRow`

Status codes (`StatusCode`):
- `0` OK
- `1` PARSE_ERROR
- `2` DIMENSION_OUT_OF_RANGE
- `3` INVALID_VALUE_COUNT
- `4` EXEC_ERROR
- `255` INTERNAL_ERROR

## 4. Single-Matrix Usage

```python
from wrapper_v1 import MatrixJob, RunConfig, run_one

job = MatrixJob(matrix_id=3, matrix="3#4,13/2,1/2,5,11/2,3")
cfg = RunConfig(include_candidates=True, exact=False, full_support=False)

result = run_one(job, cfg, run_id="demo_single")
print(result.summary.success, result.summary.ess_count)

for cand in result.candidates:
    print(cand.candidate_id, cand.is_ess, cand.stability)
```

## 5. Sequential Batch Usage

```python
from wrapper_v1 import MatrixJob, RunConfig, run_many

jobs = [
    MatrixJob(matrix_id=1, matrix="2#0,1,0"),
    MatrixJob(matrix_id=2, matrix="2#3,3/2,4"),
]

for result in run_many(jobs, RunConfig(include_candidates=False), run_id="seq_batch"):
    print(result.summary.matrix_id, result.summary.status, result.summary.ess_count)
```

## 6. Multiprocessing: Batch Iterator

Use this when your matrix set is already available as a list/iterator.

```python
from wrapper_v1 import MPConfig, RunConfig, load_jobs_from_verification_json, run_jobs_mp

jobs = load_jobs_from_verification_json("fracessa/python/verification/verification_matrices.json")
mp_cfg = MPConfig(workers=8, ordered_results=True, start_method="spawn")
run_cfg = RunConfig(include_candidates=False)

for result in run_jobs_mp(jobs, run_cfg, mp_cfg, run_id="mp_batch"):
    print(result.summary.matrix_id, result.summary.status, result.summary.elapsed_us)
```

Notes:
- `ordered_results=True` preserves submission order.
- `ordered_results=False` yields results as workers finish.
- `start_method` should usually stay `"spawn"` for portability.
- Batch submission is bounded: at most `min(queue_maxsize, workers * chunk_size)` submitted-but-not-yielded jobs are kept ahead of the consumer. This avoids queue backpressure deadlocks on very large job streams.

## 7. Multiprocessing: Queue Runner (On-the-Fly Producers)

Use this when matrices are generated dynamically.

```python
from wrapper_v1 import MPQueueRunner, MPConfig, RunConfig, MatrixJob

runner = MPQueueRunner(
    run_config=RunConfig(include_candidates=False),
    mp_config=MPConfig(workers=8, ordered_results=False),
    run_id="mp_queue",
)

try:
    for i in range(100):
        # Produce/generate matrix strings on the fly.
        runner.submit(MatrixJob(matrix_id=i, matrix="2#0,1,0"))

    runner.close_input()

    for result in runner.iter_results(expected_results=runner.submitted):
        print(result.summary.matrix_id, result.summary.status)
finally:
    runner.shutdown()
```

## 8. Output Sinks

Available sinks:
- `StreamSink` (JSON-lines to stdout or custom stream)
- `CsvSink` (summary CSV + candidates CSV)
- `JsonSink` (summary JSON array + candidates JSON array)
- `ArrowSink` (Arrow IPC stream files; requires `pyarrow`)
- `ParquetSink` (Parquet files; requires `pyarrow`)
- `MultiSink` (fan-out to multiple sinks)
- `create_sink(kind, output_dir, run_id)` helper

### 8.1 Quick sink creation

```python
from wrapper_v1 import create_sink

sink = create_sink("csv", "./out", run_id="bench01")
# writes:
# ./out/bench01_summary.csv
# ./out/bench01_candidates.csv
```

### 8.2 Write sequential results to sink

```python
from wrapper_v1 import MatrixJob, RunConfig, run_many_to_sink, create_sink

jobs = [MatrixJob(matrix_id=1, matrix="2#0,1,0")]
sink = create_sink("json", "./out", run_id="seq01")
run_many_to_sink(jobs, sink, RunConfig(include_candidates=True), run_id="seq01")
```

### 8.3 Write multiprocessing results to sink

```python
from wrapper_v1 import load_jobs_from_verification_json, run_jobs_mp_to_sink, RunConfig, MPConfig, create_sink

jobs = load_jobs_from_verification_json("fracessa/python/verification/verification_matrices.json")
sink = create_sink("parquet", "./out", run_id="mp01")
run_jobs_mp_to_sink(jobs, sink, RunConfig(include_candidates=True), MPConfig(workers=8), run_id="mp01")
```

## 9. JSON Job Loaders

### Verification schema loader

```python
from wrapper_v1 import load_jobs_from_verification_json

jobs = load_jobs_from_verification_json("fracessa/python/verification/verification_matrices.json")
```

### Generic JSON loader

```python
from wrapper_v1 import load_jobs_from_json

jobs = load_jobs_from_json("my_jobs.json")
```

Supported shapes:
- `{"matrices": [ ... ]}`
- `[ ... ]`

If a row `matrix` value has no `"dimension#"` prefix, the row must include a `dimension` field.

## 10. Performance and Behavior Notes

- Wrapper calls are in-process via pybind (no subprocess launch overhead).
- Parallelism is process-based in Python (`multiprocessing`), one matrix job per worker call.
- Candidate payload can be large; disable with `include_candidates=False` if you only need ESS counts/timings.
- For millions of matrices, prefer `run_jobs_mp_to_sink(...)` or consume `run_jobs_mp(...)` continuously so result payloads are not accumulated in Python memory.

## 11. Common Errors

- `ModuleNotFoundError: fracessa_core`
  - Run `./fracessa/build.sh` first.
  - Ensure `PYTHONPATH=fracessa/python` when launching scripts.

- Multiprocessing semaphore permission errors in restricted environments
  - Run outside restricted sandboxes or use sequential APIs.

## 12. API Reference Files

- High-level API contract: `fracessa/python/wrapper_v1/API_CONTRACT.md`
- Public exports: `fracessa/python/wrapper_v1/__init__.py`
- Native adapter: `fracessa/python/wrapper_v1/core.py`

## 13. Wrapper Test Suite

Run wrapper tests directly:

```bash
PYTHONPATH=fracessa/python python3 -m unittest discover -s fracessa/python/wrapper_v1/tests -p "test_*.py"
```

These tests include:
- pure unit tests (I/O, sink formats, API mapping),
- native single-process integration checks,
- native multiprocessing checks (skipped automatically if process primitives are blocked).
