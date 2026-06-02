# Wrapper v1 API Contract (Pybind-based)

This wrapper uses the native `fracessa_core` pybind module and does not use `python/fracessa_py.py`.

## Core types
- `MatrixJob(matrix_id: int, matrix: str, metadata: dict | None = None)`
- `RunConfig(exact=False, full_support=False, include_candidates=False, include_timing=True, unsafe=False, enable_logging=False)`
- `MPConfig(workers: int, chunk_size=128, queue_maxsize=4096, ordered_results=True, start_method='spawn'|'forkserver'|'fork')`
- `MatrixResult(summary: SummaryRow, candidates: list[CandidateRow], metadata: dict | None)`

## Single-process API
- `run_one(job, config=None, run_id=None) -> MatrixResult`
- `run_many(jobs, config=None, run_id=None) -> Iterator[MatrixResult]`
- `run_many_to_sink(jobs, sink, config=None, run_id=None) -> int`

## Multiprocessing API

### Batch style
- `run_jobs_mp(jobs, run_config, mp_config, run_id=None) -> Iterator[MatrixResult]`
- `run_jobs_mp_to_sink(jobs, sink, run_config, mp_config, run_id=None) -> int`

Batch execution streams jobs through a bounded in-flight window:
`max_buffered = min(mp_config.queue_maxsize, mp_config.workers * mp_config.chunk_size)`.
This prevents bounded queue deadlocks when processing large or generated job streams.

### Queue style (for on-the-fly matrix generation)
`MPQueueRunner`:
1. `runner = MPQueueRunner(run_config, mp_config, run_id=None)`
2. `runner.submit(job)` repeatedly
3. `runner.close_input()`
4. iterate `runner.iter_results(expected_results=runner.submitted)`
5. `runner.shutdown()`

## Input loading helpers
- `load_jobs_from_verification_json(path) -> list[MatrixJob]`
- `load_jobs_from_json(path, ...) -> list[MatrixJob]`

If `MatrixJob.matrix` is not in `dimension#values` format, `metadata['dimension']` must be present.

## Output backends (sinks)
- `CsvSink(summary_path, candidates_path)`
- `JsonSink(summary_path, candidates_path)`
- `ArrowSink(summary_path, candidates_path)` (requires `pyarrow`)
- `ParquetSink(summary_path, candidates_path)` (requires `pyarrow`)
- `StreamSink(stream=None, flush_each=False)`
- `create_sink(kind, output_dir, run_id, stream=None)` with kind in `{csv,json,arrow,parquet,stream}`

## Native module entrypoint
- `load_native_module()` loads `fracessa_core` from default Python path or repo build paths.
- `compute_job(job, config, run_id) -> MatrixResult` runs one matrix via `fracessa_core.compute_matrix(...)`.

## Status codes
- `0 OK`
- `1 PARSE_ERROR`
- `2 DIMENSION_OUT_OF_RANGE`
- `3 INVALID_VALUE_COUNT`
- `4 EXEC_ERROR`
- `255 INTERNAL_ERROR`
