from __future__ import annotations

from collections.abc import Iterable, Iterator

from .core import compute_job, new_run_id
from .types import MatrixJob, MatrixResult, RunConfig


def run_one(job: MatrixJob, config: RunConfig | None = None, run_id: str | None = None) -> MatrixResult:
    cfg = config if config is not None else RunConfig()
    rid = run_id or new_run_id("single")
    return compute_job(job=job, config=cfg, run_id=rid)


def run_many(
    jobs: Iterable[MatrixJob],
    config: RunConfig | None = None,
    run_id: str | None = None,
) -> Iterator[MatrixResult]:
    cfg = config if config is not None else RunConfig()
    rid = run_id or new_run_id("single")
    for job in jobs:
        yield compute_job(job=job, config=cfg, run_id=rid)


def run_many_to_sink(
    jobs: Iterable[MatrixJob],
    sink,
    config: RunConfig | None = None,
    run_id: str | None = None,
) -> int:
    count = 0
    try:
        for result in run_many(jobs=jobs, config=config, run_id=run_id):
            sink.write_result(result)
            count += 1
    finally:
        sink.close()
    return count
