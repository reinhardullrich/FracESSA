from __future__ import annotations

from collections.abc import Iterable, Iterator

from .core import compute_matrix, new_run_id
from .sinks import _consume_to_sink
from .types import Matrix, RunConfig


def _run_matrices(
    matrices: Iterable[Matrix],
    config: RunConfig,
    run_id: str,
) -> Iterator[dict]:
    for matrix in matrices:
        yield compute_matrix(matrix=matrix, config=config, run_id=run_id)


def run(
    matrices: Matrix | Iterable[Matrix],
    config: RunConfig | None = None,
    run_id: str | None = None,
    sink=None,
) -> dict | Iterator[dict] | int:
    cfg = config if config is not None else RunConfig()
    rid = run_id or new_run_id("run")

    if isinstance(matrices, Matrix):
        if sink is None:
            return compute_matrix(matrix=matrices, config=cfg, run_id=rid)
        results = _run_matrices((matrices,), cfg, rid)
    else:
        results = _run_matrices(matrices, cfg, rid)

    if sink is not None:
        return _consume_to_sink(results, sink)
    return results
