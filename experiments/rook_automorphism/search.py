#!/usr/bin/env python3
"""Search compact rook matrices over expanding primitive integer triples."""

from __future__ import annotations

import argparse
from contextlib import closing
from dataclasses import dataclass
from itertools import product
import json
from math import gcd
import os
from pathlib import Path
import sys
import tempfile
from time import monotonic


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "python"))
sys.path.insert(0, str(REPOSITORY_ROOT))

from fracessa import Matrix, RunConfig, run, run_multiprocessing  # noqa: E402


@dataclass(frozen=True, slots=True)
class RookTemplate:
    dimension: int
    rows: int
    columns: int
    rank: int
    pattern: tuple[str, ...]

    @property
    def group(self) -> str:
        return f"{self.rows}x{self.columns}"


TEMPLATES = (
    RookTemplate(6, 2, 3, 1, ("C", "A", "B")),
    RookTemplate(10, 2, 5, 1, ("C", "A", "C", "A", "B")),
    RookTemplate(12, 3, 4, 2, ("C", "C", "A", "B", "C", "A")),
    RookTemplate(14, 2, 7, 1, ("C", "A", "C", "A", "C", "A", "B")),
    RookTemplate(15, 3, 5, 1, ("C", "C", "A", "C", "B", "A", "C")),
    RookTemplate(18, 2, 9, 1, ("C", "A", "C", "A", "C", "A", "C", "A", "B")),
    RookTemplate(20, 4, 5, 2, ("C", "C", "C", "A", "B", "C", "C", "A", "C", "B")),
    RookTemplate(21, 3, 7, 1, ("C", "C", "A", "C", "C", "A", "B", "C", "A", "C")),
    RookTemplate(22, 2, 11, 1, ("C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "B")),
    RookTemplate(24, 3, 8, 2, ("C", "C", "A", "C", "C", "A", "C", "B", "A", "C", "C", "A")),
    RookTemplate(26, 2, 13, 1, ("C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "B")),
    RookTemplate(28, 4, 7, 2, ("C", "C", "C", "A", "C", "C", "B", "A", "C", "C", "C", "A", "C", "B")),
    RookTemplate(30, 2, 15, 1, ("C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "C", "A", "B")),
    RookTemplate(30, 3, 10, 2, ("C", "C", "A", "C", "C", "A", "C", "C", "A", "B", "C", "A", "C", "C", "A")),
    RookTemplate(30, 5, 6, 3, ("C", "C", "C", "C", "A", "B", "C", "C", "C", "A", "C", "B", "C", "C", "A")),
)

SEED = (7, 7, 15)
DEFAULT_STATE = Path(__file__).with_name("state.json")


def _stage_value(stage_index: int) -> int:
    if stage_index == 0:
        return 0
    magnitude = (stage_index + 1) // 2
    return magnitude if stage_index % 2 else -magnitude


def _stage_index(value: int) -> int:
    if value == 0:
        return 0
    return 2 * abs(value) - (1 if value > 0 else 0)


def _stage_triples(stage_index: int) -> list[tuple[int, int, int]]:
    values = tuple(_stage_value(index) for index in range(stage_index + 1))
    position = {value: index for index, value in enumerate(values)}
    newest = values[-1]
    triples = []
    for triple in product(values, repeat=3):
        a, b, c = triple
        if newest not in triple or position[a] > position[b]:
            continue
        if triple != (0, 0, 0) and gcd(gcd(abs(a), abs(b)), abs(c)) != 1:
            continue
        triples.append(triple)
    return triples


def _matrix_string(template: RookTemplate, triple: tuple[int, int, int]) -> str:
    a, b, c = triple
    values = {"A": a, "B": b, "C": c}
    return f"{template.dimension}#{','.join(str(values[token]) for token in template.pattern)}"


def _new_work() -> dict:
    return {
        "active_stage": None,
        "last_received": None,
        "completed_through": None,
        "received_in_stage": 0,
        "expected_in_stage": 0,
        "best_unverified_fast": None,
    }


def _ensure_state(state: dict) -> None:
    for template in TEMPLATES:
        dimension = state.setdefault(str(template.dimension), {})
        group = dimension.setdefault(template.group, {})
        work = group.setdefault("work", _new_work())
        for key, value in _new_work().items():
            work.setdefault(key, value)
        group.setdefault("result", None)


def _load_state(path: Path) -> dict:
    if not path.exists():
        state = {}
    else:
        state = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(state, dict):
            raise ValueError(f"{path}: state root must be a JSON object")
    _ensure_state(state)
    return state


def _write_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp")
    temporary.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _group_state(state: dict, template: RookTemplate) -> dict:
    return state[str(template.dimension)][template.group]


def _frontier(work: dict) -> int:
    completed = work["completed_through"]
    return -1 if completed is None else int(completed["job_index"])


def _progress_point(
    job_index: int,
    stage_index: int,
    stage_value: int,
    triple: tuple[int, int, int],
) -> dict:
    a, b, c = triple
    return {
        "job_index": job_index,
        "value_stage_index": stage_index,
        "value_stage": stage_value,
        "A": a,
        "B": b,
        "C": c,
    }


def _prepare_stage(work: dict, stage_index: int, stage_value: int, stage_start: int, expected: int) -> None:
    completed = min(expected, max(0, _frontier(work) - stage_start + 1))
    work["active_stage"] = {"index": stage_index, "value": stage_value}
    work["received_in_stage"] = completed
    work["expected_in_stage"] = expected


def _record_completion(
    work: dict,
    returned_above_frontier: set[int],
    job_index: int,
    stage_index: int,
    stage_value: int,
    stage_start: int,
    triples: list[tuple[int, int, int]],
) -> None:
    if job_index <= _frontier(work) or job_index in returned_above_frontier:
        raise RuntimeError(f"duplicate completed job index {job_index}")

    triple = triples[job_index - stage_start]
    work["last_received"] = _progress_point(job_index, stage_index, stage_value, triple)
    work["received_in_stage"] += 1
    returned_above_frontier.add(job_index)

    next_index = _frontier(work) + 1
    while next_index in returned_above_frontier:
        returned_above_frontier.remove(next_index)
        next_triple = triples[next_index - stage_start]
        work["completed_through"] = _progress_point(next_index, stage_index, stage_value, next_triple)
        next_index += 1


def _verified_count(group: dict) -> int:
    result = group["result"]
    return -1 if result is None else int(result["ess_count"])


def _check_native_result(result: dict) -> None:
    if int(result["status"]) != 0:
        metadata = result.get("metadata") or {}
        label = f"{metadata.get('dimension', '?')}:{metadata.get('group', '?')}"
        raise RuntimeError(f"{label}: {result['error_message'] or 'native computation failed'}")


def _pending_fast_result(result: dict) -> dict:
    metadata = result["metadata"]
    return {
        "job_index": int(result["matrix_id"]),
        "A": int(metadata["A"]),
        "B": int(metadata["B"]),
        "C": int(metadata["C"]),
        "fast_ess_count": int(result["ess_count"]),
    }


def _verify_pending(template: RookTemplate, group: dict, pending: dict) -> bool:
    triple = (int(pending["A"]), int(pending["B"]), int(pending["C"]))
    matrix = _matrix_string(template, triple)
    matrix_id = int(pending["job_index"])
    safe = run("safe", Matrix(matrix_id=matrix_id, matrix=matrix), RunConfig(include_candidates=False))
    _check_native_result(safe)
    candidate_count = int(safe["candidate_count"])
    ess_count = int(safe["ess_count"])
    candidate_structure = {str(size): int(count) for size, count in safe["candidate_structure"].items()}
    ess_structure = {str(size): int(count) for size, count in safe["ess_structure"].items()}

    changed = ess_count > _verified_count(group)
    if changed:
        a, b, c = triple
        group["result"] = {
            "matrix": matrix,
            "A": a,
            "B": b,
            "C": c,
            "candidate_count": candidate_count,
            "ess_count": ess_count,
            "candidate_structure": candidate_structure,
            "ess_structure": ess_structure,
            "safe_elapsed_ns": int(safe["elapsed_ns"]),
            "safe_fallback": safe["safe_fallback"],
        }
        print(
            f"new verified maximum {template.dimension}:{template.group}: "
            f"ESS={ess_count}, candidates={candidate_count}, (A,B,C)={triple}",
            flush=True,
        )
    return changed


def _consider_fast_result(template: RookTemplate, group: dict, fast: dict, checkpoint) -> None:
    fast_count = int(fast["ess_count"])
    if fast_count <= _verified_count(group):
        return

    pending = _pending_fast_result(fast)
    group["work"]["best_unverified_fast"] = pending
    checkpoint(force=True)
    _verify_pending(template, group, pending)
    group["work"]["best_unverified_fast"] = None
    checkpoint(force=True)


def _resume_pending_verifications(templates: tuple[RookTemplate, ...], state: dict, checkpoint) -> None:
    for template in templates:
        group = _group_state(state, template)
        pending = group["work"]["best_unverified_fast"]
        if pending is None:
            continue
        if int(pending["fast_ess_count"]) > _verified_count(group):
            _verify_pending(template, group, pending)
        group["work"]["best_unverified_fast"] = None
        checkpoint(force=True)


def _seed_missing_results(templates: tuple[RookTemplate, ...], state: dict, checkpoint) -> None:
    missing = tuple(template for template in templates if _group_state(state, template)["result"] is None)
    if not missing:
        return

    def matrices():
        for template_index, template in enumerate(missing, start=1):
            a, b, c = SEED
            yield Matrix(
                matrix_id=-template_index,
                matrix=_matrix_string(template, SEED),
                metadata={
                    "dimension": template.dimension,
                    "group": template.group,
                    "A": a,
                    "B": b,
                    "C": c,
                },
            )

    templates_by_key = {(template.dimension, template.group): template for template in missing}
    results = run_multiprocessing("fast", matrices(), config=RunConfig(include_candidates=False))
    with closing(results):
        for fast in results:
            _check_native_result(fast)
            metadata = fast["metadata"]
            template = templates_by_key[(int(metadata["dimension"]), str(metadata["group"]))]
            _consider_fast_result(template, _group_state(state, template), fast, checkpoint)


def _run_stage(
    templates: tuple[RookTemplate, ...],
    state: dict,
    stage_index: int,
    stage_start: int,
    triples: list[tuple[int, int, int]],
    checkpoint,
) -> None:
    stage_value = _stage_value(stage_index)
    stage_end = stage_start + len(triples) - 1
    pending = {template: set() for template in templates}

    for template in templates:
        work = _group_state(state, template)["work"]
        if _frontier(work) < stage_start - 1:
            raise RuntimeError(f"{template.dimension}:{template.group}: earlier value stage is incomplete")
        _prepare_stage(work, stage_index, stage_value, stage_start, len(triples))

    if triples == [(0, 0, 0)]:
        for template in templates:
            work = _group_state(state, template)["work"]
            if _frontier(work) < stage_start:
                _record_completion(work, pending[template], stage_start, stage_index, stage_value, stage_start, triples)
        checkpoint(force=True)
        return

    def matrices():
        for template in templates:
            frontier = _frontier(_group_state(state, template)["work"])
            for offset, triple in enumerate(triples):
                job_index = stage_start + offset
                if job_index <= frontier:
                    continue
                a, b, c = triple
                yield Matrix(
                    matrix_id=job_index,
                    matrix=_matrix_string(template, triple),
                    metadata={
                        "dimension": template.dimension,
                        "group": template.group,
                        "rank": template.rank,
                        "value_stage_index": stage_index,
                        "value_stage": stage_value,
                        "job_index": job_index,
                        "A": a,
                        "B": b,
                        "C": c,
                    },
                )

    templates_by_key = {(template.dimension, template.group): template for template in templates}
    results = run_multiprocessing("fast", matrices(), config=RunConfig(include_candidates=False))
    with closing(results):
        for fast in results:
            _check_native_result(fast)
            metadata = fast["metadata"]
            template = templates_by_key[(int(metadata["dimension"]), str(metadata["group"]))]
            work = _group_state(state, template)["work"]
            _record_completion(
                work,
                pending[template],
                int(metadata["job_index"]),
                stage_index,
                stage_value,
                stage_start,
                triples,
            )
            checkpoint()
            _consider_fast_result(template, _group_state(state, template), fast, checkpoint)

    for template in templates:
        if _frontier(_group_state(state, template)["work"]) != stage_end:
            raise RuntimeError(f"{template.dimension}:{template.group}: value stage {stage_value:+d} did not complete")
    checkpoint(force=True)


def _print_stage_report(templates: tuple[RookTemplate, ...], state: dict, stage_value: int) -> None:
    print(f"value stage {stage_value:+d} complete", flush=True)
    for template in templates:
        result = _group_state(state, template)["result"]
        if result is None:
            continue
        ess_count = int(result["ess_count"])
        gamma = ess_count ** (1.0 / template.dimension) if ess_count else 0.0
        print(
            f"  {template.dimension}:{template.group} ESS={ess_count} gamma={gamma:.9f} "
            f"(A,B,C)=({result['A']},{result['B']},{result['C']})",
            flush=True,
        )


def run_search(through_value: int, max_dimension: int, checkpoint_seconds: float, state_path: Path) -> None:
    templates = tuple(template for template in TEMPLATES if template.dimension <= max_dimension)
    if not templates:
        raise ValueError("max dimension excludes every compact rook template")

    state = _load_state(state_path)
    last_checkpoint = monotonic()

    def checkpoint(*, force: bool = False) -> None:
        nonlocal last_checkpoint
        now = monotonic()
        if force or now - last_checkpoint >= checkpoint_seconds:
            _write_state(state_path, state)
            last_checkpoint = now

    checkpoint(force=True)
    try:
        _resume_pending_verifications(templates, state, checkpoint)
        _seed_missing_results(templates, state, checkpoint)

        next_job_index = 0
        for stage_index in range(_stage_index(through_value) + 1):
            triples = _stage_triples(stage_index)
            stage_start = next_job_index
            stage_end = stage_start + len(triples) - 1
            next_job_index = stage_end + 1
            if all(_frontier(_group_state(state, template)["work"]) >= stage_end for template in templates):
                continue
            _run_stage(templates, state, stage_index, stage_start, triples, checkpoint)
            _print_stage_report(templates, state, _stage_value(stage_index))
    finally:
        checkpoint(force=True)


def _self_check() -> None:
    assert len({(template.dimension, template.group) for template in TEMPLATES}) == len(TEMPLATES)
    assert all(len(template.pattern) == template.dimension // 2 for template in TEMPLATES)
    assert [_stage_value(index) for index in range(5)] == [0, 1, -1, 2, -2]
    assert [_stage_index(value) for value in (0, 1, -1, 2, -2)] == list(range(5))
    assert _stage_triples(0) == [(0, 0, 0)]
    for stage_index in range(1, 5):
        newest = _stage_value(stage_index)
        triples = _stage_triples(stage_index)
        positions = {value: index for index, value in enumerate(_stage_value(i) for i in range(stage_index + 1))}
        assert len(triples) == len(set(triples))
        assert all(newest in triple for triple in triples)
        assert all(gcd(gcd(abs(a), abs(b)), abs(c)) == 1 for a, b, c in triples)
        assert all(positions[a] <= positions[b] for a, b, _ in triples)
    assert _matrix_string(TEMPLATES[0], (1, 2, 3)) == "6#3,1,2"

    state = {}
    _ensure_state(state)
    assert set(state["6"]["2x3"]) == {"work", "result"}
    work = state["6"]["2x3"]["work"]
    stage_zero = _stage_triples(0)
    _prepare_stage(work, 0, 0, 0, 1)
    _record_completion(work, set(), 0, 0, 0, 0, stage_zero)
    stage_one = _stage_triples(1)
    _prepare_stage(work, 1, 1, 1, len(stage_one))
    returned = set()
    _record_completion(work, returned, 2, 1, 1, 1, stage_one)
    assert _frontier(work) == 0
    _record_completion(work, returned, 1, 1, 1, 1, stage_one)
    assert _frontier(work) == 2

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "state.json"
        _write_state(path, state)
        assert _load_state(path)["6"]["2x3"]["work"]["completed_through"]["job_index"] == 2


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, help="final value stage, for example -2")
    parser.add_argument("--max-dimension", type=int, default=30)
    parser.add_argument("--checkpoint-seconds", type=float, default=10.0)
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if not args.self_check and args.through is None:
        parser.error("--through is required unless --self-check is used")
    if args.max_dimension < 1 or args.max_dimension > 30:
        parser.error("--max-dimension must be between 1 and 30")
    if args.checkpoint_seconds <= 0:
        parser.error("--checkpoint-seconds must be positive")
    return args


def main() -> int:
    args = _parse_args()
    if args.self_check:
        _self_check()
        print("rook search self-check passed")
        return 0
    try:
        run_search(args.through, args.max_dimension, args.checkpoint_seconds, args.state)
    except KeyboardInterrupt:
        print("interrupted; state saved", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
