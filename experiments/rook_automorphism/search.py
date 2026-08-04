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
import re
import sys
import tempfile
from time import monotonic


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "python"))
sys.path.insert(0, str(REPOSITORY_ROOT))

from fracessa import MPConfig, Matrix, RunConfig, run_multiprocessing  # noqa: E402


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

DEFAULT_CPUS = tuple(range(3, 9))
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
        "last_received": None,
        "completed_through": None,
    }


def _migrate_variables(item: dict | None) -> None:
    if item is not None and "variables" not in item and all(name in item for name in ("A", "B", "C")):
        variables = [item["A"], item["B"], item["C"]]
        migrated = {}
        for key, value in item.items():
            if key == "A":
                migrated["variables"] = variables
            elif key not in ("B", "C"):
                migrated[key] = value
        item.clear()
        item.update(migrated)


def _ensure_state(state: dict) -> None:
    for template in TEMPLATES:
        group = state.setdefault(f"{template.dimension}:{template.group}", {})
        work = group.setdefault("work", _new_work())
        for obsolete in ("active_stage", "received_in_stage", "expected_in_stage", "best_unverified_fast"):
            work.pop(obsolete, None)
        for key, value in _new_work().items():
            work.setdefault(key, value)
        group.setdefault("result", None)
        for item in (work["last_received"], work["completed_through"], group["result"]):
            _migrate_variables(item)
        result = group["result"]
        if result is not None and "safe_elapsed_ns" in result:
            result.pop("safe_elapsed_ns")
            result["fast_elapsed_ns"] = None
        if result is not None:
            result.pop("safe_fallback", None)


def _load_state(path: Path) -> dict:
    if not path.exists():
        state = {}
    else:
        state = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(state, dict):
            raise ValueError(f"{path}: state root must be a JSON object")
    for dimension, groups in tuple(state.items()):
        if ":" in dimension:
            continue
        for group, value in groups.items():
            state[f"{dimension}:{group}"] = value
        del state[dimension]
    _ensure_state(state)
    return state


def _write_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp")
    text = json.dumps(state, indent=2, ensure_ascii=False)
    text = re.sub(
        r'"variables": \[\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\s*\]',
        lambda match: f'"variables": [{match.group(1)}, {match.group(2)}, {match.group(3)}]',
        text,
    )

    def compact_structure(match: re.Match) -> str:
        entries = re.findall(r'"([^"]+)": (-?\d+)', match.group(2))
        contents = ", ".join(f'"{key}": {value}' for key, value in entries)
        return f'"{match.group(1)}": {{{contents}}}'

    text = re.sub(
        r'"(candidate_structure|ess_structure)": \{\n((?:\s+"[^"]+": -?\d+,?\n)*)\s*\}',
        compact_structure,
        text,
    )
    temporary.write_text(text + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _group_state(state: dict, template: RookTemplate) -> dict:
    return state[f"{template.dimension}:{template.group}"]


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
        "variables": [a, b, c],
    }


def _search_jobs(through_value: int) -> list[tuple[int, int, tuple[int, int, int]]]:
    jobs = []
    for stage_index in range(_stage_index(through_value) + 1):
        stage_value = _stage_value(stage_index)
        jobs.extend((stage_index, stage_value, triple) for triple in _stage_triples(stage_index))
    return jobs


def _search_job_stream(through_value: int):
    for stage_index in range(_stage_index(through_value) + 1):
        stage_value = _stage_value(stage_index)
        yield from ((stage_index, stage_value, triple) for triple in _stage_triples(stage_index))


def _record_completion(
    work: dict,
    returned_above_frontier: set[int],
    job_index: int,
    jobs: list[tuple[int, int, tuple[int, int, int]]],
) -> None:
    if job_index <= _frontier(work) or job_index in returned_above_frontier:
        raise RuntimeError(f"duplicate completed job index {job_index}")

    stage_index, stage_value, triple = jobs[job_index]
    work["last_received"] = _progress_point(job_index, stage_index, stage_value, triple)
    returned_above_frontier.add(job_index)

    next_index = _frontier(work) + 1
    while next_index in returned_above_frontier:
        returned_above_frontier.remove(next_index)
        next_stage_index, next_stage_value, next_triple = jobs[next_index]
        work["completed_through"] = _progress_point(next_index, next_stage_index, next_stage_value, next_triple)
        next_index += 1


def _record_count(group: dict) -> int:
    result = group["result"]
    return -1 if result is None else int(result["ess_count"])


def _check_native_result(result: dict) -> None:
    if int(result["status"]) != 0:
        metadata = result.get("metadata") or {}
        label = f"{metadata.get('dimension', '?')}:{metadata.get('group', '?')}"
        raise RuntimeError(f"{label}: {result['error_message'] or 'native computation failed'}")


def _consider_fast_result(template: RookTemplate, group: dict, fast: dict) -> None:
    ess_count = int(fast["ess_count"])
    if ess_count <= _record_count(group):
        return

    metadata = fast["metadata"]
    triple = [int(metadata["A"]), int(metadata["B"]), int(metadata["C"])]
    group["result"] = {
        "matrix": _matrix_string(template, tuple(triple)),
        "variables": triple,
        "candidate_count": int(fast["candidate_count"]),
        "ess_count": ess_count,
        "candidate_structure": {str(size): int(count) for size, count in fast["candidate_structure"].items()},
        "ess_structure": {str(size): int(count) for size, count in fast["ess_structure"].items()},
        "fast_elapsed_ns": int(fast["elapsed_ns"]),
    }
    print(
        f"new maximum {template.dimension}:{template.group}: "
        f"ESS={ess_count}, candidates={fast['candidate_count']}, (A,B,C)={tuple(triple)}",
        flush=True,
    )


def _run_pipeline(
    templates: tuple[RookTemplate, ...],
    state: dict,
    through_value: int,
    checkpoint,
    mp_config: MPConfig,
) -> None:
    jobs = []
    job_source = _search_job_stream(through_value)

    def ensure_job(job_index: int) -> bool:
        while len(jobs) <= job_index:
            try:
                jobs.append(next(job_source))
            except StopIteration:
                return False
        return True

    assert ensure_job(0)
    pending = {template: set() for template in templates}

    # The all-zero game is known analytically and needs no native run.
    for template in templates:
        work = _group_state(state, template)["work"]
        if _frontier(work) < 0:
            _record_completion(work, pending[template], 0, jobs)
    checkpoint(force=True)

    first_job_index = min(_frontier(_group_state(state, template)["work"]) + 1 for template in templates)
    if not ensure_job(first_job_index):
        print(f"search through {jobs[-1][1]:+d} already complete", flush=True)
        return

    print(
        f"pipelined search started: first outstanding job={first_job_index}, target={through_value:+d}",
        flush=True,
    )

    def matrices():
        job_index = first_job_index
        while ensure_job(job_index):
            stage_index, stage_value, triple = jobs[job_index]
            for template in templates:
                frontier = _frontier(_group_state(state, template)["work"])
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
            job_index += 1

    templates_by_key = {(template.dimension, template.group): template for template in templates}
    results = run_multiprocessing("fast", matrices(), config=RunConfig(include_candidates=False), mp_config=mp_config)
    returned = 0
    with closing(results):
        for fast in results:
            _check_native_result(fast)
            metadata = fast["metadata"]
            template = templates_by_key[(int(metadata["dimension"]), str(metadata["group"]))]
            work = _group_state(state, template)["work"]
            _record_completion(work, pending[template], int(metadata["job_index"]), jobs)
            _consider_fast_result(template, _group_state(state, template), fast)
            returned += 1
            if checkpoint():
                print(f"search progress: {returned} matrix runs returned in this invocation", flush=True)

    for template in templates:
        if _frontier(_group_state(state, template)["work"]) != len(jobs) - 1:
            raise RuntimeError(f"{template.dimension}:{template.group}: search did not complete")
    checkpoint(force=True)
    _print_stage_report(templates, state, jobs[-1][1])


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
            f"(A,B,C)={tuple(result['variables'])}",
            flush=True,
        )


def run_search(
    through_value: int,
    max_dimension: int,
    dimensions: tuple[int, ...] | None,
    checkpoint_seconds: float,
    state_path: Path,
    cpu_ids: tuple[int, ...],
    mp_config: MPConfig,
) -> None:
    templates = tuple(
        template
        for template in TEMPLATES
        if (template.dimension in dimensions if dimensions is not None else template.dimension <= max_dimension)
    )
    if not templates:
        raise ValueError("dimension selection excludes every compact rook template")

    selection = f"dimensions={','.join(map(str, dimensions))}" if dimensions is not None else f"max_dimension={max_dimension}"
    print(
        f"rook search starting: CPUs={','.join(map(str, cpu_ids))}, workers={mp_config.workers}, "
        f"target={through_value:+d}, {selection}, state={state_path}",
        flush=True,
    )

    state = _load_state(state_path)
    last_checkpoint = monotonic()

    def checkpoint(*, force: bool = False) -> bool:
        nonlocal last_checkpoint
        now = monotonic()
        if force or now - last_checkpoint >= checkpoint_seconds:
            _write_state(state_path, state)
            last_checkpoint = now
            return True
        return False

    checkpoint(force=True)
    try:
        _run_pipeline(templates, state, through_value, checkpoint, mp_config)
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
    assert _matrix_string(TEMPLATES[0], (1, 2, 3)) == "10#3,1,3,1,2"
    jobs = _search_jobs(-1)
    assert [job[1] for job in jobs] == [0] + [1] * 5 + [-1] * 12

    state = {}
    _ensure_state(state)
    first_key = f"{TEMPLATES[0].dimension}:{TEMPLATES[0].group}"
    assert set(state[first_key]) == {"work", "result"}
    work = state[first_key]["work"]
    _record_completion(work, set(), 0, jobs)
    returned = set()
    _record_completion(work, returned, 2, jobs)
    assert _frontier(work) == 0
    _record_completion(work, returned, 1, jobs)
    assert _frontier(work) == 2

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "state.json"
        _write_state(path, state)
        assert _load_state(path)[first_key]["work"]["completed_through"]["job_index"] == 2


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, default=1_000_000, help="final value stage; default: 1000000")
    dimension_selection = parser.add_mutually_exclusive_group()
    dimension_selection.add_argument("--max-dimension", type=int, default=30)
    dimension_selection.add_argument("--dimension", type=int, nargs="+", metavar="N", help="run only these dimensions")
    parser.add_argument("--checkpoint-seconds", type=float, default=10.0)
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--cpu", type=int, nargs="+", default=DEFAULT_CPUS, metavar="ID", help="CPU IDs to use; default: 3 4 5 6 7 8")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.max_dimension < 1 or args.max_dimension > 30:
        parser.error("--max-dimension must be between 1 and 30")
    if args.dimension is not None:
        if len(set(args.dimension)) != len(args.dimension):
            parser.error("--dimension values must be unique")
        available_dimensions = {template.dimension for template in TEMPLATES}
        unavailable_dimensions = set(args.dimension) - available_dimensions
        if unavailable_dimensions:
            parser.error(f"dimensions without a compact rook template: {sorted(unavailable_dimensions)}")
    if args.checkpoint_seconds <= 0:
        parser.error("--checkpoint-seconds must be positive")
    if not args.self_check:
        if not hasattr(os, "sched_getaffinity") or not hasattr(os, "sched_setaffinity"):
            parser.error("--cpu requires Linux CPU-affinity support")
        if len(set(args.cpu)) != len(args.cpu):
            parser.error("--cpu IDs must be unique")
        unavailable = set(args.cpu) - set(os.sched_getaffinity(0))
        if unavailable:
            parser.error(f"unavailable --cpu IDs: {sorted(unavailable)}")
    return args


def main() -> int:
    args = _parse_args()
    if args.self_check:
        _self_check()
        print("rook search self-check passed")
        return 0
    cpu_ids = tuple(args.cpu)
    previous_affinity = set(os.sched_getaffinity(0))
    os.sched_setaffinity(0, set(cpu_ids))
    try:
        run_search(
            args.through,
            args.max_dimension,
            None if args.dimension is None else tuple(args.dimension),
            args.checkpoint_seconds,
            args.state,
            cpu_ids,
            MPConfig(workers=len(cpu_ids)),
        )
    except KeyboardInterrupt:
        print("interrupted; state saved", file=sys.stderr)
        return 130
    finally:
        os.sched_setaffinity(0, previous_affinity)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
