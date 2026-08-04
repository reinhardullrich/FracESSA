#!/usr/bin/env python3
"""Search controlled S3 x D8 perturbations of the 24-strategy rook record."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations, product
import json
from math import gcd
import os
from pathlib import Path
import random
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = ROOT / "cpp" / "build" / "fracessa"
DEFAULT_RESULTS = Path(__file__).with_name("results.csv")

PARAMETERS = ("u1", "u2", "u3", "u4", "v0", "v1", "v2", "v3", "v4")
OFFSET_PARAMETERS = ("v1", "v2", "u3", "v4", "v3", "u2", "v1", "v0", "u1", "v2", "v3", "u4")
BASE = (700, 700, 700, 700, 700, 1500, 1500, 1500, 1500)
BASELINE_ESS = 15_120
FIELDNAMES = (
    "timestamp",
    "search_index",
    "phase",
    "label",
    "coefficients",
    "matrix",
    "status",
    "ess_count",
    "elapsed_ns",
    "safe_fallback",
    "wall_seconds",
)


@dataclass(frozen=True, slots=True)
class Result:
    status: str
    ess_count: int | None
    elapsed_ns: int | None
    safe_fallback: str | None
    wall_seconds: float


def matrix_string(coefficients: tuple[int, ...]) -> str:
    values = dict(zip(PARAMETERS, coefficients, strict=True))
    return "24#" + ",".join(str(values[name]) for name in OFFSET_PARAMETERS)


def canonical(coefficients: tuple[int, ...]) -> tuple[int, ...]:
    common = 0
    for value in coefficients:
        common = gcd(common, abs(value))
    if common <= 1:
        return coefficients
    return tuple(value // common for value in coefficients)


def parse_output(stdout: str, wall_seconds: float) -> Result:
    lines = stdout.splitlines()
    if len(lines) < 3:
        raise ValueError(f"unexpected FracESSA output: {stdout!r}")
    ess_count = int(lines[0])
    elapsed_ns = int(lines[1])
    safe_fallback = None if lines[2] == "null" else lines[2]
    return Result("ok", ess_count, elapsed_ns, safe_fallback, wall_seconds)


def run_matrix(binary: Path, cpu: int, matrix: str, timeout: float) -> Result:
    command = ["taskset", "-c", str(cpu), str(binary), "fast", matrix, "--timing"]
    started = time.monotonic()
    try:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return Result("timeout", None, None, None, time.monotonic() - started)
    wall_seconds = time.monotonic() - started
    if completed.returncode:
        message = completed.stderr.strip().replace("\n", " ")
        return Result(f"error:{completed.returncode}:{message}", None, None, None, wall_seconds)
    return parse_output(completed.stdout, wall_seconds)


def initial_variants():
    yield "baseline", "baseline", BASE
    for index, name in enumerate(PARAMETERS):
        for delta in (-1, 1):
            values = list(BASE)
            values[index] += delta
            yield "one_coordinate", f"{name}{delta:+d}", tuple(values)


def extension_variants():
    for step in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
        for index, name in enumerate(PARAMETERS):
            for delta in (-step, step):
                values = list(BASE)
                values[index] += delta
                yield "coordinate_sweep", f"{name}{delta:+d}", tuple(values)

    for step in (1, 4, 16, 64, 256):
        for first, second in combinations(range(len(PARAMETERS)), 2):
            for sign_first, sign_second in product((-1, 1), repeat=2):
                values = list(BASE)
                values[first] += sign_first * step
                values[second] += sign_second * step
                label = f"{PARAMETERS[first]}{sign_first * step:+d},{PARAMETERS[second]}{sign_second * step:+d}"
                yield "pair_sweep", label, tuple(values)


def random_deltas(seed: int):
    generator = random.Random(seed)
    amplitudes = (1, 1, 2, 2, 4, 4, 8, 16, 32, 64, 256, 2048)
    index = 0
    while True:
        amplitude = amplitudes[index % len(amplitudes)]
        changed = generator.randint(2, len(PARAMETERS))
        positions = generator.sample(range(len(PARAMETERS)), changed)
        deltas_by_position = [0] * len(PARAMETERS)
        deltas = []
        for position in positions:
            delta = 0
            while delta == 0:
                delta = generator.randint(-amplitude, amplitude)
            deltas_by_position[position] = delta
            deltas.append(f"{PARAMETERS[position]}{delta:+d}")
        yield "random", ",".join(deltas), tuple(deltas_by_position)
        index += 1


def read_progress(path: Path) -> tuple[set[tuple[int, ...]], int, tuple[int, ...], int]:
    tested: set[tuple[int, ...]] = set()
    best_count = BASELINE_ESS
    best_coefficients = BASE
    search_index = 0
    if not path.exists():
        return tested, best_count, best_coefficients, search_index
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                coefficients = tuple(json.loads(row["coefficients"]))
                count = int(row["ess_count"]) if row["ess_count"] else None
                index = int(row["search_index"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                continue
            tested.add(canonical(coefficients))
            search_index = max(search_index, index + 1)
            if count is not None and count >= best_count:
                best_count = count
                best_coefficients = coefficients
    return tested, best_count, best_coefficients, search_index


def write_row(writer: csv.DictWriter, search_index: int, phase: str, label: str, coefficients: tuple[int, ...], result: Result) -> None:
    writer.writerow(
        {
            "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "search_index": search_index,
            "phase": phase,
            "label": label,
            "coefficients": json.dumps(coefficients, separators=(",", ":")),
            "matrix": matrix_string(coefficients),
            "status": result.status,
            "ess_count": "" if result.ess_count is None else result.ess_count,
            "elapsed_ns": "" if result.elapsed_ns is None else result.elapsed_ns,
            "safe_fallback": "" if result.safe_fallback is None else result.safe_fallback,
            "wall_seconds": f"{result.wall_seconds:.9f}",
        }
    )


def self_check() -> None:
    assert matrix_string(BASE) == "24#1500,1500,700,1500,1500,700,1500,700,700,1500,1500,700"
    assert canonical(BASE) == (7, 7, 7, 7, 7, 15, 15, 15, 15)
    variants = list(initial_variants())
    assert len(variants) == 19
    assert len({values for _, _, values in variants}) == 19
    print("self-check passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hours", type=float, default=0.0, help="stop after this many hours; zero runs until interrupted")
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--initial-only", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return 0
    if not args.binary.is_file():
        parser.error(f"binary not found: {args.binary}")

    tested, best_count, best_coefficients, search_index = read_progress(args.results)
    args.results.parent.mkdir(parents=True, exist_ok=True)
    new_file = not args.results.exists() or args.results.stat().st_size == 0
    deadline = None if args.hours == 0 else time.monotonic() + args.hours * 3600
    last_sync = time.monotonic()
    completed_since_start = 0

    with args.results.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        if new_file:
            writer.writeheader()

        def evaluate(phase: str, label: str, coefficients: tuple[int, ...]) -> None:
            nonlocal best_count, best_coefficients, search_index, last_sync, completed_since_start
            key = canonical(coefficients)
            if key in tested or (deadline is not None and time.monotonic() >= deadline):
                return
            result = run_matrix(args.binary, args.cpu, matrix_string(coefficients), args.timeout_seconds)
            write_row(writer, search_index, phase, label, coefficients, result)
            tested.add(key)
            completed_since_start += 1
            search_index += 1

            if result.ess_count is not None and result.ess_count > best_count:
                best_count = result.ess_count
                best_coefficients = coefficients
                print(f"NEW RECORD ess={best_count} coefficients={coefficients} matrix={matrix_string(coefficients)}", flush=True)
            elif result.ess_count == best_count:
                best_coefficients = coefficients

            now = time.monotonic()
            if now - last_sync >= 10:
                stream.flush()
                os.fsync(stream.fileno())
                last_sync = now
            if completed_since_start % 25 == 0:
                remaining = "unlimited" if deadline is None else f"{max(0.0, deadline - now) / 3600:.2f}h"
                print(f"tested={len(tested)} best={best_count} remaining={remaining}", flush=True)

        try:
            for phase, label, coefficients in initial_variants():
                evaluate(phase, label, coefficients)
            if not args.initial_only:
                for phase, label, coefficients in extension_variants():
                    evaluate(phase, label, coefficients)
                    if deadline is not None and time.monotonic() >= deadline:
                        break
                for phase, label, deltas in random_deltas(args.seed):
                    coefficients = tuple(value + delta for value, delta in zip(best_coefficients, deltas, strict=True))
                    evaluate(phase, label, coefficients)
                    if deadline is not None and time.monotonic() >= deadline:
                        break
        except KeyboardInterrupt:
            print("interrupted; checkpointing", flush=True)
        finally:
            stream.flush()
            os.fsync(stream.fileno())

    print(f"finished tested={len(tested)} best={best_count} coefficients={best_coefficients}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
