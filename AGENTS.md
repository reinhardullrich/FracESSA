# AGENTS.md

**Use the Ponytail skill for every code analysis, plan, modification, review, and summary.**

## Startup

1. Read `aidocs/PROJECT.md` before work that depends on current project behavior, architecture, structure, or workflows.
2. Use `aidocs/INDEX.md` to locate only the detailed documents needed for the task.
3. Search `aidocs/CHANGES.md` only when investigating history, drift, an old name or behavior, or a previous experiment. Never read it
   routinely or sequentially.

## Worktree

1. Work in the main worktree at `/home/reinhard/projects/fracessa` unless Reinhard explicitly approves another worktree first.
2. Do not switch to, run project commands in, or modify another worktree without that prior approval.

## Source-Code Approval

1. Do not modify any C++ or Python source file (`*.cpp`, `*.hpp`, or `*.py`) without Reinhard's explicit approval for the exact files
   and intended changes.
2. Before editing, work read-only, name every source file that would change, and describe the smallest necessary patch. A general
   request to investigate or fix an issue does not approve additional files.
3. Approval is limited to the named files and scope. Stop and ask again if another source file or broader change becomes necessary.
4. Make only the approved change. Do not include opportunistic cleanup, refactoring, formatting, or adjacent fixes.
5. After editing, show the actual source changes and verification results for Reinhard's direct review before making further source
   changes.

## Priorities

1. Correctness is absolute. Speed is second. Other concerns are secondary.
2. FracESSA may inspect a `2^n` support space and execute exact numerical kernels millions or billions of times. Optimize for many
   operations on small matrices, not a few operations on large matrices.
3. Keep validation at input boundaries. Do not add checks, allocations, abstractions, or branches to proven hot paths without a
   correctness need and a representative benchmark.
4. Intentionally unchecked bitset operations exist for speed. Do not add defensive hot-path checks for states excluded by validated
   input and algorithm contracts.
5. Keep an allocation or memory optimization only when end-to-end benchmarks show a repeatable speed improvement.
6. Performance evidence must cover the affected methods, small and large games, and circular and non-circular matrices. Include
   dimensions around 20 and 23 when feasible. A synthetic kernel benchmark alone is insufficient.
7. Keep dimension-2 matrices in correctness data, but exclude them from performance runs, tables, and aggregate statistics.
8. Keep numerical code human-readable. Do not add custom allocators, pools, generic workspace frameworks, or extensive plumbing for
   small theoretical gains.
9. FracESSA orchestration code uses C++ numerical wrappers. Raw FLINT types and `fmpz_*`/`fmpq_*` calls belong only in thin `linalg`
   wrappers and specialized numerical kernels.

## Style

1. Use 140 columns as a soft line-width target for C++, Python, comments, and Markdown.
2. Do not wrap lines merely to satisfy the traditional 80-column limit.
3. Exceed 140 columns when splitting a formula, command, URL, matrix, or readable expression would make it harder to understand.

## Build And Verification

Run the normal Release build from the repository root:

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build --parallel
```

`BUILD_TESTING` defaults to `ON`; use `-DBUILD_TESTING=OFF` only when tests are intentionally unnecessary. Do not reproduce CMake's
optimization flags by hand.

Select the project-local optimized FLINT 3.6 installation with:

```bash
cmake -S cpp -B cpp/build-flint36 -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
```

Use project-local `.local-tmp/` for substantial temporary build data because this machine's `/tmp` is memory-backed. If sandboxing
blocks the normal ccache directory, request the required filesystem access instead of disabling or redirecting ccache.

Run the maintained C++/CLI and Python suites with:

```bash
ctest --test-dir cpp/build --output-on-failure --parallel
PYTHONPATH=python python3 -m unittest discover -s python/tests -p 'test_*.py'
```

Regenerate the combined C++/Python documentation with the project environment:

```bash
doxygen docs/Doxyfile
python/.venv/bin/sphinx-build -W --keep-going -b html docs site
```

Check the main database with:

```bash
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA integrity_check;'
sqlite3 testdata/fracessa_testdata.sqlite3 'PRAGMA foreign_key_check;'
```

## Performance Measurements

Unless Reinhard explicitly requests a compiler experiment, use this canonical build:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S cpp -B cpp/build-benchmark -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DFRACESSA_NATIVE_ARCH=ON
cmake --build cpp/build-benchmark --target fracessa fracessa_core --parallel
```

Build both compared revisions afresh with the same compiler, CMake, generator, dependencies, and command. Canonical local timing uses
CPU 2, one persistent Pybind process, stored per-matrix calibrations, a 0.5-second default target, and median native nanoseconds.
Record and label every deviation. Never compare CLI process timing with persistent-Pybind timing. Follow `aidocs/pyfracessa/README.md`
for the maintained timing command and report format.

## Runtime And Release

1. Do not add a per-matrix computation timeout to the CLI, native API, Python API, or worker-liveness handling. A valid matrix may run
   for hours.
2. Follow `aidocs/RELEASING.md` exactly for releases; do not duplicate or improvise the release procedure.

## Documentation

1. Keep project instructions and essential commands in this file, current project facts in `aidocs/PROJECT.md`, and every detailed
   agent document listed in `aidocs/INDEX.md`.
2. Add detailed documents only when needed. Do not duplicate their contents in `aidocs/PROJECT.md` or `aidocs/INDEX.md`.
3. Update the owning document in the same task as a documented behavior change, and update `aidocs/INDEX.md` whenever its inventory
   or navigation changes.
4. Add a `CHANGES.md` entry only when rationale, results, or evidence would not be clear from the Git diff and commit message. Routine
   edits belong only in Git.
5. Do not store session transcripts, generated source concatenations, stale profiling output, historical benchmark tables, or
   database row inventories in `AGENTS.md` or `aidocs/PROJECT.md`.
