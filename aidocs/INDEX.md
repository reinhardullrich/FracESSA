# Agent Documentation Index

## Startup

1. Read repository-root `AGENTS.md`.
2. Read `KNOWLEDGE.md`.
3. Read `reviews/CPP_REVIEW.md` before reviewing or changing C++ correctness,
   hot paths, parsing, CMake, CTest, or release behavior.
4. Read `reviews/PYTHON_REVIEW.md` before reviewing or changing the Python API,
   multiprocessing, sinks, verification, or benchmarks.
5. Read `CHANGES.md` only when history matters, when diagnosing drift, or when
   the user asks what changed.
6. Read only the task-specific references below.

## Maintained Documents

- `KNOWLEDGE.md`: canonical current architecture, workflows, constraints, and
  project policy.
- `CHANGES.md`: append-only human-readable history of meaningful project work.
- `reviews/CPP_REVIEW.md`: unresolved C++ correctness, speed, build, and release
  findings.
- `reviews/PYTHON_REVIEW.md`: unresolved Python correctness, speed, API, and
  verification findings.
- `python-wrapper/README.md`: public Python and multiprocessing API.
- `reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`: compact production-only
  bit-scanning call chain requested for hot-path work.

## Fixed References

- `experiments/`: dated benchmark reports. These are immutable snapshots, not
  statements about the current worktree.
- `reference/papers/`: mathematical source papers used by the project.

## Maintenance Rules

- Current code and tests override documentation when they disagree.
- Update the owning maintained document in the same task as a behavior change.
- Append to `CHANGES.md` when a human-readable historical record is useful;
  do not rewrite its previous entries.
- Keep current facts out of dated experiment reports.
- Do not add session histories, generated source dumps, completed reviews,
  or stale profiling output.
