# Agent Documentation Index

## Startup

1. Read repository-root `AGENTS.md`.
2. Read `KNOWLEDGE.md`.
3. Read `reviews/CPP_REVIEW.md` before reviewing or changing C++ correctness,
   hot paths, parsing, CMake, CTest, or release behavior.
4. Read `reviews/PYBIND_REVIEW.md` before reviewing or changing the native
   `fracessa_core` boundary, result/status conversion, GIL behavior, or native
   timing.
5. Read `reviews/PYTHON_REVIEW.md` before reviewing or changing the Python API,
   multiprocessing, sinks, verification, or benchmarks.
6. Read `CHANGES.md` only when history matters, when diagnosing drift, or when
   the user asks what changed.
7. Read only the task-specific references below.

## Maintained Documents

- `KNOWLEDGE.md`: canonical current architecture, workflows, constraints, and
  project policy.
- `CHANGES.md`: append-only human-readable history of meaningful project work.
- `reviews/CPP_REVIEW.md`: unresolved C++ correctness, speed, build, and release
  findings.
- `reviews/PYBIND_REVIEW.md`: unresolved native Python binding and boundary
  findings.
- `reviews/PYTHON_REVIEW.md`: unresolved Python correctness, speed, API, and
  verification findings.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md`: exact derivation of the removed
  double-PD certificate bug, its regression games, and arbitrary-small-
  perturbation counterexamples.
- `correctness/CERTIFIED_CANDIDATE_FILTER.md`: unimplemented design for exact
  affine normalization, rigorous numerical rejection certificates, and exact
  fallback in candidate enumeration.
- `architecture/UNSAFE_CANDIDATE_FILTER.md`: implemented temporary no-flag
  unsafe filter and parser-rename design, including normalization, its cheap
  danger veto, exact fallback, known limits, and verification requirements.
- `architecture/CHOICE_ONE_CANDIDATE_FILTER.md`: future-phase mathematical plan
  for adding certified Choice 1 after the unsafe phase is accepted; its source
  scope must be re-audited against that accepted implementation.
- `architecture/SUPPORT_GENERATOR_HANDOVER.md`: durable handover of the agreed
  generator API, callback rationale, candidate lifecycle, circular multiplier,
  experimental V2 status, and deferred decisions.
- `python-wrapper/README.md`: public Python and multiprocessing API.
- `plans/SUPPORT_GENERATORS.md`: implemented DFS/FKM support-generator design,
  circular bracelet reduction, orbit-aware exact-candidate pruning, correctness
  proofs, experimental evidence, and future benchmark gates.
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
