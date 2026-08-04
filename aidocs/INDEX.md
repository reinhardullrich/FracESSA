# Agent Documentation Index

Use each document according to its role below. A file under `plans/` or
`architecture/` may describe implemented behavior or retained historical
rationale; its status line and this index determine which.

## Startup

1. Read repository-root `AGENTS.md`.
2. Read `KNOWLEDGE.md`.
3. Read `CHANGES.md` only when history matters, when diagnosing drift, or when
   the user asks what changed.
4. Read only the task-specific references below.

## Current Truth

- `KNOWLEDGE.md`: canonical current architecture, workflows, constraints, and
  project policy.
- `pyfracessa/README.md`: public PyFracESSA and multiprocessing API.
- `../testdata/README.md`: canonical SQLite test-data and timing schema.

## Current Design And Technical Reference

- `architecture/SUPPORT_GENERATOR_HANDOVER.md`: concise current generator API,
  callback rationale, candidate lifecycle, circular multiplier contract, and
  retained design decisions.
- `plans/SUPPORT_GENERATORS.md`: implemented DFS/direct-bracelet generator
  design, correctness arguments, rejected alternatives, measurements, and
  future benchmark gates.
- `plans/CHEAP_CYCLIC_SYMMETRY_FILTER.md`: implemented always-on affine cyclic
  symmetry reduction plus its retained implementation and promotion record.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md`: exact derivation of the removed
  double-PD certificate bug, its regression games, and arbitrary-small-
  perturbation counterexamples.
- `correctness/FAST_CANDIDATE_FALSE_REJECTION.md`: exact ESS counterexamples for
  former fast per-support rejection rules and the current mitigations.
- `reference/ELIMINATING_THE_BORDERED_CANDIDATE_SYSTEM.md`: introductory
  mathematical derivation of the equivalent symmetric $(k-1)\times(k-1)$
  candidate system, including reconstruction and singularity equivalence.
- `reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`: compact production-only
  bit-scanning call chain requested for hot-path work.

## Research And Open Options

- `plans/MAJOR_SINGLE_CORE_PERFORMANCE_OPPORTUNITIES.md`: research review of
  theorem-backed and algorithmic routes to material single-core speedups; no
  listed idea is implemented or benchmarked by that document.
- `../research/AUTOMORPHISM_GROUPS_AND_GAMMA.md`: exact symmetry analysis of the published dimension-24 high-ESS rook game and
  related search directions.
- `../research/ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md`: rectangular-rook group orders, compact symbolic matrices, orbit ceilings,
  empirical ESS counts, and Sperner bounds through dimension 30.
- `../research/THREE_BY_K_ROOK_ESS_SEQUENCES.md`: derivation, exact finite verification, support counts, and gamma values for the
  two compact $3\times k$ rook ESS sequences.

## Historical Records And Fixed Evidence

- `CHANGES.md`: append-only human-readable history of meaningful project work.
- `history/CPP_PERFORMANCE_EXPERIMENTS_2026-08-03.md`: completed C++ allocation and implementation experiments, especially
  measured changes that were rejected.
- `history/FAST_PIPELINE_EXPERIMENTS_2026-08-03.md`: completed and rejected binary64 candidate-path experiments.
- `architecture/UNSAFE_CANDIDATE_FILTER.md`: failed temporary normalized-filter
  phase, including why it was removed and the measurements worth retaining.
- `reference/MATRIX_GENERATOR_CATALOGUE_AUDIT.md`: immutable catalogue-import
  audit snapshot; its counts describe that audit, not the live database.
- `experiments/`: dated benchmark reports. These are immutable snapshots, not
  statements about the current worktree.
- `../research/papers/`: mathematical source papers used by the project.

## Maintenance Rules

- Current code and tests override documentation when they disagree.
- Update the owning maintained document in the same task as a behavior change.
- Append to `CHANGES.md` when a human-readable historical record is useful;
  do not rewrite its previous entries.
- Dated files under `history/` contain no open tasks and are never current architecture or validation records.
- Keep current facts out of dated experiment reports.
- Do not add session histories, generated source dumps, duplicate completed
  reviews, or stale profiling output.
