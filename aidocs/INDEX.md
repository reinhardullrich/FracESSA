# Agent Documentation Index

This file routes to every other Markdown document under `aidocs/`. Each entry states the document's role and status without
duplicating its contents.

## Current Overview And Operations

- `PROJECT.md`: current FracESSA structure, product behavior, architecture, dependencies, test-data roles, Python boundary, and
  release model.
- `RELEASING.md`: current calendar-version release procedure, supported artifacts, and manual GitHub Actions workflow.
- `pyfracessa/README.md`: current public PyFracESSA API, multiprocessing behavior, result schema, sinks, timing tool, and examples.

## Current Design And Technical Reference

- `architecture/SUPPORT_GENERATOR_HANDOVER.md`: implemented generator API, callback rationale, candidate lifecycle, circular
  multiplier contract, and retained design decisions.
- `plans/SUPPORT_GENERATORS.md`: implemented DFS and direct-bracelet algorithms, correctness arguments, rejected alternatives,
  measurements, and future benchmark gates.
- `plans/CHEAP_CYCLIC_SYMMETRY_FILTER.md`: implemented affine cyclic-symmetry reduction, verification, and promotion evidence.
- `plans/EXACT_STABILITY_SCHUR_COMPLEMENT.md`: implemented exact Schur-complement stability reduction and verification plan.
- `plans/MULTIWORD_BIT_ARRAY.md`: implemented dimension-64 and multiword support architecture, staged verification, and performance
  evidence for preserving the one-word path.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md`: exact analysis of the removed binary64 positive-definiteness certificate bug.
- `correctness/FAST_CANDIDATE_FALSE_REJECTION.md`: exact ESS counterexamples for former fast rejection rules and current mitigations.
- `reference/ELIMINATING_THE_BORDERED_CANDIDATE_SYSTEM.md`: mathematical derivation of the reduced symmetric candidate system.
- `reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`: current production-only bit-scanning call chain.
- `reference/INTEGER_HADELER_COPOSITIVITY.md`: current proof, exact algorithm, branch tests, and benchmark evidence for the optimized
  integer Hadeler checker.

## Research And Open Options

- `plans/MAJOR_SINGLE_CORE_PERFORMANCE_OPPORTUNITIES.md`: unimplemented research review of possible material single-core speedups.

## Historical Decisions And Audits

- `CHANGES.md`: append-only searchable history of meaningful decisions, results, and evidence.
- `architecture/UNSAFE_CANDIDATE_FILTER.md`: failed temporary normalized-filter design and retained measurements.
- `history/CPP_PERFORMANCE_EXPERIMENTS_2026-08-03.md`: completed C++ allocation and performance experiments, especially rejected
  changes.
- `history/FAST_PIPELINE_EXPERIMENTS_2026-08-03.md`: completed and rejected binary64 candidate-path experiments.
- `history/INTEGER_STABILITY_MATRIX_MIGRATION_ONLY_2026-08-06.md`: completed integer-storage migration plan.
- `reference/MATRIX_GENERATOR_CATALOGUE_AUDIT.md`: immutable matrix-generator catalogue import audit; its counts are not live database
  counts.

## Frozen Experiments

- `experiments/HISTORICAL_DEFAULT_VERY_UNSAFE_COMPARISON_2026-07-31.md`: historical default, very-unsafe, unsafe, verified, and exact
  timing comparison.
- `experiments/SUPPORT_FRONTIER_2026-07-29.md`: support-generation and pruning alternatives, measurements, and conclusions.
- `experiments/speed_comparison_2026-07-26/COMPUTER_TIMING_COMPARISON.md`: frozen cross-computer timing comparison.
- `experiments/speed_comparison_2026-07-26/EXACT_PD_ONLY_COMPARISON.md`: frozen exact-positive-definiteness experiment.
- `experiments/speed_comparison_2026-07-26/EXACT_VS_DEFAULT_COMPARISON.md`: frozen exact-versus-default benchmark.
- `experiments/speed_comparison_2026-07-26/MICROBENCHMARK_COMPARISON.md`: frozen persistent-process exact-PD microbenchmark.
