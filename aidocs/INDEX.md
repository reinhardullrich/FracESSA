# Agent Documentation Index

Use this file only to locate the document that owns a topic. It lists every other Markdown file under `aidocs/` and labels whether
the file describes current behavior, open research, or history.

## Current Overview And Operations

- `PROJECT.md`: concise current project structure, product contract, architecture, dependencies, interfaces, and verification model.
- `RELEASING.md`: current release procedure, artifact set, and manual GitHub Actions workflow.
- `pyfracessa/README.md`: current internal persistent-Pybind timing command, calibration rules, provenance, and report contract.

## Current Architecture And Mathematics

- `architecture/SUPPORT_GENERATORS.md`: implemented non-circular DFS and circular bracelet generators, pruning proof, callback
  interface, rejected alternatives, and retained performance evidence.
- `architecture/CYCLIC_AFFINE_SYMMETRY.md`: implemented exact affine symmetry detection and reduction for circular games.
- `architecture/MULTIWORD_SUPPORT_MASKS.md`: implemented one-word/multiword boundary, representation, dispatch, and verification.
- `correctness/FAST_CANDIDATE_FALSE_REJECTION.md`: exact counterexamples and the non-certifying correctness boundary of `fast`.
- `reference/ELIMINATING_THE_BORDERED_CANDIDATE_SYSTEM.md`: mathematical derivation of the reduced symmetric candidate system.
- `reference/EXACT_STABILITY_SCHUR_COMPLEMENT.md`: exact Schur-complement derivation and implemented reduced-B stability path.

## Open Research

- `plans/MAJOR_SINGLE_CORE_PERFORMANCE_OPPORTUNITIES.md`: active plan of unimplemented hypotheses for material single-core speedups,
  with completed foundations clearly separated.

## Historical Decisions And Implementation Records

- `CHANGES.md`: append-only searchable history of decisions, results, and evidence not obvious from Git.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md`: proof and counterexamples for the removed binary64 positive-definiteness certificate.
- `history/CPP_PERFORMANCE_EXPERIMENTS_2026-08-03.md`: measured allocation/performance experiments, especially rejected changes.
- `history/FAST_PIPELINE_EXPERIMENTS_2026-08-03.md`: promoted and rejected binary64 candidate-path experiments.
- `history/INTEGER_STABILITY_MATRIX_MIGRATION_ONLY_2026-08-06.md`: dated rational-to-integer stability migration record.
- `history/LOGGING_AND_MATRIX_FORMATTING_2026-08-09.md`: dated diagnostic-log redesign and verification record.
- `history/UNSAFE_CANDIDATE_FILTER_2026-07-27.md`: retired normalized unsafe-filter design, measurements, and failure outcome.
- `plans/COPOSIT_INTEGRATION_AND_TYPE_OWNERSHIP.md`: completed coposit/FracESSA ownership migration and measured implementation record.
- `reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`: historical call chain through the removed FracESSA-owned Hadeler checker.
- `reference/INTEGER_HADELER_COPOSITIVITY.md`: proof, implementation, branch tests, and evidence for the removed FracESSA-owned
  exact integer Hadeler checker.
- `reference/MATRIX_GENERATOR_CATALOGUE_AUDIT.md`: immutable provenance and sampling audit for the matrix-generator import; its counts
  are not live database counts.

## Frozen Benchmark Reports

- `experiments/HISTORICAL_DEFAULT_VERY_UNSAFE_COMPARISON_2026-07-31.md`: historical unsafe/default/verified/exact timing and
  correctness comparison.
- `experiments/SUPPORT_FRONTIER_2026-07-29.md`: support-frontier, Gosper, DFS, streaming, pruning, timing, and memory experiments.
- `experiments/speed_comparison_2026-07-26/COMPUTER_TIMING_COMPARISON.md`: frozen cross-computer timing comparison.
- `experiments/speed_comparison_2026-07-26/EXACT_PD_ONLY_COMPARISON.md`: frozen exact-positive-definiteness experiment.
- `experiments/speed_comparison_2026-07-26/EXACT_VS_DEFAULT_COMPARISON.md`: frozen exact-versus-default benchmark.
- `experiments/speed_comparison_2026-07-26/MICROBENCHMARK_COMPARISON.md`: frozen persistent-process exact-PD microbenchmark.
