# CHANGES

Append-only release-notes style log. Add new entries; do not edit historical entries.

## 2026-02-11
1. Optimized hot-path support iteration in `fracessa/cpp/src/findeq.cpp`:
used fixed stack buffers (`uint8_t[64]`) to extract support and non-support indices once per support, then iterated index arrays directly in row-sum loops.
2. Removed repeated inner-loop bitset membership checks in candidate filtering:
`find_candidate_dbl` and `find_candidate_frc` now avoid branchy `is_set_at_pos` checks inside matrix accumulation loops.
3. Optimized matrix extraction routines in `fracessa/cpp/include/fracessa/matrix_server.hpp`:
`get_linear_system_frc`, `get_linear_system_dbl`, `get_bee_matrix_frc`, and `get_bee_matrix_dbl` now use pre-extracted fixed-size index arrays instead of full-range scans with per-element bit checks.
4. Updated project knowledge:
documented fixed-index-buffer hot-path pattern and added `CHANGES.md` append-only workflow note in `PROJECT_KNOWLEDGE.md`.
5. Full-solution construction now explicitly uses extracted fixed stack support partitions:
added helpers in `fracessa/cpp/src/findeq.cpp` to build full solutions from support/non-support index buffers for both double and rational paths.
6. Centralized support-index extraction in `bs64`:
moved `extract_set_indices` into `fracessa/cpp/include/fracessa/bitset64.hpp` (`bs64` namespace) and reused it from `findeq.cpp` and `matrix_server.hpp` to keep hot-path indexing logic in one place.
7. Inlined full-solution fill logic in candidate functions:
removed `fill_full_solution_dbl` / `fill_full_solution_frc` helper functions and embedded the same fixed-buffer fill code directly in `find_candidate_dbl` and `find_candidate_frc` per hot-path style preference.
8. Baseline regeneration script normalization fix:
updated `fracessa/python/verification/create_baselines.py` to normalize legacy reason code values to canonical current strings (`T_pd_dbl` instead of `T_pd_double`) while staying backward-compatible.
9. Baseline regeneration scope change:
updated `fracessa/python/verification/create_baselines.py` to always process all matrices from `verification_matrices.json`, including entries with `in_use=false`.
10. Unused-code quarantine and test cleanup:
commented out production-unused C++ helpers in `fracessa/cpp/include/fracessa/bitset64.hpp`, `fracessa/cpp/include/linalg/fraction.hpp`, `fracessa/cpp/include/linalg/matrix_fraction.hpp`, `fracessa/cpp/include/fracessa/supports.hpp`, and commented unused members in `fracessa/cpp/include/fracessa/fracessa.hpp` (+ matching ctor init cleanup in `fracessa/cpp/src/fracessa.cpp`).
11. Removed tests that only covered non-production-only APIs:
deleted `clear_bit_at_pos`/`next_bitset_with_same_popcount`/`hash` test cases from `fracessa/cpp/tests/test_bitset64.cpp` and `is_one`/`negate_inplace`/`abs_inplace` test cases from `fracessa/cpp/tests/test_fraction.cpp`.
12. Minor test/code hygiene:
removed an unused local in `fracessa/cpp/tests/test_integration.cpp` and removed an unused exception variable in `fracessa/cpp/src/main.cpp`.
13. GitHub branch synchronization:
force-pushed local `main` to `origin/main`, switched remote URL to SSH, deleted remote `refactor-remove-templates`, and deleted local `refactor-remove-templates` so only `main` remains.
14. Script relocation to tracked repo root:
moved `build.sh`, `test.sh`, and `python.sh` from workspace root into `fracessa/` so they are versioned on GitHub, and refactored them to run from script directory with repo-local paths (`cpp/`, `python/`, `build/`).
15. Branding and README alternative:
replaced `fracessa/README.md` with a clearer workflow-oriented version and added a new GitHub-ready vector logo at `fracessa/logo.svg` (README now references `logo.svg`).
16. Logo naming clarification:
updated `fracessa/logo.svg` subtitle text to explicitly expand FracESSA as `Fractional ESS Analyzer`.
17. Final logo switch to user-provided artwork:
set `fracessa/logo.jpeg` (uploaded image) as the active README logo by changing `fracessa/README.md` to reference `logo.jpeg` instead of `logo.svg`.
18. Removed unused function code from C++:
deleted dead/disabled function blocks from `fracessa/cpp/include/fracessa/bitset64.hpp`, `fracessa/cpp/include/linalg/fraction.hpp`, `fracessa/cpp/include/linalg/matrix_fraction.hpp`, and cleanup remnants in `fracessa/cpp/include/fracessa/supports.hpp`, `fracessa/cpp/include/fracessa/fracessa.hpp`, `fracessa/cpp/src/fracessa.cpp`, plus removed disabled test-function blocks in `fracessa/cpp/tests/test_fraction.cpp` and `fracessa/cpp/tests/test_matrix.cpp`.
19. Expanded explanatory comments across all C++ headers/sources:
added or refined project-level and algorithm-level comments in `fracessa/cpp/include/fracessa/*.hpp`, `fracessa/cpp/include/linalg/*.hpp`, `fracessa/cpp/src/*.cpp`, and `fracessa/cpp/tests/*.cpp` to make pipeline flow, mathematical criteria, and performance-oriented design choices easier to understand for mathematically experienced contributors new to this codebase.
20. Added fast/full verification workflow to project test scripts:
`fracessa/test.sh` now accepts `--full`; default mode runs build + CTest + Python verification in fast mode (`--fast-cutoff 2.0`), while `--full` runs all verification matrices.
21. Updated Python verification runner behavior for CTest-style gating:
`fracessa/python/run_matrices.py` now supports `--full` and `--fast-cutoff`, uses timing-based matrix selection in fast mode (prefers latest local run timings, falls back to baseline timings), includes deterministic `candidate_id` in candidate comparison keys, and exits non-zero on verification failure.
22. Fixed Python wrapper executable path detection after script relocation:
`fracessa/python/fracessa_py.py` now prefers nested repo binary path `fracessa/build/fracessa` and keeps workspace-level `build/fracessa` as compatibility fallback.
23. Script passthrough update:
`fracessa/python.sh` now forwards CLI arguments to `python/run_matrices.py`.
24. Refactored correctness gating into CTest:
added per-matrix correctness verifier `fracessa/python/verification/ctest_verify_matrix.py` and registered one CTest per verification matrix (`VerificationMatrix_<id>`) from `cpp/tests/CMakeLists.txt`, enabling parallel matrix correctness execution.
25. Added matrix selection utility for fast/full correctness modes:
introduced `fracessa/python/verification/matrix_selection.py` to emit test matrix IDs and fast subsets (`<=2s` reference timing with latest-results preference and baseline fallback), and to expose all IDs for full mode.
26. Reworked `fracessa/test.sh` flow:
default now runs (a) core C++ CTests excluding matrix-verification tests and (b) fast matrix correctness CTests selected via `matrix_selection.py`, both with `ctest -j`; `--full` runs all matrix correctness CTests including `in_use=false` matrices.
27. Converted Python verification runner to speed-only role:
rewrote `fracessa/python/run_matrices.py` so it benchmarks and reports timing artifacts only; correctness pass/fail is no longer enforced there and is handled in CTest matrix verification tests.
28. Correctness comparison strictness details:
CTests compare candidate rows against baseline with deterministic `candidate_id` included and normalized floating formatting for `payoff_double` to avoid false mismatches from string formatting differences.
29. Fixed solver positivity bug affecting zero/negative-payoff candidates:
in `fracessa/cpp/include/linalg/linear_solver.hpp`, positivity checks now apply only to support-probability unknowns (`i < n-1`), while the final unknown (payoff) is allowed to be zero/negative for both double and rational solvers.
30. Renamed rational-PD reason label:
changed runtime stability label from `T_pd_rat` to `T_pd_frc` in `fracessa/cpp/src/checkstab.cpp` and updated baseline conversion logic in `fracessa/python/verification/create_baselines.py` (including legacy mapping from `T_pd_rat`).
31. Baseline candidate data normalization:
updated `fracessa/python/verification/baseline_candidates.csv` to replace `T_pd_rat` with `T_pd_frc`.
32. Reduced false negatives in correctness CTests for PD branch selection:
`fracessa/python/verification/ctest_verify_matrix.py` now treats `T_pd_dbl` and `T_pd_frc` as equivalent (`T_pd_*`) during comparison, while keeping all other reason labels strict.
33. Validation outcome after above fixes:
`./fracessa/test.sh` fast mode now passes 33/33 verification matrix tests (plus 7/7 core C++ tests).
34. Removed deprecated `in_use` matrix metadata:
deleted `in_use` from both verification datasets (`fracessa/python/verification/verification_matrices.json` and `fracessa/cpp/tests/test_data/verification_matrices.json`) so fast/full behavior is driven only by explicit matrix-ID policy.
35. Switched speed runner selection to static fast/full matrix IDs:
`fracessa/python/run_matrices.py` no longer uses timing cutoff selection (`--fast-cutoff` removed); fast mode now runs all matrices except hardcoded exclusions from `verification/matrix_selection.py` (currently IDs `32` and `34`), and full mode runs all matrices.
36. Baseline script message cleanup:
updated `fracessa/python/verification/create_baselines.py` wording to reflect all-matrix processing without referencing removed `in_use` state.
37. Canonical context update:
updated `PROJECT_KNOWLEDGE.md` to document the new static fast/full selection policy and removal of `in_use` from verification metadata.
38. Improved Python speed wrapper CLI UX:
`fracessa/python.sh` now mirrors `test.sh`-style help text (`-h/--help`), explicitly documents fast/full behavior, rejects unknown options with usage output, and forwards only supported mode (`--full`) to `python/run_matrices.py`.
39. Extracted matrix parsing into reusable module:
moved safe/unsafe matrix parsing logic out of `fracessa/cpp/src/main.cpp` into `fracessa/cpp/src/matrix_parser.cpp` with API in `fracessa/cpp/include/fracessa/matrix_parser.hpp`, and linked it through `fracessa_lib`.
40. Added parser unit tests:
introduced `fracessa/cpp/tests/test_matrix_parser.cpp` and `MatrixParserTests` CTest target to cover safe parser acceptance/rejection cases plus unsafe-parser dimension-bypass behavior.
41. Expanded linear solver edge-case coverage:
extended `fracessa/cpp/tests/test_linear_solver.cpp` with explicit regressions for (a) negative payoff acceptance, (b) non-positive support rejection, and (c) singular-system rejection for both double and fraction solvers.
42. Test suite status update:
non-verification C++ test suite now runs 8/8 tests (added `matrix_parser` test binary) while fast verification CTests remain 33/33 passing.
43. Added direct support-container tests:
introduced `fracessa/cpp/tests/test_supports.cpp` and `SupportsTests` CTest target to validate (a) support-bucket cardinality counts, (b) circular-symmetry canonicalization behavior, and (c) `remove_supersets` pruning semantics.
44. Added CS orbit/metadata integration regression:
extended `fracessa/cpp/tests/test_integration.cpp` with a circular-symmetric case (`5#1,3`) asserting deterministic orbit-expanded candidates and consistent `shift_reference` propagation.
45. Test suite status update:
non-verification C++ test suite now runs 9/9 tests (added `supports` test binary) with fast verification CTests still passing 33/33.
46. Strengthened bitset hot-path test coverage:
extended `fracessa/cpp/tests/test_bitset64.cpp` with targeted regressions for `extract_set_indices`, exact `rot_one_right` mapping behavior, masking of out-of-dimension bits during rotation, and top-bit/iteration boundary handling.
47. Strengthened matrix/LU negative-path test coverage:
extended `fracessa/cpp/tests/test_matrix.cpp` with double PD rejection cases (indefinite and non-square), and extended `fracessa/cpp/tests/test_lu.cpp` with singular-matrix regressions covering `isSingular`, zero determinant, `inverse()` throw, and `solve()` division-by-zero throw behavior.
48. Added CLI parser black-box CTest coverage:
introduced `fracessa/cpp/tests/cli_parser_blackbox.py` and wired `CliParserBlackBoxTests` in `fracessa/cpp/tests/CMakeLists.txt` to execute the built `fracessa` binary directly, covering safe parser success/failure behavior and `--unsafe` routing.
49. Test suite status update:
non-verification CTest suite now runs 10/10 tests (9 C++ test binaries + 1 CLI parser black-box script), with fast verification CTests still passing 33/33.
50. Removed unused duplicate matrix fixture file:
deleted `fracessa/cpp/tests/test_data/verification_matrices.json` because active test and verification flows use `fracessa/python/verification/verification_matrices.json` as the single source.
51. Added native pybind bridge:
introduced `fracessa/cpp/src/pybind_module.cpp` with `fracessa_core.compute_matrix(...)`, native status codes, candidate serialization, and GIL release during core compute.
52. Integrated pybind build into CMake:
added `pybind11` FetchContent dependency, built `fracessa_core` module from `cpp/CMakeLists.txt`, and linked it against `fracessa_lib`.
53. Enabled shared-module-safe linking for wrapper builds:
updated CMake math-library discovery to prefer shared GMP/MPFR on Unix, and enabled PIC (`CMAKE_POSITION_INDEPENDENT_CODE`) so `fracessa_core` links cleanly as a shared Python extension.
54. Implemented pybind-based Python wrapper package:
added `fracessa/python/wrapper_v1/` modules for typed jobs/config/results, native loader/adapter (`core.py`), single-run APIs (`single.py`), multiprocessing batch+queue execution (`mp.py`), and JSON/CSV/Arrow/Parquet/stream sinks.
55. Added wrapper contract document:
created `fracessa/python/wrapper_v1/API_CONTRACT.md` with concrete v1 API signatures, status codes, multiprocessing workflow, and sink/output options.
56. Updated canonical project context:
extended `PROJECT_KNOWLEDGE.md` to include pybind/wrapper components, dependency changes, wrapper smoke-test reality, and pybind-specific build caveats.
57. Added full wrapper user documentation:
created `fracessa/python/wrapper_v1/README.md` with build/import instructions, API walkthrough, single and multiprocessing examples (batch + queue), sink usage (`csv`, `json`, `arrow`, `parquet`, `stream`), and operational caveats.
58. Removed wrapper timeout API surface:
deleted `RunConfig.timeout_s`, removed `STATUS_TIMEOUT` from `fracessa_core` + status mapping, and updated wrapper docs/contracts so long-running matrices are explicitly supported without internal timeout controls.
59. Added wrapper Python test suite:
introduced `fracessa/python/wrapper_v1/tests/` (stdlib `unittest`) with unit coverage for types/core/io/sinks/single-run flow plus native integration checks for single and multiprocessing wrapper paths.
60. Integrated wrapper tests into project test script:
updated `fracessa/test.sh` to always run `python3 -m unittest discover -s python/wrapper_v1/tests -p \"test_*.py\"` between core CTests and matrix verification CTests.
61. Wrapper docs/test guidance update:
extended `fracessa/python/wrapper_v1/README.md` with a dedicated wrapper test-suite section and direct invocation command.
62. Canonical context update:
updated `PROJECT_KNOWLEDGE.md` to document wrapper-test execution in `test.sh`, direct wrapper-test command, no-timeout policy for wrapper jobs, and current wrapper unittest outcome.
63. Canonical context file rename:
renamed `PROJECT_KNOWLEDGE.md` to `KNOWLEDGE.md`, updated `AGENTS.md` startup/update references accordingly, and refreshed the knowledge document header/date to reflect the new canonical filename.
64. Documentation policy restructuring:
refactored `AGENTS.md` into a reusable non-project-specific template, explicitly required reading/updating both `KNOWLEDGE.md` and `CHANGES.md`, and moved documentation-ownership clarification into `KNOWLEDGE.md` so project specifics stay there.
65. Removed stale CMake test-data copy:
deleted the obsolete `cpp/tests/test_data` copy step from `fracessa/cpp/tests/CMakeLists.txt`; active matrix verification data remains centralized in `fracessa/python/verification/`, and `KNOWLEDGE.md` now states that policy explicitly.
66. Cleaned stale project-structure note:
removed the deleted `cpp/tests/test_data/verification_matrices.json` fixture from `analysis/PROJECT_STRUCTURE.md`; the documented active verification dataset remains `fracessa/python/verification/verification_matrices.json`.
67. Saved full repository review:
added `fracessa/REPO_REVIEW_2026-06-02.md` with prioritized findings from the full code/repo review, including multiprocessing backpressure, timeout policy gaps, release-build CPU targeting, CI test coverage, Windows integer parsing, tracked documentation placement, stale backup files, and multiprocessing config validation.
68. Fixed wrapper multiprocessing batch backpressure:
changed `python/wrapper_v1/mp.py` so `run_jobs_mp` interleaves job submission with result draining through a bounded submitted-but-not-yielded window (`min(queue_maxsize, workers * chunk_size)`), preventing bounded queue deadlocks for large matrix streams; added `MPConfig` value validation and unit coverage for the streaming submission behavior.
69. Removed active Python computation timeouts:
removed per-matrix subprocess timeouts from `python/fracessa_py.py`, `python/run_matrices.py`, `python/verification/ctest_verify_matrix.py`, and `python/verification/create_baselines.py`; updated baseline metadata and project knowledge to state that long-running matrices are externally managed rather than internally timed out.
70. Made release CPU targeting portable:
added `FRACESSA_NATIVE_ARCH` CMake option (default `ON` for local raw-speed builds) and configured GitHub release workflow builds with `-DFRACESSA_NATIVE_ARCH=OFF` so released Linux/macOS binaries do not require native CPU features from the runner.
71. Added release workflow test gate:
updated `.github/workflows/release.yml` to set up Python 3.12 and run core CTests, wrapper unittests, and fast verification matrix CTests after building and before publishing release artifacts.
72. Fixed platform-width rational parsing risk:
added a string-based `linalg::fraction` constructor using `fmpq_set_str`, changed safe and unsafe matrix parsers to construct CLI values from exact rational text instead of `long` casts, and added parser regressions for values above the Windows 32-bit `long` range.
73. Kept project memory files local-only:
confirmed `AGENTS.md`, `KNOWLEDGE.md`, and `CHANGES.md` remain workspace-root local files outside the nested `fracessa/` Git repo, removed the project-memory finding from the active repo review as an intentional policy decision, and documented that ownership policy in `KNOWLEDGE.md`.
74. Removed stale tracked backup files:
deleted `.github/workflows/release.yml.dynamic.old` and `cpp/CMakeLists.txt.dynamic.old` from the nested repo after confirming they were obsolete review findings.
75. Preserved unsafe parser speed while fixing integer width:
kept safe parser values on FLINT string parsing, but restored `parse_matrix_string_unsafe` to direct digit scanning with `long long` fraction construction so `--unsafe` avoids substring allocation/text parsing while no longer casting large values through platform-width `long`.
76. Refreshed review validation notes:
updated `KNOWLEDGE.md` and `REPO_REVIEW_2026-06-02.md` after the latest `./fracessa/test.sh --full` run and release-style no-native build check.
77. Fixed numeric fraction negative-denominator handling:
normalized signed denominators in `linalg::fraction` numeric constructors before calling FLINT's unsigned-denominator setter, and added constructor regression coverage for negative denominators.
78. Fixed top-bit next-set-bit undefined shift:
updated `bs64::find_pos_next_set_bit` to return the sentinel before shifting by 64 when called at position 63, matching the intended iteration contract and fixing the GitHub Linux/macOS test failure seen on the `v0.33` release run.
79. Added bitset first-set-bit call-chain documentation:
created `fracessa/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md` to trace every active production and test call path upward from `bs64::find_pos_first_set_bit`, including direct callers, analyzer entry points, and hot-path summaries.
80. Narrowed first-set-bit call-chain documentation:
updated `fracessa/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md` to exclude test-only call paths and show only the main production code paths through CLI, pybind, analyzer search, candidate detection, stability, MatrixServer, and copositivity.
81. Optimized set-index extraction:
changed `bs64::extract_set_indices` to mask once for `[0, dimension)` and iterate by clearing the lowest set bit (`bits &= bits - 1`) instead of calling `find_pos_first_set_bit`/`find_pos_next_set_bit`; updated the first-set-bit call-chain document and project knowledge to match the new hot path.
82. Clarified first-set-bit usage analysis:
updated `fracessa/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md` so it opens with the current direct production locations where `bs64::find_pos_first_set_bit` is still used and explicitly notes that `extract_set_indices` is no longer a caller.
83. Added correctness regression matrices:
extended `fracessa/python/verification/verification_matrices.json` with IDs 36 and 37 for the original and minimized floating-positive-definiteness false positives, plus ID 38 for the `10^-20`-scaled two-strategy game rejected by the double candidate filter; added their mathematically correct candidate rows to `baseline_candidates.csv` so the cases expose the two open correctness findings through CTest.
84. Expanded verification edge-case coverage:
added IDs 39-44 for common-payoff translation precision loss, dimension one, a negative-payoff mixed ESS, an ESS decided by the second condition after a tie, a fully neutral game, and a zero-payoff mixed ESS; recorded exact candidate baselines, with ID 39 intentionally red in the current default floating path.
85. Consolidated the repository root and source layout:
moved the existing Git metadata from the former nested `fracessa/` directory to the workspace root while preserving history and `origin`; moved active C++ code, tests, local builds, and profiling under `cpp/`, kept Python at top-level `python/`, and placed the README/logo plus project analysis, experiments, papers, and memory files at the repository root.
86. Refactored relocated build and automation paths:
updated root build/test/speed scripts, release workflow, CTest verification data paths, pybind/legacy Python discovery, Callgrind tools, experiment scripts, README, and analysis snapshot tooling for `cpp/build` and top-level `python/`; expanded `.gitignore` to exclude generated builds, results, logs, profiles, caches, local editor state, and experiment source/build copies.
87. Unified CMake Python discovery:
made CMake select the project script's `python3`/`python` command unless explicitly overridden, enabled pybind11's modern FindPython mode, and reused the same interpreter for the extension and Python-driven CTests to prevent mixed Python headers and extension ABI suffixes.
88. Recorded project-local agent workflow policies:
documented normal elevated `ccache` use when sandbox permissions block its cache and made the Ponytail skill mandatory for all FracESSA code modification, analysis, review, design, and code summarization.
89. Normalized project documentation under `aidocs/`:
moved canonical knowledge/history, analyses, repository review, bitset call-chain notes, wrapper documentation, profiling summary, benchmark reports, and mathematical papers under topic folders in `aidocs/`; added `aidocs/INDEX.md`, reduced root `AGENTS.md` to a pointer, retained root `README.md` for GitHub, and updated live documentation/script paths.
90. Hardened CMake portability:
raised the declared minimum to CMake 3.18 to match `Python3 Development.Module`, delegated LTO to CMake's successful IPO capability check instead of forcing `-flto`, and preserved explicit GMP/MPFR library overrides.
91. Updated release action runtimes:
moved checkout, Python setup, cache, artifact upload/download, and release publication to their current Node 24-compatible major versions; validated the workflow with `actionlint` 1.7.12.
92. Revalidated clean native and release-style builds:
rebuilt every target from empty native and `FRACESSA_NATIVE_ARCH=OFF` trees, passed 10/10 core CTests and 23/23 wrapper tests in both, and confirmed that only the known correctness regressions 36-39 fail the 42-case fast matrix gate.

## 2026-07-26
93. Simplified and revalidated agent documentation with Ponytail:
replaced stale point-in-time analyses, generated source dumps, session history, duplicate wrapper contracts, the dated review copy, and old profiling output with a maintained documentation set; consolidated current architecture and workflow facts in `aidocs/KNOWLEDGE.md`, moved all unresolved review findings to `aidocs/reviews/OPEN_ISSUES.md`, merged the wrapper contract into `aidocs/python-wrapper/README.md`, retained the requested compact bitset call chain and dated benchmark reports, and explicitly kept `aidocs/CHANGES.md` as the append-only human-readable project history.
94. Split and revalidated the active repository review:
re-audited the C++ core/build/release paths and the Python wrapper/verification/benchmark paths independently with Ponytail; replaced the single review ledger with `aidocs/reviews/CPP_REVIEW.md` and `aidocs/reviews/PYTHON_REVIEW.md`; preserved every still-valid unresolved finding; added source-backed findings for safe-parser zero-denominator crashes, FLINT deallocation, strict verifier parsing, empty/failing sink output, and ineffective timing suppression; and updated all maintained documentation links. No production code changed in this review task.
95. Double-checked the split reviews independently and with Ponytail:
re-ran the Release build, 10 core/CLI tests, 23 wrapper tests, and all 44 matrix tests; reproduced the parser crash, bitset undefined behavior, queue deadlock, and sink/verifier failures; corrected overstated subprocess, timeout, parser, atomicity, and index-reuse claims; removed speculative packaging and duplicate omnibus-test findings; and added the missed CMake Debug-configuration issue. No production code changed.
96. Added an explicit source-code approval gate:
documented that every `*.cpp`, `*.hpp`, and `*.py` change requires Reinhard's prior approval of the exact files and scope, must remain minimal, requires renewed approval if scope expands, and must be presented afterward for direct review and acceptance.
97. Removed floating positive-definiteness as an ESS certificate:
deleted the double Bee/Cholesky shortcut from `check_stability()` so exact rational positive-definiteness is the only PD decision; verification IDs 36 and 37 now return zero ESS in both default and `--exact` modes, while the complete matrix suite improves from 40/44 to 42/44 with only the separate IDs 38-39 candidate-filter issue remaining.
98. Documented the double-PD counterexample mathematically:
added `aidocs/correctness/DOUBLE_PD_FALSE_POSITIVES.md` with both regression games, their exact singular Bee matrices, observed rounded Cholesky pivots, an explicit `eta = 10^-20` indefinite game, a scaled construction, and the indistinguishable-rounding proof that no double-only tolerance can certify exact positive-definiteness.
99. Added a normal-number double-PD counterexample:
documented the all-integer game `4*A37`, whose payoffs range only from `-18` to `10` but whose exact singular Bee matrix leaves a `2.13e-14` rounded Cholesky pivot; verified that the pre-fix program reports two false ESS even with a hypothetical `10^-20` absolute or relative tolerance, while the corrected exact path returns zero.
100. Constructed a normal-input `10^-20` rounded-pivot example:
scaled matrix 37 by the exact binary factor `4^-9`, yielding payoffs only around `5.7e-6` to `1.7e-5` but a spurious double Cholesky pivot of `2.0328790734103208e-20`; verified two false pre-fix ESS versus zero exact ESS and documented why an order-one diagonal cannot store a nonzero `10^-20` subtraction result in binary64.
101. Strengthened the `10^-20` rounded-pivot example without small-looking payoffs:
added one to every payoff in the scaled game, putting every displayed input between `0.99998283386230469` and `1.0000095367431641` while preserving the exact singular Bee matrix and its spurious `2.0328790734103208e-20` double pivot; verified two false pre-fix ESS versus zero with exact positive-definiteness.
102. Recorded a certified candidate-filter design:
documented that the affine normalization constant is an exact input matrix entry such as `A(0,0)`, not a solved equilibrium payoff; specified exact translation/scale normalization, a rigorous residual/error-bound rejection certificate with exact fallback, the limits of input-level "good matrix" warnings, and the required correctness/performance experiment; expanded the active C++ finding to include the rounded-probability rejection and the exact evidence scope of verification IDs 38-39. No production code changed.
103. Expanded and sourced the certified candidate-filter design:
identified Rump's Theorem 10.2 as the direct source of the residual/error certificate; added the Neumann-series derivation, exact FracESSA rejection inequalities, rigorous rational-to-binary64 enclosure and rounding requirements, interval and round-to-nearest implementation routes, failure semantics, cost model, correctness/performance validation plan, and primary references to Rump, Ogita-Rump-Oishi, Krawczyk, and IEEE 1788.1. No production code changed.
104. Fixed safe zero-denominator parsing and the dimension boundary contract:
made the shared safe string fraction constructor reject a zero denominator before FLINT canonicalization, added fraction-level and CLI crash regressions, confirmed the safe parser accepts dimensions 1 through 63 and rejects 0/64, changed the unchecked parser test to dimension 63 without adding validation, and documented that support masks provide 64 storage bits while complete analyzer search requires `1 <= n < 64`. The Release build, 10/10 core/CLI CTests, 23/23 wrapper tests, and all 40 expected-green fast verification matrices passed; long-run IDs 32/34 and known candidate-filter regressions 38/39 were excluded.
105. Restricted zero-denominator validation to safe parsing:
removed the check and regression from the general fraction string constructor, moved the check into `parse_fraction_token()` in `matrix_parser.cpp`, and replaced the fraction-level test with a safe-parser unit regression. Numeric fraction construction, the unsafe parser, and hot arithmetic remain unchecked and unchanged by design.
106. Closed alternate zero-denominator spellings in the safe parser:
after the final boundary audit reproduced crashes for whitespace forms such as `1/0 ` and `1/ 0`, made `parse_fraction_token()` require an optional denominator sign followed only by decimal digits before FLINT construction and added regressions for signed, zero-padded, trailing-space, and leading-space zero denominators.
107. Split native binding findings into a third review category:
created `aidocs/reviews/PYBIND_REVIEW.md` for the `fracessa_core` boundary; moved pybind error-status reclassification, ESS-count narrowing, and native timing overhead out of the C++ and Python ledgers; kept malformed dimension parsing with the shared C++ parser; and updated review scopes, startup routing, project knowledge, and cross-references without duplicating findings.
108. Corrected FLINT string deallocation:
replaced both plain `free()` calls for `fmpq_get_str(nullptr, ...)` results with FLINT's matching `flint_free()` and removed the obsolete `<cstdlib>` include; the focused fraction test and complete Release build passed, as did all core/CLI and wrapper tests, while only the already-known verification matrices 38-39 remained red.
109. Made safe dimension parsing strict and allocation-free:
replaced partial-prefix `std::stoull` parsing with complete-range `std::from_chars` parsing, so malformed dimensions such as `2x`, alphabetic prefixes, whitespace, and leading plus signs are rejected without constructing a temporary substring or throwing; added a focused parser regression, left the intentionally unchecked parser unchanged, and removed the resolved review finding.
110. Corrected the CLI candidate table header:
removed the trailing semicolon from `candidate::header()` so its 11 fields match candidate rows, and extended the real-executable black-box test to compare every emitted row's field count with the header; existing CLI consumers index the 11 data fields and require no compatibility changes.
111. Defined and accelerated lowest-bit mask extraction:
made `ctz64()` explicitly require a nonzero mask and removed the MSVC zero branch, replaced the `ctz`-plus-shift implementation of `lowest_set_bit_as_bit()` with the documented unsigned `x & -x` expression, retained defined zero-to-zero behavior without a branch, and removed only tests that demanded a nonexistent first-bit position for zero. On ARM64 with GCC 16.1 and `-O3 -march=native`, repeated scalar measurements were approximately 0.512 ns versus 0.551 ns per operation, about 7% faster; the focused bitset, core/CLI, and wrapper suites passed, with only known matrix regressions 38-39 remaining.
112. Made `--fullsupport` avoid eager support enumeration:
constructed the full-support mask directly and delayed `Supports::initialize()` until normal or fallback enumeration is actually needed; added compact integration coverage for both immediate full-support success and failed-full-support fallback. A dimension-24 exact full-support smoke case completed in under 0.01 seconds at approximately 12.3 MiB RSS instead of first generating 16,777,215 masks; all core/CLI and wrapper tests passed, with only known matrix regressions 38-39 remaining.

## 2026-07-27
113. Benchmarked default candidate filtering against full `--exact` execution:
used the persistent-process three-second protocol on verification IDs 1-33 and 35, excluding matrix 34 by request; the summed primary medians increased from 5.129 seconds to 171.323 seconds (`33.41x`), while a second pass measured `32.50x`; all ESS and candidate counts agreed in the detailed pass, and the per-matrix results and single-sample limits are preserved in the dated experiment report and raw CSV/JSON files.
