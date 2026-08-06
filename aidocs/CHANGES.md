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

## 2026-07-28
114. Planned a clean Choice 1 candidate-filter reimplementation:
created `choice-one-candidate-filter` from exact `main` in an isolated worktree and recorded the minimal mathematical, source, comment, and verification plan without merging the experimental branch or changing production source.
115. Double-checked the Choice 1 reimplementation plan with Ponytail:
verified its minimal nine-file C++ scope against `main` and the frozen Choice 1 implementation, made the proof-critical enclosure and bound formulas explicit, documented the dimension and floating-environment contracts, and corrected IDs 45-47 from an open decision to required regression data already committed on `certified-filter`.
116. Added the only useful source excerpt to the Choice 1 plan:
recorded the exact C++17 adjacent-binary64 stepping functions because their signed-zero and bit-order behavior is language-specific; left the remaining kernel as formulas instead of duplicating the full implementation.
117. Closed the final Choice 1 source-scope ambiguity:
confirmed that the legacy double game, bordered-system, Bee, and solver symbols have no production callers outside the candidate filter, then named the six exact `MatrixServer` removals in the plan rather than permitting a broader cleanup.
118. Planned the separate unsafe numerical mode and parser rename:
read the implemented `certified-filter` branch, selected the shorter public names `--unsafe-parser` and `unsafe_parser`, specified `--unsafe` as an opt-in normalized historical double filter with the same constant-cost danger veto and exact fallback, recorded the exact C++/pybind/Python source scope, and made the preserved IDs 45-47 explicit counterexamples rather than claiming the heuristic is certified.
119. Corrected the candidate-filter implementation order:
revised the current plan to contain only temporary unsafe-default and exact numerical modes, reused the existing `exact` routing instead of adding a third mode or Choice 1 scaffolding, kept `--unsafe-parser` independent, deferred the numerical Python `unsafe` argument and strict proof kernel to the later Choice 1 phase, and postponed activation of known unsafe counterexamples 45-47 until Choice 1 becomes the default. No production source changed.
120. Re-audited the temporary unsafe-filter plan with Ponytail:
moved normalization to one conditional pre-enumeration initialization so no availability branch enters the per-support path, required the fused solver to reuse one support extraction, removed only double-Bee state made invalid by conditional initialization, added a missing subnormal-conversion fallback and the unresolved-review update, and separated small-matrix performance acceptance from aggregate timing after recording the reference mode's large per-matrix variation. No production source changed.
121. Corrected the unsafe-filter plan after a second full audit:
required the CLI warning to name both false-negative risk and exact-fallback cost, reused existing payoff and constant-matrix fixtures instead of inventing an unpreserved random-corpus harness, covered both zero and subnormal conversion fallback in one focused MatrixServer test, included the plan and index in implementation-time documentation updates, and treated earlier benchmark results as historical screening data rather than falsely paired measurements. No production source changed.
122. Closed the final unsafe-plan verification ambiguity:
separated the observable `--exact` warning test from the structural source inspection that proves unsafe initialization is never called and double buffers remain unallocated, avoiding pointless allocator instrumentation while making the no-allocation requirement an explicit acceptance condition. No production source changed.
123. Implemented the temporary unsafe-default candidate filter and parser rename:
made no flag and `--unsafe` select the same explicitly warned uncertified numerical mode, renamed unchecked parsing to `--unsafe-parser`/`unsafe_parser`, rejected conflicting `--exact --unsafe`, and kept exact mode free of double initialization. Added one-time exact affine normalization with zero/subnormal/underflow fallback, fused the small bordered double solve into `find_candidate_dbl`, reused fixed support arrays and scratch storage, added the cheap danger veto with exact fallback, and removed the now-unused generic double solver and double-Bee storage. The Release and ASan/UBSan builds each passed 10/10 core/CLI, 24/24 wrapper, and 44/44 matrix checks; default and exact candidate output also matched for all 42 matrices permitted to run exact, excluding IDs 33-34 by instruction. The existing three-second benchmark matched all ESS/candidate counts and changed summed medians by `+0.052%` for IDs 1-12 and `+0.629%` overall versus the saved unsafe run; known reference IDs 45-47 still prove this heuristic is not certified.
124. Removed only C++ helpers with no callers anywhere:
deleted three unused const `data()` accessors and the unused binary fraction subtraction operator; retained every test-only helper for separate review and changed no tests or runtime behavior.
125. Removed C++ functionality maintained only by its tests:
deleted the retired double positive-definiteness implementation, the unused double column accessor, and six unused fraction conveniences; removed their dedicated tests while preserving production coverage by testing equality through `operator==` and normalization fallback through `div_inplace()`.
126. Removed redundant C++ includes:
deleted nine unused include directives from the analyzer and matrix headers, including the obsolete double-PD dependencies, and replaced the two removed bitset dependencies with direct `<cstddef>` and `<utility>` includes required by the matrix containers.
127. Rewrote native algorithm comments for human readers:
documented the support-to-candidate-to-ESS flow, bordered equilibrium equations, extended support, exact-versus-unsafe boundary, equilibrium-based superset pruning, circular-symmetry reduction, Bee construction, exact positive-definiteness, Bomze reduction, and Hadeler copositivity test across the native core; corrected misleading filter, pruning, no-initialization, and bit-string comments without changing executable tokens.
128. Clarified the remaining production C++ boundaries and fixed trailing-comma parsing:
documented the shared parser contracts, FLINT rational ownership and unchecked denominator precondition, minimal double scratch matrix, independent CLI modes, safe-versus-unsafe scanning, and pybind GIL/result boundary; made the safe parser expose and reject a final empty token after a trailing comma, added its focused regression, and left the unchecked parser unchanged. The focused parser test, 10/10 core/CLI tests, 24/24 wrapper tests, and 42/42 fast matrix checks passed.
129. Standardized in-place fraction addition and subtraction notation:
replaced the Bee builder's named `add_inplace()`/`sub_inplace()` calls with `operator+=`/`operator-=`, made those operators call FLINT's alias-safe `fmpq_add`/`fmpq_sub` functions directly, removed the now-unused named helpers, and updated their focused tests. The rebuilt Release CLI was byte-for-byte identical before and after the notation change; 10/10 core/CLI tests, 24/24 wrapper tests, and 42/42 fast matrix checks passed.

## 2026-07-29
130. Re-audited the complete active C++ runtime with Ponytail:
updated every stale reference in `CPP_REVIEW.md`, retained the known unsafe-filter false-negative finding, and added confirmed findings for signed overflow in large `--unsafe-parser` values, eager `2^n` support-frontier storage, non-monotonic CLI timing, and two uncovered exact proof branches. The audit reproduced a safe-versus-unsafe exact ESS-count difference for `2#9223372036854775808,0,1`, confirmed the overflow with UBSan, measured 8,388,607 stored masks and about 69 MiB RSS at `n=23`, and passed the current Release build, 10/10 core/CLI tests, 24/24 wrapper tests, and 44/44 full matrix checks. No source code changed.
131. Unified and accelerated matrix parsing:
removed the unchecked parser implementation, `--unsafe-parser` CLI option, pybind `unsafe_parser` argument, and Python `RunConfig.unsafe_parser` field. The single validating parser now constructs numerator and denominator components of at most 18 digits directly and sends larger exact values to FLINT text conversion without signed overflow. It reserves symmetric storage once and validates value syntax, counts, and zero denominators while constructing the fractions; the redundant second-separator scan was removed because the value scanner already rejects it. Production medians were approximately 134 ns for 2-by-2 integers, 278 ns for 4-by-4 integers, and 419 ns for 4-by-4 fractions; a compact-first reserve experiment was rejected because its 2% circular gain caused 4-23% symmetric regressions. The Release build, 10/10 core/CLI tests, 23/23 wrapper tests, all 44 matrix checks, and focused ASan/UBSan parser tests passed.
132. Benchmarked and deferred reusable exact-solver output:
compared the current per-solve output allocation with a pre-sized scratch vector using fixed exact systems and grouped support sweeps. Reuse saved roughly 32-40 ns on the smallest successful solves but only 0.9-1.8% across grouped sweeps of 255-4,095 supports; all 5,373 compared exact solutions matched. Removed the finding from `CPP_REVIEW.md` because the measured aggregate gain does not justify a three-file production change now. No production code changed.
133. Explored on-demand support generation and superset block jumps:
benchmarked seven isolated alternatives to eager support storage and repeated `remove_if`, including fixed-cardinality Gosper generation, direct superset marking, byte-indexed block jumps, split generation, DFS hypergraph pruning, and an adaptive combination. The final adaptive executable was byte-identical to production on all 44 official and 35 generated non-circular matrices; paired strong-pruning cases improved by 31-74%, matrix-34 peak RSS fell by 48.7%, and weak cases ranged from small wins to roughly 4-5% regressions. Rejected a dimension-only DFS gate and the 235-line adaptive design as unnecessarily complex; recorded the 157-line byte-indexed Gosper variant as the simpler production candidate. No production source changed.
134. Found the missing prime-circular byte-jump regression:
benchmarked official dimension-19 IDs 28 and 29 after noticing that the representative table omitted them. Both are circular-symmetric at prime dimension, but byte-indexed Gosper was 52.07% slower on ID 28 and 2.30% slower on ID 29. ID 28's 1,843 rotated size-7/8 candidates trigger approximately 6.73 million direct superset-table updates even though production stores only one canonical mask per 19-member orbit. Corrected the experiment and review conclusion: the byte design is not production-ready until direct marking becomes circular-aware. No production source changed.
135. Completed the byte-indexed Gosper comparison for every verification matrix:
ran paired three-second production/byte benchmarks for all 44 fixtures, split into 17 ordinary full matrices and 27 circular-symmetric matrices, and confirmed byte-identical candidate output afterward. The composite circular dimensions 18, 20, 21, 22, and 24 produced the five major wins of 30-73%, while prime circular IDs 22, 26, 28, and 33 regressed by 2.8-52.4%; ID 29 at the same prime dimension 19 regressed only 2.2%, proving that primality controls orbit reduction but not the matrix-specific pruning workload. No production source changed.
136. Added larger non-circular verification matrices:
added active IDs 48-54 at dimensions 15-24 using exact Hilbert, Sylvester Hadamard, symmetric Paley conference, MINIJ, Fiedler, deterministic mixed-integer, and deterministic positive-integer families; preserved IDs 45-47 for the certified-filter regressions; and added all 133 production candidate rows to the baseline. Paired three-second production/byte-Gosper measurements ranged from a 65.70% MINIJ win to a 6.44% mixed-random regression, with byte-identical candidate output. The Release build, 10/10 core/CLI tests, 23/23 wrapper tests, and all 51 active matrix checks passed. No production source changed.
137. Added an isolated general binary-bracelet generator:
implemented fixed-weight necklace generation followed by reflection reduction for dimensions 1 through 63, and generated complete representative listings for dimensions 8 and 10. Orbit-size checks cover all 255 and 1,023 nonempty raw supports respectively, and an independent exhaustive comparison matched every representative: 29 bracelets for dimension 8 and 77 for dimension 10. No production source changed.
138. Implemented and measured orbit-aware circular bracelet pruning as an isolated experiment:
generated one fixed-weight representative per rotation/reflection orbit, registered every transformed exact-candidate support as a future subset-pruning rule, and reconstructed all equivalent mathematical candidates from one solve. All 27 circular fixtures and 104 generated circular games matched production after excluding only enumeration metadata, all 24 non-circular fixtures remained byte-identical, and sanitizers passed. Paired three-second timings improved 24/27 circular fixtures with a `4.28x` geometric-mean speedup; matrix 34 fell from 22.132 seconds to 0.650 seconds and peak RSS from 99,520 KiB to 13,008 KiB. Candidate order, IDs, and shift references remain an explicit compatibility decision before production use. No production source changed.
139. Compared direct fixed-content bracelet generation with the existing FKM-plus-reflection method:
adapted the paper's optimized `BraceletFC` recursion in a standalone C++ experiment while retaining the current generator unchanged; all support-size orbit sets matched without duplicates through dimension 24 under Release and ASan/UBSan builds. Generator-only three-second measurements found direct generation `17.34%` slower by geometric mean at dimensions 8-14 but `12.19%` faster at dimensions 19-24, reaching a `22.31%` win at dimension 24. The result does not justify an unconditional replacement for small-matrix workloads; no production source changed.
140. Re-ran support-frontier alternatives only on non-circular matrices:
compared production, byte-indexed Gosper, and pure fixed-cardinality DFS on all 24 active non-circular fixtures, including dimensions 15-24; all complete candidate outputs were byte-identical. Pure DFS improved large IDs 49, 51, 53, and 54 by `40.68%`, `67.50%`, `33.19%`, and `14.52%`, but regressed ID 52 by `2.67%`; byte-Gosper's remaining wins were confined to tiny absolute runtimes. Saved all raw and consolidated three-second results under the support-frontier experiment; no production source changed.
141. Added a non-circular unique-full-support ESS control matrix:
constructed the 20-by-20 game `diag(-1,-2,...,-20)`, proved that every proper support is invaded by an omitted zero-payoff strategy and that the unique full-support equilibrium is ESS by negative definiteness, and confirmed one candidate with support mask `2^20-1`, support size 20, and exact `T_pd_frc` stability. Saved the matrix, output, and one-shot timing under the support-frontier experiment; no production source changed.
142. Replaced the sparse full-support ESS control with a dense one:
used `A = 11^T - diag(2,3,...,21)`, giving nonzero diagonal entries `-1` through `-20` and off-diagonal entries all equal to one. Every proper-support equilibrium is invaded because an omitted strategy earns one instead of `1-c_S`, while the unique full-support equilibrium is ESS because the matrix is strictly negative on the simplex tangent space. Normal enumeration and exact full-support execution agreed on one ESS candidate with all 20 strategies; no production source changed.
143. Varied every off-diagonal entry in the dense full-support ESS control:
replaced the constant off-diagonal construction with the weighted complete-graph game `A_ij=i+j` and `A_ii=-(18i+211)` for one-based indices. Its 37 distinct off-diagonal values range from 3 through 39, it has no zero entries, every row sums to `-1`, and `A=-I-L` is strictly negative definite. Normal enumeration and exact full-support execution again agreed on exactly one candidate: the uniform 20-strategy ESS; no production source changed.
144. Promoted the dense full-support ESS control to verification matrix 55:
added its compact and full representations plus the single uniform candidate baseline, registered it automatically as a fast CTest, and passed the focused correctness check. Two independent pinned three-second medians averaged `427.008 ms` for production, `426.733 ms` for byte-Gosper (`-0.064%`), and `426.025 ms` for pure DFS (`-0.230%`); candidate outputs were byte-identical, so neither alternative produced a meaningful gain on this no-pruning case. No production source changed.
145. Benchmarked a table-free Gosper candidate-list scan:
created one isolated experimental support header that retains only the current cardinality layer and exact-candidate masks, orders candidates by descending lowest bit, and jumps over matching forbidden Gosper blocks without a `2^n` bitmap or byte table. Complete candidate output was byte-identical to production on all 25 non-circular matrices. Pinned three-second measurements won 18/25 rows; large IDs 49, 51, 53, and 54 improved by `38.22%`, `66.60%`, `31.81%`, and `12.47%`, while the four other dimensions 15-24 changed by `+0.19%` to `+0.97%`. The all-row geometric mean improved `9.01%` and summed medians improved `2.92%`; no production source changed.
146. Tested a streaming support-generator interface:
replaced the experimental one-cardinality vector with `start_cardinality()` and `get_next_support()`, retaining only one Gosper mask plus exact-candidate pruning state. Complete candidate output was byte-identical to production on all 52 active matrices. Streaming preserved the candidate-scan gains but changed summed non-circular medians by `+0.24%` versus retaining one layer and was `1.02%` slower by geometric mean on dimensions 15-24. On a dimension-24 no-pruning control, peak RSS fell from 34,576 KiB to 13,632 KiB (`-60.57%`) while elapsed time changed from 9.37 to 9.39 seconds; no production source changed.
147. Created a staged SQLite test-matrix database:
defined a strict three-table data model plus metadata under `testdata/`, imported all 52 active matrices and their candidate baselines exactly, and added IDs 56-66 as non-circular complete-multipartite many-ESS benchmark matrices. The database contains 63 matrices, 65,962 candidates, and 63,369 ESS rows; each matrix stores counts, support-size structures, qualitative tags, and required origin/purpose text instead of a name. Integrity, foreign keys, source hashes, exact source round-trips, candidate IDs, support masks, vector supports, structures, and constructed-family counts all passed. Existing JSON/CSV files remain the active Python and CTest inputs, and benchmark-run storage remains deferred; no production source changed.
148. Simplified test-matrix tags:
moved all 114 qualitative tag assignments into a checked JSON array on each `matrices` row and removed the separate `matrix_tags` table. Every association was preserved in sorted order; no production source changed.
149. Removed test-database metadata bookkeeping:
dropped the `dataset_metadata` table and its source hashes, creation date, and schema-version fields. Matrix provenance remains in each row's required `origin`, while source history remains in Git and project documentation; the database now contains only matrices and candidates. No production source changed.
150. Recovered provenance for the historical test matrices:
matched IDs 1-35 against archived REF/EFR regression suites, documentation examples, matrix exchanges, manuscript material, and Werner Schachinger's extended comparison corpus; replaced generic origins with the strongest factual group or individual history supported by those records. No production source changed.
151. Refined two historical test-matrix origins with external sources:
recorded the 2013 note identifying ID 4's original nonsymmetric game with Bomze's 1986 classification paper and matched ID 9 exactly to Example 2.1 of Bomze's 1992 ESS paper. In-house REF/EFR matrices retain only their project provenance rather than references to Reinhard Ullrich's own publications. No scientific data or production source changed.
152. Excluded generated experiment check executables from version control:
kept the reproducible check source and benchmark results while ignoring local architecture-specific test binaries under experiment `checks/` directories.

## 2026-07-30
153. Recorded the deferred DFS/FKM support-generation architecture:
documented the current eager frontier, fixed-cardinality DFS branch-pruning proof, Gosper tradeoffs, FKM-style necklace and bracelet generation, circular orbit reconstruction, the failure of canonical-only subset tests with explicit counterexamples, correct orbit-storage alternatives, all retained timing and memory evidence, literature boundary, staged implementation experiment, and correctness/performance acceptance gates. Clarified that enumeration metadata may change only after order-insensitive mathematical equivalence is proven. No production source changed.
154. Replaced the eager support frontier with production support generators:
added callback-driven fixed-cardinality DFS for non-circular matrices and fixed-content FKM bracelet generation for circular matrices, fused lowest-bit exact-candidate pruning into both recursions, reconstructed complete circular rotation/reflection candidate orbits, and split candidate analysis from orbit recording. Each generator now owns the complete cardinality sweep and passes both the mask and support size to the analyzer callback; the analyzer selects the concrete generator with two direct branches and no generic wrapper lambda. An independent order-insensitive comparison matched every mathematical candidate and ESS result across all 52 active matrices; regenerated candidate IDs and shift references for 23 circular fixtures in the CSV and SQLite baselines; all 62 CTests passed in both Release and ASan/UBSan builds. Renamed the maintained design to `plans/SUPPORT_GENERATORS.md`; the compact bit-parallel orbit representation remains an optional benchmark rather than a third production class.
155. Compressed circular candidate output to one bracelet representative:
removed rotated/reflected candidate-object reconstruction and `shift_reference`, stored only the smallest support in each circular dihedral orbit, and added its distinct-orbit `multiplier`; non-circular rows use null. Circular generators still queue every compact orbit mask required for correct pruning, and weighted ESS totals are unchanged. Propagated the nullable field through CLI, pybind, both Python wrappers, CSV/JSON/Arrow/Parquet sinks, baseline regeneration, and CTest verification. The active CSV fell from 38,044 expanded rows to 1,196 representatives, while the SQLite snapshot fell from 65,962 rows to 29,114 representatives and retained weighted totals of 65,962 candidates and 63,369 ESS. All 62 CTests passed in Release and ASan/UBSan builds, all 25 wrapper tests passed, and database integrity, foreign keys, multiplier bounds, canonical supports, weighted counts, and support-size structures matched.
156. Added the compact circular-generator comparison implementation and design handover:
kept production on expanded dihedral forbidden masks, added test-only `CircularSupportGeneratorV2` with one forbidden representative and bit-parallel rotation/reflection alignment, derived its multiplier directly from the FKM period and reflection relationship, and moved arbitrary-width left rotation into the shared bitset operations. Recorded the callback, analyzer/finalization, inheritance, output, validation, and deferred V2 decisions in `architecture/SUPPORT_GENERATOR_HANDOVER.md`; production wiring remains unchanged.
157. Consolidated local REF/EFR history:
moved the four preserved REF/EFR trees under the ignored `zzz_legacy/` collection and renamed them by language and historical period without changing their contents; retained nested Git metadata, including the existing dirty legacy REF worktree. A projects-wide filename, source-signature, and ZIP-content scan found no additional REF/EFR implementation roots. The collection remains local-only and is not active project source.
158. Reduced redundant legacy storage without removing unique source:
permanently removed the two extracted Boost 1.71 trees and seven byte-identical duplicate Boost 1.62 trees after verifying both retained Boost 1.71 archives and the complete canonical Boost 1.62 tree hash; moved five byte-identical nested 2018 project copies and 24 generated `build*`/`obj` directories to the desktop Trash for recovery. Retained all Git metadata, unique source variants, binary output, dependency packages, both Boost 1.71 distributions, and one Boost 1.62 tree. The active collection fell from 3.19 GB to 576 MB allocated, while approximately 119 MB remains recoverable in Trash.
159. Removed nested Git control from the legacy collection:
moved the two internal `.git` directories to the desktop Trash and confirmed that no internal `.git`, `.gitignore`, `.gitattributes`, `.gitmodules`, `.gitkeep`, or `.github` path remains. The four legacy folders now occupy 538,533,888 allocated bytes (about 514 MiB) together; the main FracESSA repository continues to ignore the local-only `zzz_legacy/` collection.
160. Consolidated the Werner 2016 C++ archive:
renamed the retained 2019-imported copy to `REF_CPP_2016_Version_Werner`, preserved the only remaining `CMakeLists.txt` and Boost 1.62 tree inside it, and moved the now-redundant `REF_history_mixed_2008-2019/REF_version_werner` copy to the desktop Trash.
161. Clarified the local 2019 REF archive name:
renamed `REF_CPP_2019_modified_2025` to `REF_CPP_2019` because its core implementation is the substantial 2019 local rewrite; the December 2025 differences are limited to optional timing output and compiler settings.
162. Created a direct-move EFR timeline:
created `zzz_legacy/EFR/` and moved the four selected complete C# project trees into `EFR_2016-04`, `EFR_2016-09`, `EFR_2018-03`, and `EFR_2019-08`. Directory inode checks confirmed that each tree was moved rather than copied; the standalone `NewRational.2.1` experiments, duplicate EFR trees, and skipped intermediate September 2016 port were left untouched.
163. Removed the skipped intermediate EFR snapshot from the active archive:
after confirming identical source and project content, moved both the primary intermediate September 2016 project and its nested `REF_von_2018` duplicate to the desktop Trash. The four selected `zzz_legacy/EFR/` timeline projects were unchanged.
164. Removed duplicate April and September 2016 EFR snapshots:
after verifying that the April tree was byte-identical and that all September program and project files were identical, moved both duplicate trees to the desktop Trash. The sole September difference was MonoDevelop active-document and cursor state in `EFR.userprefs`; the four selected `zzz_legacy/EFR/` timeline projects were unchanged.
165. Consolidated the NewRational archive:
verified that all three `NewRational.2.1` trees were byte-identical, moved one directly to `zzz_legacy/EFR/NewRational`, and moved the two redundant copies to the desktop Trash. The versioned ZIP remains inside the retained project.
166. Consolidated the legacy EFR R test driver:
retained the newer October 7, 2016 `efr_tester.R` at `zzz_legacy/REF_history_mixed_2008-2019/zzz_old/efr_tester.R`, then moved `EFR_CSharp_2008-2016` and both remaining `_EFR` containers to the desktop Trash. Their older or duplicate tester scripts, identical `.RData` workspaces, and trivial `.Rhistory` files were intentionally removed from the active archive.
167. Replaced the mixed REF archive with date-named milestones:
created source-only `REF_2016-10-06` and `REF_2016-11-16` milestones, renamed the retained Werner and 2019 trees to `REF_2016-11-20` and `REF_2019-09-20`, consolidated one copy of each research script under `REF_R`, and retained the actual newest EFR tester from October 21, 2016. Moved both Boost 1.71 distributions into `REF_2019-09-20/dependencies/`, verified the retained source and duplicate relationships, and moved the remaining 63 MiB `REF_history_mixed_2008-2019` tree to the desktop Trash. The discarded remainder comprised superseded intermediate source snapshots, exact duplicates, binaries, objects, profiles, logs, IDE state, R workspace/history, the unrelated `EFCsimple` tutorial, and the duplicated travel document.
168. Marked the Werner REF milestone explicitly:
renamed the local legacy folder `REF_2016-11-20` to `REF_2016-11-20-Werner` and updated its project documentation references.
169. Flattened the Werner REF milestone:
moved the contents of its redundant inner `REF/` directory directly into `REF_2016-11-20-Werner/` and updated the documented Boost path.
170. Removed the redundant Boost 1.71 ZIP distribution:
moved `REF_2019-09-20/dependencies/boost_1_71_0.zip` to the desktop Trash and retained the smaller verified `boost_1_71_0.tar.gz`; the separate extracted Boost 1.62 tree remains available for the November 2016 versions.
171. Published the curated legacy collection:
removed the repository-wide `zzz_legacy/` exclusion, added the complete retained REF/EFR history to `main`, and tracked the 118.5 MB Boost 1.71 archive with Git LFS because it exceeds GitHub's normal 100 MB file limit. Generated builds, caches, and other existing ignored output remain local.
172. Removed Python result-row wrappers:
deleted `SummaryRow`, `CandidateRow`, `MatrixResult`, `RunStats`, and `result_to_dict()`; made sequential, multiprocessing, and sink paths use one flat result dictionary with dictionary candidate rows; updated wrapper tests and documentation while preserving CSV fields and native status/result values.
173. Replaced the public multiprocessing queue runner:
removed `MPQueueRunner` and kept one private bounded-streaming helper behind `run_jobs_mp()`; made job serialization fail synchronously, detected worker exits while waiting, cancelled workers when consumers stop early, used a monotonic shutdown deadline, and added focused regression coverage plus native and spawn validation.
174. Made SQLite the canonical test-data store and removed the obsolete active pipelines:
retained the compressed `testdata/fracessa_testdata.sqlite3` snapshot with its strict schema and documentation; removed the JSON/CSV baselines, archived executables, baseline generator, matrix verifier/selector, legacy subprocess benchmark, wrapper, launcher, and JSON-fed Callgrind scripts plus their Python, CMake, test, release, and documentation integration while retaining the dated historical experiment archives. No replacement verification or benchmark runner was added.
175. Made multiprocessing result serialization fail visibly:
serialized complete worker results before handing them to the queue feeder and decoded them in the parent, so an unpicklable result now exits its worker and reaches the existing liveness error instead of being silently dropped; made the input-serialization test independent of interpreter-specific exception wording and added a return-path regression. All 25 wrapper tests pass on Python 3.14, while Python 3.12 passes its 23 available tests and skips only the two locally unavailable native-extension tests.
176. Prevented file sinks from overwriting prior runs:
routed CSV, JSON, Arrow, and Parquet output pairs through atomic exclusive creation, so a reused run ID raises `FileExistsError` while preserving existing files; added regression coverage for complete and partial existing output, retained readable Arrow/Parquet round trips, and removed the resolved output-ID truncation finding. All 26 wrapper tests pass on Python 3.14, while Python 3.12 passes its 24 available tests and skips only the two locally unavailable native-extension tests.
177. Completed result-sink correctness and batching fixes:
made CSV, JSON, Arrow, and Parquet create exclusive summary/candidate/metadata triplets; streamed non-null per-matrix metadata into format-specific JSON sidecars; gave every empty dataset a readable stable schema; made stream, multi-sink, and file finalization attempt all closes before propagating the first error; and buffered Arrow/Parquet rows in bounded 1,024-row windows, reducing a 1,100-result regression from 1,100 batches or row groups to two. All 31 wrapper tests pass on Python 3.14, while Python 3.12 passes its 28 available tests and skips only two unavailable native-extension tests and one unavailable-Pyarrow test. Parallel analyzer logging remains separate and unchanged.
178. Rejected native logging in multiprocessing:
made `run_jobs_mp()` raise `ValueError` for `RunConfig.enable_logging=True` before constructing its queue runner, preventing independent native workers from writing and rotating the shared `log/fracessa.log`; sequential logging remains unchanged. The focused multiprocessing tests pass on Python 3.12 and 3.14, all 32 wrapper tests pass on Python 3.14, and Python 3.12 passes its 29 available tests with only two native-extension and one Pyarrow skip.
179. Simplified multiprocessing to completion-order shared queues:
removed `MPConfig.ordered_results` and its sequence/pending-result machinery; renamed the misleading `chunk_size` setting to `prefetch_per_worker` while retaining its default of 128; and kept one matrix per shared job-queue item plus one shared result queue, so every idle worker can claim the next matrix without IPC batching. All 10 native tests and 32 wrapper tests pass on Python 3.14, while Python 3.12 passes its 29 available wrapper tests with only two native-extension and one Pyarrow skip.
180. Standardized the Python toolchain on Python 3.14:
changed release CI from Python 3.12 to 3.14, required Python 3.14 or newer during CMake configuration, and upgraded pybind11 from v2.13.6 to v3.0.4 for official Python 3.14 support. A clean Python 3.14.6 configure and rebuild produced the CPython 3.14 extension; all 10 native tests, all 32 wrapper tests, and a 32-job spawn-mode multiprocessing check pass.
181. Simplified Pybind invalid-input reporting:
deleted the second parser and the detailed dimension/value-count status codes; every safe-parser rejection now returns `PARSE_ERROR` with `Invalid matrix string`. One native regression covers malformed structure, invalid dimension, invalid value count, and an invalid rational; all 10 native tests and all 33 wrapper tests pass on Python 3.14.
182. Preserved the analyzer's ESS-count width through Pybind:
changed `NativeResult::ess_count` from 32-bit `int` to the analyzer's `size_t` and removed the narrowing cast, so Win64 and other 64-bit builds pass the full count directly to Python. All 10 native tests and all 33 wrapper tests pass.
183. Returned shared parser diagnostics through Pybind:
changed the safe parser from a boolean, `stderr`-writing API to a throwing API using `std::invalid_argument`; made the CLI catch and print the diagnostic while Pybind catches it into the existing `PARSE_ERROR` dictionary; removed the parser's internal catch-and-print path and added native coverage for detailed returned messages. All 10 native/CLI tests and all 33 wrapper tests pass on Python 3.14.
184. Standardized elapsed timing on monotonic nanoseconds:
changed both the CLI and Pybind analyzer timers to `std::chrono::steady_clock`; renamed the complete native/Python/sink schema from `elapsed_us` to `elapsed_ns`; made CLI `--timing` print integer nanoseconds; and removed `RunConfig.include_timing` so the native duration is always returned. All 10 native/CLI tests and all 32 wrapper tests pass on Python 3.14.
185. Completed the native Pybind boundary cleanup:
widened matrix IDs from 32-bit `int` to signed `std::int64_t` through the CLI, analyzer, and Pybind boundary; added maximum-value CLI and wrapper regressions; added an exact value-and-Python-type contract for all 11 converted native candidate fields, including nullable `multiplier`; and removed the unused `pybind11/stl.h` include. All 10 native/CLI tests and all 33 wrapper tests pass on Python 3.14.
186. Made generic JSON schema failures explicit:
required the configured matrix key in top-level objects, required its selected rows to be a list, and required each row to be a JSON object; added one table-driven regression for a misspelled key, non-list container, and scalar row. Valid empty lists remain supported. All 10 native/CLI tests and all 34 wrapper tests pass on Python 3.14.
187. Reduced output sinks to disk formats only:
deleted `StreamSink`, `MultiSink`, and `ArrowSink`; removed their public exports, factory branches, implementation, and tests; and retained only CSV, JSON, and Parquet with their exclusive summary/candidate/metadata output triplets. The factory now rejects `stream` and `arrow`. All 10 native/CLI tests and all 32 wrapper tests pass on Python 3.14.
188. Made sink initialization transactional:
converted the shared exclusive-output opener into a rollback context used by CSV, JSON, and Parquet; partial initialization now closes raw files and any constructed Parquet writers, removes only paths created by that attempt, and re-raises the original error so the same run ID can be retried. A forced second-JSON-writer regression verifies cleanup, error propagation, and successful retry; a direct partial-Parquet-writer probe also passed. All 10 native/CLI tests and all 33 wrapper tests pass on Python 3.14.
189. Required Parquet coverage in release CI:
installed PyArrow after Python setup on every release matrix operating system so the wrapper suite cannot silently skip retained Parquet tests. The workflow parses as YAML, all 10 focused sink tests pass with PyArrow present, and `git diff --check` passes.
190. Removed the unused job-iterator pass-through:
deleted the unexported, uncalled `iter_jobs()` generator and the collection imports used only by it. The focused JSON-loader tests and the complete wrapper suite pass, and `git diff --check` passes.
191. Extended sink rollback through the complete run:
retained each exclusive output triplet's existing `ExitStack` until successful finalization; CSV, JSON, and Parquet now abort and delete their three new paths after initialization, computation, write, or close failures; and sequential and multiprocessing sink consumers share one cleanup path that closes active result iterators without allowing cleanup errors to replace the original exception. Regressions cover direct serialization failure in every format, pipeline computation failure in every format, finalization failure, same-ID retry, and competing write/close errors. All 10 native/CLI tests and all 37 wrapper tests pass on Python 3.14.
192. Closed the remaining Python/Pybind correctness findings:
validated `MatrixJob` IDs, matrix strings, metadata, and JSON dimensions without lossy coercion; made the multiprocessing error fallback safe for invalid mutated IDs; encoded every non-finite JSON float as `null` under strict JSON; required the native binding in integration tests; ran the three-platform build/test workflow for ordinary pushes and pull requests while limiting packaging and publication to release tags; and serialized logging-enabled threaded Pybind calls with one native mutex while leaving non-logging execution concurrent. All 10 native/CLI tests and all 43 wrapper tests pass on Python 3.14.
193. Reconciled the Python/Pybind work with the production support-generator redesign:
merged the current `main` history, retained its compressed SQLite snapshot and one nullable `multiplier` candidate contract, ported the flat dictionaries and transactional CSV/JSON/Parquet sinks without a compatibility layer, and removed the obsolete deferred generator plan. The shared parser still returns detailed exceptions through CLI/Pybind, matrix IDs remain signed 64-bit, and analyzer timing remains monotonic nanoseconds; the public README now documents the actual unit. The Release build, all 10 C++/CLI tests, all 45 Python tests with native and PyArrow coverage, a 64-job spawn run with early close and serialization failure, SQLite integrity/foreign keys, and weighted totals of 65,962 candidates and 63,369 ESS all passed.
194. Simplified the maintained Python execution interface:
renamed the public input object to `Matrix`; consolidated sequential execution into `run` and process-based execution into `run_multiprocessing`; made both accept one matrix or an iterable and an optional sink; retained `compute_matrix` as the sole public low-level native adapter; renamed the JSON loader around matrices; and removed the redundant native status-map export. `run_multiprocessing` now mirrors the `run` parameter order and adds only a final optional `MPConfig`, whose default uses the CPUs available to the Python process. All 47 wrapper tests pass on Python 3.14, including the native default-spawn integration path.
195. Documented the complete production Python wrapper:
added standard docstrings to every production module, class, function, and method, including private multiprocessing and transactional sink lifecycle helpers; documented inputs, return modes, failures, cleanup behavior, and configuration fields; and added `pydoc`, `pdoc`, and Sphinx `autodoc` generation guidance without adding a runtime dependency.
196. Renamed the maintained Python surface to PyFracESSA:
moved the public package and tests from `python/wrapper_v1/` to `python/fracessa/`, changed imports to `fracessa`, moved its public documentation to `aidocs/pyfracessa/`, and updated live scripts, CI, tests, mock targets, architecture references, and generated-documentation commands. The private Pybind extension remains `fracessa_core`; no compatibility package preserves the obsolete `wrapper_v1` name.
197. Added sequential database-backed timing:
added a single-build runner that reads the canonical SQLite matrices, pins execution to one selected Linux CPU, records current Pybind or legacy CLI analyzer timings as nanoseconds, retains observed ESS counts, and groups separately invoked builds by session. A pilot at or above the configurable one-second target is the only iteration; faster cases use the pilot as warmup and average an adaptively sized measured batch, with a measured-duration correction when the initial estimate is short. The same database stores build labels, moving source references, immutable revisions, binary hashes, backend, mode, CPU, comments, target and measured wall durations, iteration counts, and average native durations; a compact report prints correctness and optional ratios to a named baseline. Rejected fixed-five samples were deleted. The current Pybind build then completed all 63 matrices in unsafe and exact modes on CPU 2: all 126 results matched, unsafe accumulated 6,367,755 iterations over 107.308 measured seconds, and exact accumulated 6,741,645 iterations over 1,038.390 measured seconds. All 50 wrapper tests, database integrity, foreign keys, completeness, and the real adaptive Pybind/report paths pass.
198. Completed circular/non-circular matrix coverage for dimensions 2-25:
added deterministic random-integer IDs 67-79 for the 13 missing dimension/circularity combinations. Exact and unsafe Pybind runs agreed on every accepted candidate ID, support, rational vector, ESS classification, stability reason, and exact payoff before the rows were inserted. The canonical database now contains 76 matrices, 29,160 candidate representatives, 66,209 weighted candidates, and 63,434 weighted ESS. Extended the existing pinned-CPU current-build timing session to all 76 matrices and both modes: all 152 observations match, with 7,899,684 unsafe iterations over 120.423 measured seconds and 8,455,777 exact iterations over 1,053.185 measured seconds. Integrity, foreign keys, dimension coverage, weighted structures, and timing completeness pass; no production source changed.
199. Added the published ESS-growth lower bound to correctness reports:
extended every timing report row with matrix circularity, dimension, and `gamma_lower_bound = expected_ess ** (1 / dimension)`, matching the n-th-root lower bound used by Bomze, Schachinger, and Ullrich. The value is derived from the canonical expected ESS count so an unsafe mismatch cannot distort it, and the SQLite schema remains unchanged. A focused dimension-3/eight-ESS regression prints `2.000000`; all 51 wrapper tests, the real 152-row report, and `git diff --check` pass.
200. Aligned the canonical test data with both published Bomze-Schachinger-Ullrich result tables:
audited every literal and constructed matrix in Tables 1 and 2 against exact database values rather than ESS counts alone; replaced the same-count dimension-12 and dimension-17 alternatives at IDs 18 and 26 with the exact Table 1 vectors; added the missing Table 2 base matrices as IDs 80-81 and its nine constructed matrices as IDs 82-90; and removed redundant same-property alternatives at IDs 12 and 21. Unsafe and exact Pybind runs agreed candidate-for-candidate for all 13 changed or added published matrices and reproduced every published ESS total. The database now contains 85 matrices, 49,155 representatives, 86,150 weighted candidates, and 83,375 weighted ESS. Its pinned-CPU current-build session contains all 170 unsafe/exact samples with no mismatch: unsafe accumulated 7,854,428 iterations over 136.921 measured seconds, and exact accumulated 8,350,554 iterations over 1,339.616 measured seconds. SQLite integrity, foreign keys, aggregate counts, table coverage, and timing completeness pass; no production source changed.
201. Implemented the bounded-error Choice 1 candidate filter from current `main`:
made the rigorously one-sided rejection proof the no-flag default, retained the existing heuristic behind explicit `--unsafe`, and kept `--exact` as a complete numerical-filter bypass. Added a dedicated bounded-filter object that owns its lazy affine normalization, MPFR rational enclosures, retained-LU solution-error bounds, and scratch state independently of `MatrixServer`; added strict support-probability and outside-gain rejection checks, exact fallback for every inconclusive support, and an optional exact-rejection oracle; isolated the proof object from fast-math, contraction, and IPO/LTO. Centralized compiler and runtime availability checks so unsupported default mode stops before enumeration and requires an explicit exact or unsafe choice. Propagated the appended unsafe selector through the CLI, pybind, and `RunConfig` without shifting existing positional APIs, and activated exact regression IDs 45-47. Release passed 11/11 core/CLI tests, 26/26 wrapper tests, and 55/55 matrices; the oracle passed all 53 permitted matrices with IDs 33-34 excluded, and ASan/UBSan passed the core/CLI suite plus IDs 45-47. The historical pinned three-second set measured 2,108.563 ms in summed medians, 10.03x faster than the saved pre-generator Choice 1 run and 81.25x faster than saved full exact, primarily because current `main` already contains the new support generators.
202. Clarified ownership of the unsafe normalized matrix:
renamed `MatrixServer::initialize_unsafe_filter()` to `initialize_game_matrix_dbl()` so the server only prepares, stores, and returns matrix data while `fracessa` alone decides whether failed preparation selects exact fallback. The focused MatrixServer conversion regression was retained under the matrix-specific name. Release passed 66/66 CTests and 26/26 wrapper tests; the exact-rejection oracle passed 64/64 permitted tests with verification IDs 33-34 excluded, and ASan/UBSan passed the 11 core/CLI tests plus IDs 45-47.
203. Renamed the rigorous default rejection component to candidate-rejector-double:
renamed its files and build targets with `double` written out and its class, member, and method to `candidate_rejector_dbl`; updated the focused test target, oracle macro, documentation, and links consistently. The rejector method now returns true exactly when rejection is proven, while the separate unsafe `find_candidate_dbl()` retains its existing candidate-possible result. The mathematical proof and routing behavior are unchanged. Release passed 66/66 CTests and 26/26 wrapper tests; the exact-rejection oracle passed 64/64 permitted tests with verification IDs 33-34 excluded, and ASan/UBSan passed the 11 core/CLI tests plus IDs 45-47.
204. Split candidate solving into three concrete procedures:
replaced `MatrixServer`, the mixed `findeq.cpp`, and candidate-rejector naming with `find_candidate_verified`, `find_candidate_unsafe`, and `find_candidate_exact`, each implemented as a state-owning class with matching HPP/CPP files. `fracessa` now owns the exact game once; all three procedures reference it and retain only their own reusable scratch, while stability directly owns and builds its Bee matrix. Verified and unsafe `false` results both mean no candidate under their respective proof/heuristic contracts, and `true` continues to the exact procedure. Release passed all 66 CTests and 26 wrapper tests; the verified-search oracle passed all 64 permitted tests with IDs 33-34 excluded, and ASan/UBSan passed the 11 core/CLI tests plus IDs 45-47. The verified source was also confirmed to compile with strict floating-point semantics, contraction and IPO/LTO disabled; the repository's fast timing sweep completed all 53 selected matrices successfully.
205. Reconciled verified candidate search with the merged PyFracESSA main branch:
preserved signed 64-bit matrix IDs, monotonic nanosecond results, and the canonical SQLite/PyFracESSA architecture while adding verified default, explicit unsafe, and exact search consistently across C++, CLI, Pybind, and `RunConfig`. Restored regression IDs 45-47 to SQLite, removed the development exact-rejection oracle and its CMake option, and retained the ordinary database sweep as the end-to-end result check. Release passed all 11 C++/CLI tests and all 53 PyFracESSA tests; a complete verified-mode run matched the stored ESS count for all 88 matrices; ASan/UBSan passed all 11 C++/CLI tests; SQLite integrity and foreign keys passed.
206. Preserved the retired Callgrind tooling and logo concepts:
restored the four deleted JSON-fed Callgrind scripts unchanged under `archive/callgrind/` as inactive historical reference, kept generated local `cpp/callgrind/callgrind.out.*` profiles ignored, and added the original glossy and geometrically corrected SVG/PNG simplex logo concepts under `experiments/`.
207. Published the retained historical Callgrind profiles:
added the complete local `cpp/callgrind/callgrind.out.1` through `.35` set to GitHub as 35 historical Callgrind 3.15.0 profiles from the former x86-64 Linux build while retaining the ignore rule for future generated profiles.
208. Made MatrixParser header-only:
moved the unchanged validating parser and its private helpers into inline definitions in `cpp/include/fracessa/matrix_parser.hpp`, deleted `cpp/src/matrix_parser.cpp`, and removed that source from CMake. All 11 C++/CLI tests and all 53 PyFracESSA tests pass.
209. Benchmarked historical-default-very-unsafe, current unsafe, and current verified search:
built raw-double main revision `32f61679da64` and current main revision `34e003168607` as Release/native/LTO CLIs, then measured all 87 retained matrices sequentially on CPU 2 with the one-second adaptive protocol. Historical-default-very-unsafe mismatches IDs 38-39, current unsafe mismatches IDs 45-47, and current verified matches all 87. Their summed average-call times are 140.814, 59.402, and 331.353 seconds; the per-matrix median ratios for historical-default-very-unsafe versus unsafe and verified are 0.977582 and 0.431289. The existing timing session now holds all 261 new rows plus its retained 168-row Pybind baseline; SQLite integrity and foreign keys pass.
210. Removed the dimension-one canonical fixture and benchmark observations:
deleted matrix 40, its candidate, and its five timing rows from SQLite; removed its rows from the retained support-frontier and three-mode timing reports and raw benchmark tables; and removed the dimension-one direct-bracelet measurements while retaining dimension-one unit verification. The canonical store now contains 87 matrices, 49,157 representatives, 86,152 weighted candidates, 83,377 weighted ESS, and 429 timing rows.
211. Replaced mixed-backend averages with one persistent-Pybind median benchmark:
changed the timing runner to size each sample from the first returned C++ nanosecond duration, include that first value, and store the median instead of using Python wall time or an arithmetic mean. Deleted all 429 former timing rows, then measured all 87 matrices on CPU 2 in current unsafe, verified, and exact modes from one loaded native module and historical-default-very-unsafe from a second loaded module. The temporary historical binding adapter changed only its parser selection and timer unit; no historical algorithm source changed. The new 348-row session has mismatch IDs 38-39 for historical very unsafe and 45-47 for current unsafe, while verified and exact match all matrices. Summed per-matrix medians are 138.484, 59.389, 347.148, and 1,437.953 seconds respectively. All 11 C++/CLI tests and all 54 PyFracESSA tests pass; SQLite integrity and foreign keys pass.
212. Benchmarked the release supplied to Werner in November 2016:
matched the backed-up private download page, email record, source ZIP, and retained `REF_2016-11-20-Werner` algorithm sources; built a temporary persistent Pybind adapter that changed only result delivery and native nanosecond timing; and measured all 87 matrices in Werner very-unsafe mode on CPU 2. It mismatches IDs 38-39 and sums to 234.824 seconds. Werner exact was attempted only where current exact took at most 30 seconds, then limited to two minutes per matrix: 72 correct rows completed and sum to 702.719 seconds; IDs 33, 45, 47, 55, 63-66, and 89-90 were not attempted, while IDs 32, 34, 52, 62, and 70 exceeded two minutes and were skipped. The timing session now contains 507 rows, the comparison report appends Werner very-unsafe and exact columns, and SQLite integrity and foreign keys pass.
213. Adopted a wider project line-width convention:
set 120 columns as the soft target for active C++, Python, comments, and Markdown; removed the traditional 80-column expectation while allowing longer formulas, commands, URLs, matrices, and expressions when wrapping would hurt readability.
214. Added the three Bomze-Schachinger-Ullrich papers:
stored the verified local final or print copies of the 2014 cp-rank paper for dimensions seven through eleven, the 2015 cp-rank bounds paper, and the 2018 Standard Quadratic Optimization Problem paper under `aidocs/reference/papers/` alongside the existing Bomze and Hadeler references.
215. Created a top-level research-paper collection:
moved all five retained mathematical papers, including the existing 1983 Hadeler and 1992 Bomze references, from `aidocs/reference/papers/` to `research/papers/` without renaming the files.
216. Added complete machine-readable transcriptions of the foundational papers:
transcribed Bomze's 1992 ESS-detection paper and Hadeler's 1983 copositive-matrix paper into Markdown beside their source PDFs, preserving the full prose, theorem and proof structure, equations as LaTeX, examples, acknowledgements, and references while omitting only repeated journal page furniture.
217. Re-audited both foundational paper transcriptions against the rendered PDFs:
checked every page, formula, index, theorem, proof, example, and reference; restored one omitted Bomze phrase,
matched two Bomze subset symbols to the print, corrected Hadeler equation (13) from `b_{iij}` to `b_{iji}`,
documented three apparent Bomze print errors and Hadeler's printed \(\beta\) error, and independently
verified the central algebraic identities with exact rational arithmetic.
218. Reintroduced the historical raw-double candidate procedure as an explicit production mode:
added `find_candidate_very_unsafe` with the unnormalized conversion, pivot cutoff, and margins from revision `32f61679da64`; replaced the separate exact/unsafe Boolean selectors with one `analysis_mode` across the analyzer, CLI, Pybind, `RunConfig`, and timing tool; and retained verified search as the default. The CLI now uses `--mode`, while Python uses `RunConfig(mode=...)`; accepted values are `verified`, `exact`, `unsafe`, and `very_unsafe`. Release passed all 11 C++/CLI tests and all 56 PyFracESSA tests. A persistent native sweep of all 87 retained matrices matched every expected ESS count except IDs 38-39, exactly reproducing the historical raw-double failures.
219. Moved exact Gaussian elimination to its sole production owner:
fused the unchanged pivoting, elimination, back-substitution, and positive-support checks directly after bordered-system construction in `find_candidate_exact.cpp`; removed the now-single-purpose `linear_solver.hpp`; and retained its four focused regressions through the public exact-candidate procedure. Release passed all 11 C++/CLI tests and all 56 PyFracESSA tests, all 24 representative pre/post candidate-output comparisons were byte-identical across the four numerical modes, and ASan/UBSan passed all 11 C++/CLI tests.
220. Replaced the exact bordered solve and common Bee factorization with one reduced-Hessian factorization:
eliminated the payoff/normalization border into the symmetric $(k-1)$-dimensional system $Hy=r$, added exact congruence-preserving rational $LDL^T$ with symmetric permutations and mandatory 2-by-2 pivots for nonsingular zero-diagonal Schur complements, reconstructed the exact candidate, and retained the factorization inertia. Stability now rejects immediately when $H$ is not negative definite and accepts immediately when $H\prec0$ and extended support equals support; only $H\prec0$ with outside best replies constructs Bee and uses the unchanged Bomze/copositivity fallback. Detailed source comments derive every reduction, pivot, solve, permutation, and stability implication. The focused zero-diagonal 2-by-2 regression, Release 11/11 C++/CLI tests, 56/56 Python tests, and ASan/UBSan 11/11 C++/CLI tests pass. Exact candidate output matched the former solver on 240 deterministic symmetric games and all 22 retained matrices exercising larger extended support; the canonical unsafe sweep matched all 84 expected candidate contracts while preserving known unsafe failures 45-47. The complete one-second persistent-Pybind exact benchmark on CPU 2 covered 85 matrices, with IDs 33-34 omitted from that run, and matched every expected ESS count. Summed medians fall from 1,386.743 to 1,184.045 seconds (14.62%); the median matrix improves by 31.58%, and 82 of 85 matrices are faster. IDs 45 and 47 regress from 74.655 to 168.882 seconds and from 72.686 to 132.230 seconds respectively; excluding those two cases, summed time improves by 28.76%. ID 46 is 3.8 milliseconds slower (13.4%); its absolute effect is small, but one repeat is needed before calling it signal or timing noise.
221. Benchmarked integer fraction-free exact candidate screening:
added an isolated five-way comparison of the old bordered rational Gaussian solver, current reduced-Hessian rational $LDL^T$, integer bordered FFLU, a conservative complete FFLU-plus-candidate-$LDL^T$ hybrid, and rational fraction-free FFLU without changing production source. The CPU-2 one-second sweep completed all 82 matrices below the common two-minute expectation cutoff, including ordinary IDs 33-34 and excluding IDs 45, 47, 65, 66, and 90. Every candidate contract matched; the complete hybrid also matched every exact negative-definiteness decision. Its summed medians fell to 83.230 seconds from 326.256 seconds for current $LDL^T$ and 458.322 seconds for old Gaussian, while adding only 0.34% to bare FFLU. It lost all 26 dimension-2-to-8 rows, won 25 of 28 medium rows, and won 27 of 28 large rows; the complete table and raw nanoseconds are retained with the experiment.
222. Added an experimental fraction-free symmetric candidate solver:
adapted FLINT FFLU's immediate-integer and arbitrary-precision Bareiss arithmetic to one specialized FFLDLT-style factorization of the integer border-reduced system $H y=r$. The factorization handles every symmetric nonsingular system, including indefinite systems and zero-diagonal Schur complements, solves the candidate with one common integer denominator, and obtains inertia from the same pivots. Indefinite inertia never rejects a candidate, so every valid candidate still reaches support-superset pruning. The complete CPU-2 one-second sweep reused the five stored baselines and measured only FFLDLT; summed medians were 58.543 seconds, 87.23% below old Gaussian, 82.06% below current rational $LDL^T$, and below the 82.952-second candidate-only FFLU total. Complete candidate and inertia comparisons passed on five retained matrices, 357 deterministic games, 105 adversarial indefinite block systems, structured initial and late zero-diagonal cases, and 48 ASan/UBSan games including forced arbitrary-precision inputs; all 82 retained ordered support/extended-support checksums matched the exact database. Production source remains unchanged.
223. Replaced the production rational candidate factorization with fraction-free integer $LDL^T$:
converted each exact game once to one integer matrix and positive common denominator, reused integer reduced-system, right-hand-side, solution, and scratch storage across supports, and called the extracted FLINT-style FFLDLT kernel for singularity, candidate solution, and inertia. Exact sign and outside-payoff checks now remain integer; public payoff and vector fractions are materialized only after candidate validation succeeds. Release passed all 11 C++/CLI and 56 Python tests, ASan/UBSan passed all 11 C++/CLI tests, and ten balanced matrices through dimension 23 matched the canonical candidate IDs, supports, extended supports, vectors, payoffs, ESS flags, and weighted ESS counts. A CPU-2 persistent-Pybind comparison against preserved rational revision `799be715` measured nine dimension-3-to-23 matrices: the median improvement was 78.33% and the arithmetic mean improvement was 56.32%; only the 0.8-microsecond dimension-4 circular case regressed, by 5.18%.
224. Removed the normalized unsafe mode and renamed the historical raw mode to unsafe:
isolated full-support diagnostics showed that normalization fixed raw-double failures 38-39 but introduced false rejections 45-47 by shrinking bordered pivots and changing tiny positive probabilities to negative approximations. Production now exposes only `verified`, `exact`, and raw `unsafe`; the separate `very_unsafe` spelling, class, sources, routing, and tests are removed. The raw algorithm remains identical to revision `32f61679da64`, its one-time input warnings remain available, and it again uses the normal Release/LTO flags while only verified search uses strict floating-point flags. Release passed all 11 C++/CLI tests and all 55 Python tests; focused CLI and native regressions verify the removed mode is rejected, raw unsafe still misses IDs 38-39, and raw unsafe retains ID 46.
225. Converted unsafe input warnings into silent matrix-wide exact fallback:
when any of the six one-time exact-to-double conversion checks triggers, the analyzer bypasses unsafe rejection for every generated support of that matrix and lets the existing exact candidate solver decide. The historical pivot cutoff and all other per-support unsafe logic remain unchanged. Removed the unconditional CLI warning; focused CLI and native regressions confirm that IDs 38-39 now match exact analysis without warning output. Release passed all 11 C++/CLI tests and all 55 Python tests.
226. Replaced the numerical-mode surface with required fast or safe search:
removed the verified proof implementation and its strict floating-point build target; renamed raw `find_candidate_unsafe` to `find_candidate_fast` and exact `find_candidate_exact` to `find_candidate_safe`; and replaced the old mode enum, flags, defaults, and `RunConfig.mode` with a required first `fast` or `safe` argument in the CLI, Pybind, sequential Python, and multiprocessing APIs. Fast retains its six matrix-wide safe-fallback checks and exact confirmation of every surviving support; safe starts directly with exact candidate solving. The timing tool now accepts `--method`, maps those two choices onto historical binary interfaces only inside the benchmark adapter, and preserves the SQLite column name while allowing new labels. Werner default and the preserved pre-mode default rows were relabeled `fast`, and Werner exact rows were relabeled `safe`; the later three-mode snapshot retains its historical labels. Release passed all 10 C++/CLI tests and all 56 Python tests; SQLite integrity and foreign keys pass.
227. Benchmarked the committed current fast method:
rebuilt revision `8697ebaf`, kept one Pybind process pinned to CPU 2, and measured every retained matrix with dimension at least 3 using the one-second target and median native nanoseconds. All 78 observed ESS counts match the canonical database. Dimension-2 matrices were excluded by benchmark policy. The stored module SHA-256 is `dfe70370e674c5e766067ab65526a5e3b58c484f581db36ed52d280c8c169b9f`; SQLite integrity passes.
228. Documented a current fast false rejection:
constructed an exact symmetric three-strategy game whose full-support candidate $(1/3,1/3,1/3)$ is an ESS, while the raw-double elimination produces a $7.5\times10^{-13}$ pivot and rejects it against the retained $10^{-12}$ cutoff. The matrix passes all six matrix-wide input checks, separating this failure from the input-conversion regressions at IDs 38-39 and the retired normalized-heuristic regressions at IDs 45-47. The exact proof and CLI reproduction are preserved in `correctness/FAST_CANDIDATE_FALSE_REJECTION.md`; the matrix is not yet in canonical SQLite data.
229. Completed the fast false-rejection construction for the remaining candidate conditions:
after changing only small pivots to safe fallback in an isolated Release build, constructed one exact three-strategy ESS whose
strictly positive $10^{-10}$ probability is computed as negative and one exact four-strategy ESS whose outside strategy is
strictly worse by $10^{-4}$ but is computed above the fast rejection margin. Branch-isolated runs recover each ESS only when its
corresponding probability or outside-payoff rejection is changed to safe fallback. All three exact constructions and
reproductions are preserved in `correctness/FAST_CANDIDATE_FALSE_REJECTION.md`; production source and canonical SQLite data remain
unchanged.
230. Added an independent experimental test candidate search:
copied fast double matrix storage, bordered Gaussian elimination, probability test, and outside-payoff test into new
`find_candidate_test` HPP/CPP files so experiments cannot change or share production fast behavior. Test search replaces fast's
six conversion checks with the exact integer precision span $P=M/m$ after clearing the matrix's least common positive
denominator, selects matrix-wide safe fallback before double conversion when $P\geq10^9$, and sends a support to exact checking
when either bordered pivot is below $10^{-12}$. The CLI, native binding, and maintained Python
execution APIs accept required method `test`; the canonical timing schema and stored database remain unchanged. Boundary tests
cover $P<10^9$, $P=10^9$, tiny values, and close large values. The three documented fast false-rejection games produce ESS
counts fast/test/safe of 0/1/1, 0/1/1, and 1/2/2. All 10 C++/CLI tests and all 56 Python tests pass.
231. Added all three current fast false-rejection counterexamples to canonical test data:
stored the pivot-cutoff, positive-probability, and outside-payoff games as non-circular IDs 91-93 with their complete exact
candidate rows. Their fast/test/safe ESS counts are `0/1/1`, `0/1/1`, and `1/2/2`. The database now contains 90 distinct
matrices, 49,161 candidate representatives, 86,156 weighted candidates, and 83,381 weighted ESS. SQLite integrity, foreign keys,
summary counts, support-size structures, and freshly recomputed safe candidate contracts pass. No timing rows were added.
232. Removed duplicate exact matrix preparation from experimental test search:
the precision-span decision now reads the safe solver's existing integer-scaled game and common denominator instead of allocating
another rational matrix, integer matrix, and denominator and repeating FLINT's denominator clearing. Matrix 11 remains correct at
14 ESS; five alternating one-second CPU-2 batches measured fast at 5.833 microseconds and test at 6.125 microseconds, leaving a
5.00% test overhead instead of the earlier roughly 15-17%. All 10 C++/CLI tests and all 56 Python tests pass, and the three fast
false-rejection regressions still produce the expected test ESS counts 1, 1, and 2.
233. Promoted the experimental test candidate logic to production fast search:
fast now uses the safe solver's existing integer game for the exact $P\geq10^9$ matrix-wide fallback, converts directly to the
unnormalized binary64 matrix otherwise, and treats every pivot below $10^{-12}$ as inconclusive so exact arithmetic decides that
support. Test remains an independent source copy and is currently behaviorally identical to fast. All 10 C++/CLI tests and all 56
Python tests pass. Regression IDs 91-93 now produce fast/test/safe ESS counts `1/1/1`, `1/1/1`, and `2/2/2`.
234. Added the C++ boundary around exact FLINT integers:
introduced header-only `linalg::integer` and `linalg::matrix_int` owning wrappers, changed the fraction-free $LDL^T$ interface and
safe candidate solver to use those C++ types, and confined raw FLINT calls to the wrappers and specialized numerical kernel. One
focused matrix test reconstructs every rational entry after common-denominator conversion. Release passed all 10 C++/CLI tests
and all 56 Python tests; ASan/UBSan passed all 10 C++/CLI tests; complete candidate output matched the pushed checkpoint on 11
representative and regression matrices. A balanced eight-matrix CPU-2 safe benchmark through dimensions 3-23 produced a median
runtime ratio of `0.967` to the checkpoint, so the inline readability layer showed no measurable overhead.
235. Rebuilt and benchmarked the exact pre/post C++ FLINT-wrapper revisions:
purged every timing family except Werner and the July 27 pre-refactor build, renamed that preserved build `classic`, and added
paired safe runs for adjacent clean revisions `3547df5d` and `29799de8`. The CPU-2 persistent-Pybind runs use one-second native
medians for 77 dimension-3-or-larger matrices; IDs 65-66 and 90 were not attempted, and ID 47 was removed after its pre-wrapper
run exceeded 30 seconds. All ESS counts match. The wrapper build wins 59 rows, ties 4, and loses 14; its mean and median changes
are $-7.01\%$ and $-8.50\%$, while summed medians fall by $6.27\%$ from 73.083 to 68.499 seconds. A reverse-order matrix-60 Perf
check confirmed fewer instructions, branches, and cycles. Disabling strict aliasing improved only the raw-pointer build and
reduced, but did not remove, the gap, isolating compiler code generation as the source and alias optimization as one contributor;
no mathematical algorithm or storage-copy operation changed. The canonical database now contains 400 timing rows and passes
integrity and foreign-key checks.
236. Fixed one canonical compiler configuration for performance comparisons:
documented the exact Release/native/LTO Ninja build through the system ccache compiler wrappers, prohibited method-specific and
unrequested floating-point, aliasing, sanitizer, profiling, debug, or optimization flags, and required both revisions in a source
comparison to use the same clean toolchain and dependency environment. Toolchain changes now start a separately recorded benchmark
lineage; only explicitly requested and labelled experiments may deviate from the canonical configuration.
237. Promoted the equilibrated reduced-system double solver to production fast search:
copied the independent test implementation into fast without introducing shared numerical code. Fast now removes the common game
denominator, applies the exact integer $P\geq10^9$ matrix-wide safe fallback, performs one common power-of-two normalization and
one LAPACK-derived symmetric BIN equilibration $A\mapsto DAD$, eliminates the candidate border, and solves each transformed
reduced Hessian with the adapted Bunch-Kaufman $LDL^T$ factorization. Test remains an independent source copy. The complete
CPU-2 benchmark of the test implementation measured a median per-matrix change of $+1.13\%$ against unscaled reduced-system
$LDL^T$; temporarily removing the $P$ gate changed the median by only $+0.11\%$ and was reverted. Release passed all 10 C++/CLI
tests and all 56 Python tests, ASan/UBSan passed all 10 C++/CLI tests, and fast reproduced every stored exact candidate contract
and ESS classification for all 90 canonical matrices. Fast and test now also send non-finite completed probability and payoff
accumulations to exact checking, while the CLI regression compares complete fast/test/safe candidate output for all three stored
historical false rejections. MSVC Release `/fp:fast` remains tracked separately because it does not guarantee IEEE special-value
behavior.
238. Preserved non-finite fast fallback on MSVC Release:
replaced `/fp:fast` with `/fp:precise`, because the candidate paths rely on IEEE NaN and infinity propagation plus
`std::isfinite()` to select exact fallback. Linux Release tests remain unchanged; Windows compilation and regression verification
are deferred to the next Windows release build.
239. Promoted the measured FP-S02 and FP-S03 fast-path simplifications:
fast now keeps only the outside-support complement mask until factorization and probability validation succeed, then walks its set
bits directly in ascending order. Exact-to-double conversion and final equilibration scaling each process one symmetric triangle
and mirror it. The temporary shared conversion flag was removed, and fast/test source differs only in class and header names. On
the balanced 81-matrix CPU-2 panel, FP-S02 alone improved the median by 0.81%; combined FP-S02+FP-S03 improved it by 1.91%, while
all ESS counts matched. Release passed all 10 C++/CLI and 56 Python tests, and one promoted-fast pass reproduced all 81 benchmarked
matrix ESS counts including IDs 45-47 and 91-93.
240. Extended the canonical matrix catalog metadata:
added names, family and subfamily classifications, descriptions, source URLs, original formats and IDs, and approximate
creation or first-use dates to all 87 matrix rows. Legacy REF fixtures are explicitly classified as historical, with January
1 used as a documented placeholder when only the year is known; no mathematical matrix, candidate, or timing data changed.
241. Imported the small symmetric SuiteSparse corpus:
selected all 12 real, square, numerically symmetric SuiteSparse matrices of dimension at most 30; converted Matrix Market
pattern, integer, and finite real tokens to exact FracESSA rational input; skipped SuiteSparse ID 2758 because its dense matrix
duplicates existing ID 1; and added the other 11 as IDs 91-101 with complete exact candidate baselines and source metadata.
Exact and verified candidate output matched for every import. Dimension 27 ID 92 carries the new `super_large` tag; the
canonical database now contains 98 distinct matrices, 49,388 representatives, 86,383 weighted candidates, and 83,462 ESS.
242. Audited the finite NIST Matrix Market repository under the same import rule:
the official order-below-31 square search returned seven files, but none produced a new distinct symmetric game. `CAN 24`
exactly duplicates ID 91; `FIDAP005` is the same source matrix as ID 92 with fewer printed digits; four files are unsymmetric;
and NIST withdrew its incorrect `LAP 25` Matrix Market assembly. IDs 91-92 now retain the NIST pages as alternate provenance
and carry `nist_matrix_market`; parameterized generators remain outside the finite downloadable-file audit.
243. Extended the exact matrix catalog through the parser's dimension-63 limit:
re-audited all 56 real square SuiteSparse files through order 63, retained 28 exact-symmetric sources, skipped ID 2758 as an
existing exact duplicate, and added 16 new catalog-only rows as IDs 102-117. Extended the finite NIST symmetric search through
order 63; its nine downloads add only exact or rounded SuiteSparse provenance, while `LAP 25` remains withdrawn. Audited all
136 fixed QAPLIB instances, tested A and B independently, retained 182 in-range symmetric occurrences, removed 30 internal
duplicates, and added 152 exact distinct matrices as IDs 118-269 with DOI and CC BY 4.0 attribution. Made the four baseline
summary fields nullable only as an all-or-none catalog state, and made the timing runner skip such rows by default or reject
their explicit selection. The canonical database now has 266 pairwise-distinct matrices, 98 complete baselines, 168 catalog-only
rows, 49,388 representatives, 86,383 weighted candidates, and 83,462 ESS. All 55 Python tests, exact source round-trips,
SQLite quick/integrity checks, foreign keys, and global dense deduplication pass; analyzer source is unchanged.
244. Imported the explicit symmetric TSPLIB95 subset:
audited all 111 files in the official symmetric-TSP archive, selected the 17 instances declaring
`EDGE_WEIGHT_TYPE: EXPLICIT`, excluded six above dimension 63, and added the 11 remaining exact integer matrices as
catalog-only IDs 270-280. Expanded all three retained TSPLIB storage layouts and matched every edge exactly against the
independent official XML files; all imports are non-circular and globally distinct. The canonical database now contains
277 matrices, 98 complete baselines, and 179 catalog-only rows; candidate and timing data are unchanged.
245. Imported the in-range Biq Mac Library corpus:
audited all 468 files in the official archive as 343 logical instances, verified all 125 dense/sparse pairs exactly, and
retained the 33 instances through dimension 63: 10 Beasley, 13 Glover-Kochenberger-Alidaee, and 10 Rudy Max-Cut matrices.
Expanded sparse symmetric-Q files and Max-Cut edge lists exactly, matched all 46 retained individual files byte-for-byte
against the aggregate archive, and found no internal or existing dense duplicates. Added the non-circular matrices as
catalog-only IDs 281-313. The database now contains 310 matrices, 98 complete baselines, and 212 catalog-only rows;
candidate and timing data are unchanged.
246. Imported the exact symmetric subset of Magma's Hadamard databases:
verified the ordinary and skew archives against Magma's published hashes, decoded their compact binary index/data format,
and cross-checked the decoder against the exact degree-16 handbook example. Audited all 4,474 ordinary representatives
through degree 63 and all 638 skew representatives; only six ordinary matrices are symmetric, one each at degrees 1, 2,
4, 8, 16, and 32, while no skew matrix is symmetric. Added the six globally distinct `+/-1` matrices as catalog-only IDs
314-319. The database now contains 316 matrices, 98 complete baselines, and 218 catalog-only rows; candidate and timing
data are unchanged.
247. Imported the in-range QPLIB quadratic matrix corpus:
used QPLIB's official statistics to select all 35 of 453 problems with at most 63 variables, downloaded their canonical
`.qplib` files, and parsed every explicitly stored objective and quadratic-constraint lower triangle with exact rational
arithmetic. Mirrored each lower triangle to its symmetric Hessian and independently matched all 2,475 matrix occurrences
against PyQPLIB 0.1.7. Exact dense deduplication collapsed 869 repeated occurrences, leaving 1,606 globally new,
non-circular catalog rows at IDs 320-1925 with complete alternate provenance. The database now contains 1,922 matrices,
98 complete baselines, and 1,824 catalog-only rows; candidate and timing data are unchanged.
248. Imported the in-range OR-Library matrix corpus:
audited every locally maintained family in the current index plus the still-hosted urban-transit page and retained only
explicitly stored exact rational square matrices with
dimension at most 63 and exact symmetry. The eligible source occurrences include binary-quadratic, capacitated spanning-tree,
aircraft-separation, CAB hub-location, portfolio-correlation, and urban-transit-demand matrices; recovered the now-missing
official `td1` and `td2` files from their 2011 Internet Archive snapshots. Removed ten repeated `capmst` occurrences and one
repeated portfolio occurrence, rejected every asymmetric or nonnumeric table, and added 57 globally new non-circular
catalog rows as IDs 1926-1982. The database now contains 1,979 matrices, 98 complete baselines, and 1,881 catalog-only rows;
candidate and timing data are unchanged.
249. Audited the COMPl_e_ib 1.1 control benchmark library:
constructed its 168 state matrices from the official MATLAB source and data files, excluded 57 above FracESSA's dimension-63
limit, and found that none of the 111 in-range state matrices is exactly symmetric. Control-channel arrays and synthesized
identity, zero, and weighting matrices are not independent benchmark matrices, so no catalog row or database value changed.
250. Audited the SLICOT model-reduction benchmark collection:
excluded 17 of its 18 linear-system models because their orders exceed 63, then checked the sole in-range order-48 building
model and found its state matrix `A` asymmetric. Its official `build.mat` is byte-identical to COMPl_e_ib's already-audited
`lah.mat`, so no catalog row or database value changed.
251. Imported the in-range KONECT adjacency-matrix corpus:
downloaded and exactly reconstructed all 23 available unipartite network files through dimension 63, retained all nine
undirected matrices plus the directed but reciprocal `moreno_taro` matrix, and rejected the other 13 directed matrices as
asymmetric. Excluded seven rectangular bipartite files, mapped Dolphins and Zachary karate club to exact existing IDs 114-115
as alternate provenance, and added the eight globally new non-circular matrices as catalog-only IDs 1983-1990. The database
now contains 1,987 matrices, 98 complete baselines, and 1,889 catalog-only rows; candidate and timing data are unchanged.
252. Added a deterministic House of Graphs stratified sample:
audited all 23,988 canonical graphs with dimensions 1-63 against their API adjacency lists, crossed five dimension bands with
ten structural/control categories, and selected up to three per populated stratum by SHA-256 rank with seed `20260802`. One
45-63 disconnected-cyclic stratum contained only one graph and one graph overlapped two strata, leaving 147 unique exact
zero-diagonal `0/1` adjacency matrices. None duplicates an existing dense matrix; eight use compact circular storage. Added
them as catalog-only IDs 1991-2137 with source IDs, canonical graph6 values, names, selected strata, and source links. The
database now contains 2,134 matrices, 98 complete baselines, and 2,036 catalog-only rows; candidate and timing data are
unchanged.
253. Reduced the QPLIB corpus to objective matrices only:
removed all 1,576 distinct rows whose only retained role was a quadratic constraint and preserved the 30 explicitly stored
quadratic objectives. Exact comparison confirmed that the objectives are pairwise distinct and duplicate no pre-QPLIB row;
no removed constraint-only row matched the pre-QPLIB database. Existing matrix IDs were not renumbered, so the resulting gaps
between IDs 320 and 1925 are intentional. The database now contains 558 matrices, 98 complete baselines, and 460 catalog-only
rows; candidate and timing data are unchanged. Recorded objective-only import as the general policy for optimization-problem
collections while preserving independently sourced copies of the same mathematical matrix.
254. Added a deterministic Network Data Repository sample:
audited the current 6,628-row graph-category index, identified 1,241 entries reporting dimensions 1-63, crossed their source
categories with five dimension bands, and ranked candidates by SHA-256 with seed `20260802`. Inspected the downloaded source
archives directly and accepted only square exact-symmetric Matrix Market files or edge lists explicitly undirected or exactly
reciprocal; rejected rectangular and unsymmetric data, temporal streams, bipartite tables without a square adjacency matrix,
malformed files, and every case requiring forced symmetrization. Skipped exact existing and within-sample duplicates and
added 39 globally new catalog rows at IDs 2138-2176: 15 animal-social, 15 cheminformatics, six protein, two DIMACS, and one
biological matrix. All source archives round-trip exactly; one matrix uses compact circular storage. The database now contains
597 matrices, 98 complete baselines, and 499 catalog-only rows; candidate and timing data are unchanged.
255. Imported SDPLIB objective matrices only:
audited all 92 problems in the SDPLIB 1.2 mirror, selected the 30 whose complete block-diagonal dimensions are at most 63,
and retained exactly one `F0` objective coefficient matrix from each. Parsed SDPA sparse integer, decimal, and scientific
tokens as exact fractions and expanded the source blocks into complete symmetric matrices. Excluded all 1,799 `F1...Fm`
constraint matrices and did not catalog individual blocks. The 30 objectives are pairwise and globally distinct and occupy
catalog-only IDs 2177-2206; 15 carry the `super_large` tag. The database now contains 627 matrices, 98 complete baselines,
and 529 catalog-only rows; candidate and timing data are unchanged.
256. Combined the expanded matrix catalog with the candidate-search optimization branch:
retained all 627 catalog rows from `main`, moved the branch's three fast false-rejection regressions from colliding IDs 91-93
to IDs 2207-2209, remapped their candidate and timing references, and retained the branch's curated 724 timing rows instead of
restoring superseded benchmark sessions. The merged database contains 630 pairwise-distinct matrices, 101 complete baselines,
529 catalog-only rows, 49,392 candidate representatives, 86,387 weighted candidates, and 83,466 weighted ESS. SQLite integrity
and foreign-key checks pass.
257. Added exact-source Markdown transcriptions of the three Bomze-Schachinger-Ullrich papers:
located the LaTeX sources that generated the retained 2014, 2015, and 2018 PDFs and copied them beside the papers; extracted the
complete prose, theorem/proof structure, formulas, tables, figures, acknowledgments, and references into Markdown; and added seven
rendered figure assets. Audited all 77 PDF pages visually, checked 3,191 LaTeX math spans for balanced syntax, reconciled all 117
mapped theorem and equation cross-references with the exact sources and printed PDFs, restored the dropped proof labels, and
corrected the 2018 display sequence to the printed equation numbers 1-27.
258. Independently re-audited the three TeX-backed paper transcriptions:
compared their mathematical and prose streams, theorem/proof structure, equation tags, tables, figures, footnotes, citations, and
references against both the exact TeX sources and retained PDFs; removed one leaked `\cline{1-6}` marker from the 2015 Table 1
transcription; and corrected references 3-6 and their citations in the 2018 transcription after confirming that an auxiliary file
from a different draft had supplied the wrong numbering. Restored the 2014 corresponding-author details and the printed 2015
classification heading found during the title-page pass.
259. Re-audited all five research-paper transcriptions as one collection:
checked all 105 printed pages again, re-ran exact-source comparisons for the three TeX-backed papers, verified every
numbered-equation sequence, all 95 mapped citations, all 120 bibliography entries, both reconstructed numerical tables, all seven
figure assets, local links, footnotes, and Markdown parsing, and found no further transcription discrepancies.
260. Added explicit per-matrix timing calibrations:
stored fast and safe calibration durations as exact integer nanoseconds on `matrices` and made the timing runner use them instead
of its first measured sample to choose the iteration count. Fast now covers all 630 matrices with 406 positive values and 224
timeouts after retaining existing values and applying the one-second cutoff only to missing rows; safe retains 77 positive values.
Missing calibrations stop a run before measurement, the default target is 0.5 seconds, and calibration values remain manually
maintained matrix metadata rather than benchmark-history rows. A `-1` calibration records a process killed at its cutoff and
selects one benchmark iteration.
261. Completed safe calibration and retained its reusable runner:
added `scripts/calibrate_matrices.py` to fill only null fast or safe calibration fields on CPU 2 with a one-second default cutoff.
Safe mode also writes missing exact weighted counts, support-size structures, and representative candidate rows atomically. The
safe pass processed all 553 missing values, producing 364 positive calibrations and 266 cutoff sentinels overall and adding 267
complete baselines. The database now has 368 analyzed matrices, 57,012 representatives, 94,048 weighted candidates, and 85,305
weighted ESS. Pre-existing calibration, baseline, candidate, timing, and catalog data remained unchanged by value.
262. Normalized every eligible circulant matrix into compact zero-diagonal storage:
added `scripts/normalize_circular_matrices.py`, which audits exact fractions and converts a circulant source matrix by storing
`A - dJ`, where `d` is its common diagonal value, while recording that exact value in the matrix description. Converted IDs 1,
38, 39, 41, 43, 44, and 2203 with `d = 0, 0, 100000000000000000000, -1, 0, -1, 1`; documented dimension-one ID 314 as the sole
full-storage circular exception. Recomputed the seven baselines and fast/safe calibrations: weighted totals remain 94,048
candidates and 85,305 ESS, while orbit compression reduces stored representatives from 57,012 to 56,962. Removed 18 historical
timings measured against the old stored matrices and added 14 current-build fast/safe normalization measurements. The database
now contains 630 source/catalog identities, 628 distinct normalized stored matrices, 48 compact circular rows, and 720 timings.
263. Removed the newer normalized duplicate matrices:
IDs 39 and 41 both became the same compact dimension-two game as older ID 1 after subtracting their respective common diagonal
values. Removed the two newer matrix rows and their cascading candidate, calibration, and four current timing records, leaving 628
pairwise-distinct stored matrices, 56,960 representatives, 94,046 weighted candidates, 85,303 weighted ESS, 46 compact circular
rows, and 716 timings. The normalization script now identifies identical compact results and retains the oldest `created_at` row.
264. Retried every safe calibration timeout with a two-minute cutoff:
extended `scripts/calibrate_matrices.py` with a safe-only `--retry-timeouts` pass that selects `safe_calibration_ns = -1`, runs each
matrix exactly once, stores a successful exact baseline atomically when missing, leaves unsuccessful rows at `-1`, and continues
after an individual worker exits. The pass raised positive safe calibrations from 362 to 419, reduced `-1` rows from 266 to 209,
and added 56 complete baselines. The database now has 422 analyzed matrices, 61,618 candidate representatives, 98,704 weighted
candidates, and 86,507 weighted ESS. All stored candidate weights, ESS weights, support-size structures, SQLite integrity, and
foreign keys validate.
265. Imported the combined exact matrix-generator catalogue from current `main` without duplicating strategically equivalent games:
copied the 478 raw Anymatrix, TypedMatrices.jl, and Matrix Depot rows at IDs 2210-2687 and the generator tags for exact existing
IDs 48, 49, 314, 1995, 2001, and 2155. Exact circular normalization converted 53 new full matrices to compact `A - dJ` storage;
normalized ID 2215 then duplicated ID 44 and was removed, leaving 477 new catalog-only rows, including 80 compact rows and five
documented dimension-one exceptions. Added the native SQLite unique index on `(dimension, matrix)`; `matrix` alone is not unique
because compact strings omit dimension. The final database has 1,105 distinct normalized matrices, 126 compact rows, 422 analyzed
rows, and 683 catalog-only rows. Existing candidates, timings, baselines, and calibrations are unchanged; both calibration fields
are null on every new row. Exact source round-trips, global normalized deduplication, schema shape, tags, SQLite integrity, foreign
keys, and a second idempotent circular-normalization pass all validate.
266. Moved the dimension-26-through-63 classification out of tags and into `size_class`:
added `super_large` as the fourth size class, narrowed `large` to dimensions 17-25, migrated all 501 affected rows, and removed the
redundant `super_large` tag from each. The timing CLI now accepts `--size-class super_large`; its parser regression also uses the
new value. Updated the matrix-selection fixture to use distinct inputs under the new `(dimension, matrix)` uniqueness constraint.
All 10 C++/CLI tests and all 60 Python tests pass; SQLite integrity, foreign keys, and exact preservation of all nonclassification
matrix fields, candidates, and timings validate.
267. Calibrated fast search for all 477 generator-catalogue matrices:
ran the current Release Pybind module on CPU 2 with the one-second per-matrix cutoff, storing each result immediately. The pass
finished in 3.24 minutes with 337 positive native durations and 140 `-1` timeout sentinels. Fast calibration now covers all 1,105
matrices with 741 measurements and 364 timeouts; no fast calibration remains null. Safe calibrations and all candidate baselines
remain unchanged.
268. Calibrated safe search and exact baselines for all 477 generator-catalogue matrices:
ran the same Release Pybind module on CPU 2 with the one-second per-matrix cutoff and atomic per-matrix writes. The 3.60-minute pass
stored 332 positive native durations and their complete candidate rows, weighted candidate and ESS totals, and both support-size
structures; 145 timeout rows received `-1` with no partial baseline. Safe calibration now covers all 1,105 matrices with 751
measurements and 354 timeouts. The database has 754 complete baselines, 66,542 candidate representatives, 104,888 weighted
candidates, and 91,733 weighted ESS. Summary totals, support structures, candidate rows, SQLite integrity, and foreign keys all
validate exactly.
269. Removed positive rational scalar duplicates from the canonical matrix store:
compared exact stored value vectors only within equal `(dimension, is_cs)` classes and retained the lowest matrix ID. Deleted IDs
38, 44, 2212-2214, 2220, and 2254 in favor of IDs 1, 314, 2211, 68, and 2250; seven candidate rows and four normalization timing
rows cascaded with them. Negative scalar multiples remain distinct because they reverse payoff comparisons. The resulting 1,098
matrices contain 747 complete baselines, 66,535 candidate representatives, 104,875 weighted candidates, 91,720 weighted ESS, and
712 timing rows. Replaced deleted integration fixture ID 38 with precision-span regression ID 2208. Exact scalar re-audit finds no
remaining positive-scale group; the focused native regression, SQLite integrity, and foreign keys pass.
270. Consolidated dimension-one and compact all-zero corpus representatives:
retained ID 314 as the sole dimension-one game and removed equivalent IDs 2210-2211. Retained the largest compact all-zero game,
dimension-50 ID 2203, as the zero-payoff degeneracy/performance representative and removed IDs 43, 158, 1992, 2080, and 2127.
The cross-dimension zero-game rule is an explicit sampling choice, not a mathematical equivalence claim. Seven candidate rows and
two normalization timing rows cascaded with the deletions. The resulting 1,091 matrices contain 740 complete baselines, 66,528
candidate representatives, 104,826 weighted candidates, 91,718 weighted ESS, and 710 timing rows. Exact audits leave only ID 314
at dimension one and ID 2203 among compact all-zero games; SQLite integrity and foreign keys pass.
271. Consolidated diagonal benchmark-corpus representatives:
retained dimension-one ID 314, compact all-zero dimension-50 ID 2203, and nonzero dimension-60 ID 2180. Removed the other 27
diagonal matrices: 20 Strakos generator rows, six SDPLIB objective rows, and SuiteSparse ID 105. This is an intentional sampling
choice, not a mathematical-equivalence rule for diagonal games. Foreign-key cascades removed 728 candidate representatives
representing 728 candidates and 584 ESS; no timing row was affected. The resulting 1,064 matrices contain 713 complete baselines,
65,800 candidate representatives, 104,098 weighted candidates, 91,134 weighted ESS, and 710 timing rows. Exact re-audit leaves
only diagonal IDs 314, 2180, and 2203.
272. Added the paper-style ESS growth bound to matrix metadata:
`matrices.gamma_lower_bound` is a generated real column equal to `ess_count ** (1 / dimension)` and null when the exact ESS count
is unknown. The Bomze-Schachinger-Ullrich paper calls it the lower bound for $\gamma$ implied by the matrix and labels the result
column `$\gamma \geq$`. All 713 analyzed matrices expose a value without duplicating mutable data; `matrix_overview` places the
column immediately after `ess_structure`.
273. Added the eight missing payoff games from the 2014 and 2015 Bomze-Schachinger-Ullrich copositive constructions:
IDs 2688-2695 cover dimensions 6, 7, 8, 9, 11, 12, 13, and 14 and store the primitive integer transformation
`A = (dJ - S) / g`. Safe analysis records `8, 14, 20, 30, 33, 60, 108, 192` ESS; the order-8 and order-9 games have two and three
more ESS than the papers' copositive-zero counts. All eight safe and fast calibrations completed. The database now contains 1,072
matrices, 721 complete baselines, 66,195 candidate representatives, 104,563 weighted candidates, and 91,599 weighted ESS.
274. Made every whole-matrix fast-to-safe switch explicit:
added the nullable `safe_fallback` result and database field with the exact values `precision_span`, `double_conversion`, and
`equilibration`. CLI timing prints the value on line 3; Pybind, Python summaries, CSV/JSON/Parquet sinks, calibration, timing rows,
and `matrix_overview` preserve it, with the view placing it immediately after `gamma_lower_bound`. Calibration classifies before
starting its cutoff process, and timing ratio summaries exclude
fast rows that actually measured safe search. Classified all 1,072 stored matrices: 960 have no whole-matrix fallback, 45 use
`precision_span`, and 67 use `equilibration`; no current matrix hits `double_conversion`. Backfilled the field only for timing
panels that used the current precision-span/equilibration gate, leaving classic and Werner historical rows unchanged.
All 10 C++/CLI tests, all 63 Python tests, the 1,072-row native classification audit, SQLite integrity, and foreign keys pass.
275. Limited the fast/test equilibration zero-row exception to the affected coordinates:
exact all-zero rows and columns now retain scale one while the LAPACK-derived BIN iteration continues over every nonzero
coordinate with its reduced active dimension. Invalid arithmetic or failure to converge in that active block still switches the
whole game to safe search. A focused regression covers both a harmless isolated zero coordinate and a zero coordinate beside a
separately failing active block. Reclassification released 54 former zero-row fallbacks and identified 59 games that previously
used the last scaling after reaching the 100-iteration limit. The resulting matrix classifications are 955 without whole-game
fallback, 45 `precision_span`, four `equilibration`, and 68 `equilibration_non_convergence`; nine of the latter had already fallen
back generically because their zero rows prevented the active block from being examined. No database field except
`matrices.safe_fallback` changed.
276. Renamed the current invalid-equilibration fallback to `equilibration_invalid`:
current fast/test output and matrix classifications now distinguish invalid or non-finite equilibration arithmetic from the
100-iteration `equilibration_non_convergence` limit. Both reasons switch the whole game to safe search. The 13 historical timing
rows produced by older generic-fallback binaries retain `equilibration` rather than being relabelled inaccurately.
277. Removed the unreachable `double_conversion` fallback:
the preceding exact precision-span gate guarantees every nonzero power-of-two-normalized integer entry remains within roughly 30
binary exponents of the maximum when `P < 10^9`, so exponent overflow, binary64 underflow, and non-finite conversion cannot occur.
Current output and schemas now expose only `precision_span`, `equilibration_invalid`, and
`equilibration_non_convergence`; historical generic `equilibration` timing rows remain supported.
278. Made timeout calibration retries symmetric and resumable:
`scripts/calibrate_matrices.py fast|safe --retry-timeouts` now runs only the selected method's `-1` rows once and stores missing
candidate data atomically with each completed calibration. Fast results are exact only when a whole-matrix `safe_fallback` occurs;
otherwise a remaining `safe_calibration_ns = -1` marks them as unverified. The first two-minute fast pass was stopped after six
attempts: IDs 274, 2477, and 216 completed with candidate data, while IDs 2475, 2478, and 217 timed out. There are 319 fast
timeouts remaining.
279. Documented the project worktree policy:
work stays in the main worktree unless Reinhard explicitly approves using another worktree first.
280. Completed the two-minute fast-timeout retry with explicit multicore scheduling:
repeating `--cpu ID` now keeps one matrix pinned to each selected core, assigns the next matrix as that core becomes free, and
leaves SQLite writes in the main process. The pass attempted all 319 previous timeouts; two rows completed before reserving CPU 2,
and the remaining 317 attempts used performance CPUs 3 through 9. Twenty-one matrices completed and added 683 representative rows
for 841 weighted candidates and 236 ESS. Two completed through exact fallback and 19 remain unverified fast results. The other 298
rows timed out: 41 were classified for exact fallback and 257 remained on the fast path. Fast calibration now contains 774
positive durations and 298 timeout sentinels across all 1,072 matrices.
281. Completed the circular support generator V1-versus-V2 comparison:
verified identical candidates, ESS classifications, multipliers, fallback results, and order on all 81 quick-test matrices, plus
independent arbitrary-forbidden, multiplier, dimension-63, and sanitizer checks. The CPU-2 Release fast-mode benchmark timed the
33 circular quick cases with a 0.5-second native-duration target. Compact V2 was 41.48% slower at the median and 91.40% slower by
geometric mean; V1 won 28 cases and matrix 34 was 6.824 times faster with V1. Retained expanded-orbit V1 in production and recorded
the complete method and table in `experiments/circular_support_v2_2026-08-02/README.md`.
282. Implemented the 2014 direct fixed-density bracelet algorithm as experimental generator V3:
copied the source paper into `research/papers/`, extracted its mathematics and binary specialization into Markdown, and independently adapted `BraceFD` in the existing isolated generator harness. V3 matches both independent generators for every support size through dimension 24, emits canonical representatives once and in ascending order, and passes optimized plus ASan/UBSan verification. On the saved CPU-2 three-second panel for dimensions 5-24, V3 is 55.30% faster than V1 at the median and 53.31% faster by geometric mean; at dimension 24 it is 2.86 times as fast. Production remains unchanged pending an end-to-end pruning benchmark.
283. Promoted V3 to the production circular support generator while retaining V1 and V2:
integrated V1's callback, expanded forbidden-orbit pruning, and multiplier contract into direct binary `BraceFD`, then changed the
single circular runtime construction to `CircularSupportGeneratorV3`. V1/V3 generation matches through dimension 24 and V3 covers
dimension 63; all 81 quick matrices produced byte-identical fast and safe candidate output. Release and ASan/UBSan CTests plus all
63 Python tests passed. Two CPU-2 persistent-Pybind panels timed 31 circular quick cases; the conservative reverse-order panel was
19.90% faster at the median overall and 34.70% faster for dimensions 19 and above.
284. Added resumable full consistency calibration:
`matrices.calibration_timestamp` records the latest full audit and `calibration_comment` retains an append-only JSON history.
`scripts/calibrate_matrices.py audit` processes matrices by dimension and matrix ID, targets one second of native samples, uses a
two-minute minimum or 120% of the previous per-method calibration as its cutoff, validates complete candidate data, fills missing
results, and preserves mismatches in the log. Repeated CPU options run matrices concurrently while keeping a matrix's fast and safe
calls pinned to the same performance core. Fast runs precede safe runs; a fast timeout skips safe, while a whole-matrix safe fallback
supplies both calibrations.
285. Added the automorphism-group and gamma research note:
derived the published dimension-24 record as a 3-by-8 rook matrix with complete group `S_3 x S_8`, verified that its 15,120 ESS
form one full-group orbit, and recorded its quadratic form, spectrum, the general rectangular rook construction, and untested
high-symmetry inputs through dimension 30. The note keeps maximizing `gamma = ESS_count^(1/n)` as the objective and treats group
size only as a means of reducing equivalent support calculations.
286. Added the dimension-2-to-30 rectangular-rook automorphism ranking:
ranked every nondegenerate factorization `n = m * k` by the exact full group order, including the extra transposition symmetry for
square grids. Dimensions with fewer than three such factorizations are marked explicitly, and the note separates this finite rook
family from the broader problem of classifying all possible matrix automorphism groups.
287. Added exact symbolic matrix templates to the rectangular-rook ranking:
proved that two distinct off-diagonal variables are minimal, recorded every available compact circular string, and marked the
non-coprime grids that require full symmetric input. The templates use `A` for the same-row-or-column relation and `B` otherwise.
288. Added the general three-variable rook templates:
renamed the minimum representation to `P`/`Q` and added `A`/`B`/`C` compact strings that distinguish same-row, same-column, and
unrelated cells. Square grids are marked separately because three distinct values remove the row-column transposition symmetry.
289. Added the theoretical one-orbit gamma ceiling to every rook group:
computed the largest possible support orbit after the smallest unavoidable stabilizer, rather than assuming every group element
always gives a distinct support. Each row now records that orbit size and its `orbit_size^(1/n)` value.
290. Added the absolute theoretical gamma ceiling for dimensions 2 through 30:
recorded the independent Sperner bound
`gamma_max(n) = binomial(n, floor(n/2))^(1/n)` beside every group-specific one-orbit ceiling.
291. Added the largest ESS count currently found at each dimension:
queried all matrix families in the canonical database and placed `MAX(matrices.ess_count)` for dimensions 2 through 30 beside
the group-specific orbit ceiling and the absolute Sperner ceiling.
292. Implemented the opt-in cyclic symmetry filter:
added exact affine multiplier detection for circular matrices, filtered multiplier-equivalent V3 bracelets before solving, and
reused V3's existing dihedral expansion to preserve enlarged-orbit pruning and exact output multipliers. The feature is disabled
by default and exposed as CLI `--cyclic-symmetry-filter`, Python `RunConfig.cyclic_symmetry_filter`, and a final optional native
constructor parameter. All 10 C++/CLI tests, all 63 Python tests, and the ASan/UBSan CTests pass. Independent fast and safe audits
of all 33 circular quick-test matrices expanded both output forms to identical exact candidate-support classifications; 18 of the
66 method/matrix runs had extra affine multipliers. Matrix 34 preserves 15,120 candidates and 15,120 ESS while compressing 345
dihedral representative rows to 93 affine representative rows.
293. Replaced the shared linear-time support reflection with a fixed 64-bit reversal:
kept `bs64::reflect()` inline, used six constant bit-swap stages plus one alignment shift, and checked it against the direct
definition for every dimension from 0 through 64. Release and ASan/UBSan suites pass, and both canonical binaries produced
identical candidate output for all 200 cyclic-symmetry experiment matrices in fast and safe mode. Two half-second CPU-2 passes
with reversed build order measured combined per-matrix median reductions of 0.50%/0.67% on ordinary circular fast/safe cases and
1.54%/1.41% when the cyclic symmetry filter found extra multipliers; the dimension-50 case improved by about 7.6% in the
reverse-order pass. The absence of an ordinary-case regression and the repeatable extra-symmetry gain justified retaining it.
294. Made the cyclic symmetry filter automatic for circular matrices:
removed the CLI flag, Python `RunConfig` field, native binding argument, and C++ constructor configuration. Every circular matrix
now performs exact multiplier detection once; the helper disengages immediately when it finds only the identity class. The 101
no-extra cases among 120 stored circular matrices and completed fast/safe median on/off ratios of 1.0000 justified deleting the
experimental switch. Non-circular matrices remain unchanged.
295. Preserved the universal dihedral candidate-output contract under the affine filter:
each solved affine orbit now reconstructs one exact candidate row per distinct rotation/reflection orbit, and every row keeps its
own dihedral multiplier bounded by twice the dimension. The affine group remains internal to solve reduction and pruning; it is
never folded into an ambiguous matrix-specific output multiplier. All 200 generated extra-symmetry matrices match the former
filter-off candidate output byte for byte in fast and safe mode, and all 19 stored extra-symmetry database baselines remain exact.
296. Adopted calendar versioning with release-tag enforcement:
made `cpp/CMakeLists.txt` the single version source using `YEAR.MONTH.DAY.RELEASE_OF_DAY`, compiled that value into the CLI instead
of keeping a second hardcoded version, and made tagged CI reject any `v*` tag whose value differs from the built executable.
297. Aligned logged candidates and CLI floating-point output:
added a black-box affine-symmetry regression that requires the final logged candidate rows and IDs to equal sorted CLI output, and
replaced six-decimal `payoff_dbl` formatting with the standard round-trip binary64 precision in both output paths.
298. Simplified the result status contract:
removed the redundant `success` Boolean from native results, the Python wrapper, sinks, tests, and public documentation; callers and
the calibration script now use the integer `status` plus `error_message` contract. All 63 Python tests and the focused calibration
contract check pass.
299. Separated current documentation from retained history and research:
reorganized the documentation index by role, refreshed the live API, database, timing, calibration, and validation facts, and made
review ownership explicit. Failed numerical and performance approaches remain preserved as historical evidence, while obsolete
feature-branch instructions and the completed handover checklist were removed.
300. Re-audited the documentation cleanup against source, SQLite, and the mathematical references:
corrected the current V2 generator status and the $n=2$ exception to the cyclic multiplier-class bound, and recorded the inaccurate
source comment that says only ESS supports trigger pruning. The implementation already prunes after every exact equilibrium, as
justified by Bomze's support criterion; no production source changed.
301. Corrected documentation-cleanup defects found by the mathematical re-audit:
fixed the research-paper directory and indexed the two rook-symmetry notes, removed one remaining stale `test-only` label from the
unused V2 generator, and described the copositive construction precisely as producing isolated global, therefore strict local,
maximizers and ESS.
302. Recorded the final V2 support-generator decision:
V2 was a correct experiment but slower than both V1 and production V3, never entered production, and remains intentionally in the
source as evidence of a rejected design. Removed the conflicting review recommendation to delete it.
303. Corrected the C++ search-orchestration comment:
documented that every exact equilibrium support triggers strict-superset pruning, while stability classification determines only
whether the candidate is an ESS. Removed the resolved documentation-only finding from the C++ review; runtime behavior is unchanged.
304. Retired the completed review folder:
removed duplicate Pybind and Python validation ledgers, moved the unique C++ and fast-pipeline benchmark decisions into explicitly
dated `history/` records, and moved the three still-valid low-priority dead-code notes into `KNOWLEDGE.md`. Removed stale startup
routing and repaired all live cross-references; no source code changed.
305. Derived and recorded the two compact three-by-$k$ rook ESS sequences:
interpreted supports as spanning trees of $K_{3,k}$, proved the cycle obstruction for nonnegative interaction curvature,
counted the $k=3s+1$ and $k=3s+2$ support orbits, and exactly verified representative ESS through game dimension 60. The new
research note distinguishes constructed lower bounds from complete counts and records the slight dimension-30 gamma improvement.
306. Removed obsolete raw Callgrind output:
deleted the 35 tracked Callgrind 3.15.0 profiles from the former x86-64 build and the ignored local profiling log. The inactive
historical scripts remain under `archive/callgrind/`, and the ignore rule continues to exclude newly generated profiles.
307. Removed the remaining legacy Callgrind output:
deleted the two Callgrind 3.13.0 profiles stored with the 2019 REF snapshot; they described an obsolete executable and retained no
current performance evidence.
308. Documented the complete search workflow in the public README:
summarized support generation, circular symmetry reduction, exact-equilibrium superset pruning, reduced-Hessian solving and inertia,
and the retained Bomze copositivity fallback in one algorithmic overview. Added a plain-language definition of circular-symmetric
matrices and explained why their high-ESS examples motivate specialized support generation.
309. Replaced the public FracESSA logo:
promoted the V6 simplex concept with the separated origin, three axis intercepts, blue `ESS`, black final `A`, and optically balanced
wordmark spacing to `logo.png`. The README already references this file.
310. Made releases manual, cache-stable, and tag-after-test:
changed the release workflow to run only by manual dispatch from `main`, reject existing calendar-version tags, reuse default-branch
vcpkg binary caches across releases, and create the tag and GitHub release only after every standalone and Python build succeeds.
PyPI trusted publishing now follows the GitHub release, while ordinary pushes and pull requests continue to run no Actions.
311. Made research and experiment material local-only:
ignored the complete `research/` and `experiments/` directories and removed their previously tracked files from the current repository
tree while preserving every file in the local worktree. Existing historical commits and tags remain unchanged.
312. Simplified standalone release downloads:
replaced each platform archive and its duplicated documentation folder with one directly downloadable executable. The GitHub release
body now links the GPL license, exact tagged source archive, and third-party notices; Linux and macOS users receive the required
`chmod +x` instruction. Python wheel packaging is unchanged.
313. Removed redundant public maintenance wrappers:
made the complete root `scripts/` directory local-only while preserving its maintenance helpers in the worktree, and deleted the
redundant root `build.sh` and `test.sh` wrappers. Public build and test instructions now use their underlying CMake, CTest, and Python
commands directly.
314. Replaced the README's opening problem statement:
made the paper's canonical Standard Quadratic Optimization Problem the first prose sentence, displayed its numbered objective, and
defined the standard simplex separately with GitHub-compatible mathematics.
315. Fixed the README simplex delimiters:
replaced escaped dynamic braces, which GitHub reduced to an invalid `\\left{`, with portable `\\lbrace` and `\\rbrace` math
delimiters.
316. Connected the README's optimization and game interpretations:
explained that the same symmetric rational matrix defines both the StQP objective and a symmetric partnership game's payoff, and that
its ESS are precisely the strict local maximizers of the StQP objective on the simplex.
317. Removed an unnecessary README equation reference:
stated the ESS equivalence directly in terms of the quadratic form and removed the display equation's unused numeric tag.
318. Aligned the Python distribution and import names:
renamed the maintained package from `fracessa` to `pyfracessa`, updated wheel and source-distribution manifests, tests, timing
commands, examples, and active documentation, and retained `fracessa_core` as the private native implementation module. The old
`fracessa` import is not preserved as a compatibility alias.
319. Refreshed the top-level Python usage guide:
rechecked its CLI, installation, input, algorithm, and search-method statements against the current implementation and added a
complete `run_multiprocessing()` example showing `MPConfig`, completion-order results, the portable `spawn` guard, and the
separation between analysis and worker-process configuration.
320. Repaired the rendered command-line example:
placed the two matrix rows on separate Markdown lines so GitHub preserves the LaTeX row separator, displayed the exact ESS
with proper fractions, and replaced both prose value counts with properly rendered fractions instead of slash notation.
321. Clarified the CLI and Python result contracts:
documented that Python always returns weighted candidate and ESS counts plus support-size structures, while the CLI prints only the
ESS count by default and adds timing or representative CSV rows only when requested. Explained how CLI candidate columns can
reconstruct the structures without claiming that separate structure fields are emitted.

322. Clarified arbitrary-precision rational support:
explained that input entries and the canonical vector/payoff results are exact numerator/denominator pairs backed by
arbitrary-precision integers, without a fixed digit limit, while identifying `payoff_dbl` as an optional approximation.

323. Reviewed and tightened the top-level README:
kept it as an introduction and practical quick-start, reduced the algorithm description, qualified completeness as a `safe`
guarantee, defined the simplex indices, documented platform minima and the one-dimensional compact-input exception, made build
commands portable, and checked the CLI/Python result descriptions against their live outputs.

324. Removed duplicated search-method documentation:
deleted the redundant algorithm section and kept `safe` versus `fast` in one authoritative place.

325. Unified command-line and Python summary output:
made every CLI analysis emit one ten-field JSON summary on success or failure, removed the optional timing output mode, and retained
the candidate CSV only as an optional following table. Made JSON, CSV, and Parquet Python sinks use the same summary information,
while preserving nested structures in JSON and encoding them compactly in flat formats. The historical timing adapter still reads
old line-based binaries after detecting their `--timing` option.

326. Corrected the post-output-change documentation audit:
removed the last current-reference claim about three-line CLI timing output, documented that `--matrixid` also sets the JSON field,
made the build examples portable to macOS, and synchronized current matrix, calibration, fallback, candidate, and ESS counts with
the canonical SQLite database.

327. Published the FracESSA project website:
added a responsive static landing page using the existing logo, crawl metadata, `robots.txt`, and a sitemap; deployed it through a
dedicated GitHub Pages workflow; and connected the repository homepage to `https://reinhardullrich.github.io/fracessa/`.

328. Prepared release 2026.8.5.1:
updated the package homepage to the project website, made the README logo resolve on PyPI, and clarified that ordinary pushes can
deploy GitHub Pages but do not start release builds.

329. Planned direct exact stability reduction through the retained candidate factorization:
specified a narrow multi-right-hand-side solve for the already-factored negative-definite reduced Hessian and an integer-scaled Schur
complement that replaces complete rational Bee construction and recursive unrestricted-coordinate elimination. Recorded the exact
scaling proof, minimal source scope, proposed code, correctness regressions, and benchmark acceptance criteria without changing source.

330. Implemented direct exact stability reduction through the retained candidate factorization:
added a narrow multi-right-hand-side replay for an already-factored negative-definite reduced Hessian, constructed a positive integer
multiple of the smaller exact Schur complement, and removed complete rational Bee construction and recursive unrestricted-coordinate
elimination. A three-dimensional exact factor-reuse regression covers previous-pivot division and two simultaneous right-hand sides.
All C++ and Python tests passed, 7,808 canonical representative candidates representing 14,659 weighted candidates across 309
matrices matched exactly, and affected end-to-end benchmarks improved by about 9-46% without a measurable regression on the
unchanged early-rejection path. The final caller/callee cleanup passes the already-computed outside-best-reply mask into the Schur
builder and moves reference-index calculations used only by logging behind the logging condition.

331. Removed the obsolete copositivity-checker version suffix:
renamed the private `CopositivityCheckerV3` implementation class to `CopositivityChecker` without changing its algorithm or public
`is_strictly_copositive()` entry point, and updated the retained call-chain reference.

332. Corrected the reduced-Hessian rejection reason:
replaced the obsolete `kay_size`-dependent labels with `F_reduced_hessian_not_nd` whenever the exact reduced Hessian is not negative
definite. Exactly backfilled the canonical candidate database: 2,796 rows received the new label, while 344 genuine one-dimensional
Schur-complement rejections retained `F_not_pd_kay_0_1`. When comparing an older build, treat either historical label as equivalent
to the new label only when all other candidate fields match and the current build reports `F_reduced_hessian_not_nd`.

333. Aligned stability-matrix naming with Bomze's paper:
renamed the directly constructed scaled Schur result to the `scaled reduced B matrix`, matching the paper's final $B^{(r)}$ up to its
irrelevant positive scale. Source identifiers now use `build_scaled_reduced_b` and `scaled_reduced_b_`; the mathematical comments retain
that the implementation obtains this matrix through one exact Schur complement. No arithmetic or decision changed.

334. Clarified the exact stability decision flow:
added six short numbered path comments immediately above the mathematical branches in `check_stability()`. Each comment now states
the condition, its mathematical consequence, and whether the candidate is accepted, rejected, or passed to strict copositivity;
logging-only conditions remain uncommented. No executable code changed.

335. Removed a redundant reduced-$B$ alias:
made `check_stability()` use its accurately named `scaled_reduced_b_` member directly for logging, positive definiteness, and strict
copositivity instead of creating a second local reference with the same name. No matrix is copied and behavior is unchanged.

336. Standardized exact stability-path terminology:
called failed necessary conditions `early rejection`, sufficient positive results `early acceptance`, and the remaining
strict-copositivity branch the `final decision`. Removed the stale `Bee` wording from the matrix positive-definiteness comment. No
executable code changed.

337. Consolidated exact early stability decisions:
added a local research note that defines every current and proposed early acceptance, early rejection, exact reduction, and final
Hadeler optimization; records their proofs, costs, implementation status, recommended cheapest-first order, and the exact canonical
candidate audit supporting a diagonal-first experiment. Linked it from the current copositivity-flow note. No source code changed.

338. Added one shared exact copositivity sign scan:
scanned positive diagonals before one upper off-diagonal triangle, built fixed-size negative-neighbor masks, recorded row-level
positive and negative off-diagonal presence without another sign pass, and accumulated the denominator-one matrix's all-ones
quadratic value with integer additions. The stability path now rejects nonpositive diagonals and a nonpositive all-ones value before
exact positive-definiteness factorization while preserving the existing candidate status labels. Removed the conventional but
unnecessary positive factor two from the scaled reduced B construction while retaining its negation so later checks keep their
standard positive sign conventions.

339. Gave new exact early rejections certificate-specific stability reasons:
used `F_not_copos_nonpositive_diagonal` for a nonpositive diagonal witness and
`F_not_copos_nonpositive_all_ones_value` for a nonpositive all-ones quadratic witness. Historical stability labels and stored
candidate rows remain unchanged pending a separate deliberate database backfill.

340. Moved scaled reduced B logging to its construction point:
when logging is enabled, `check_stability()` now records `kay` and the complete scaled reduced B matrix immediately after building
it, before the sign scan and every early rejection. Numerical behavior is unchanged.

341. Clarified the all-ones early-rejection proof:
documented directly at Path 5 that the all-ones vector is nonzero and nonnegative, so a nonpositive exact quadratic value is an
explicit witness that the scaled reduced B matrix is not strictly copositive. No executable code changed.

342. Added exact early acceptance for nonnegative off-diagonal entries:
after positive diagonals pass, `check_stability()` now accepts immediately when the shared sign scan found no negative
off-diagonal entry, recording `T_copos_nonnegative_off_diagonal`. This is one mask comparison after the existing scan. The path also
covers the positive one-dimensional case, making the later `kay_size <= 1` branch and its popcount unreachable; both were removed.
Added an end-to-end regression for the new certificate-specific reason.

343. Shared direct exact copositivity criteria for dimensions one through three:
moved the low-dimensional Bomze/Hadeler arithmetic out of `check_stability()` into three reusable inline numerical functions.
`check_stability()` now only dispatches by outside-reply dimension and records `T_copos_k1` through `k3` or their matching
rejection labels. The existing Hadeler checker uses the same functions for principal subsets of sizes one through three without a
temporary matrix or LU factorization; its general determinant-adjugate path now begins at size four, so the obsolete special
one- and two-dimensional adjugate construction was removed. Before connecting the shared functions to Hadeler, exhaustive symmetric
integer matrices with entries from -2 through 2 produced identical decisions for all 15,755 cases. An end-to-end comparison of 165
stored games with actual small-K candidates found no mathematical or output difference beyond the intentional stability labels.

344. Made both low-dimensional dispatch sites explicit:
replaced the sequential `if` branches in `check_stability()` and the Hadeler principal-subset checker with `switch` statements over
their respective dimensions. Cases one through three still call the same shared exact criteria and `default` continues to the
general sign-scan or determinant-adjugate path. No arithmetic, output label, or decision changed.

345. Reduced the local early-stability research note to unfinished work:
removed every section describing checks and control flow already implemented in production, including obsolete implementation
plans and database evidence. The note now contains only the unimplemented Z-matrix shortcut, negative-graph decomposition, exact
row reductions, disconnected-subset shortcut, and one-solve Hadeler replacement, with the mathematics and experiment boundaries
needed to evaluate those proposals. Current behavior remains documented in the separate current-flow note.

346. Restored the complete cheapest-first stability diagram:
kept the implemented routing stages in one Mermaid diagram so the remaining exact experiments retain their proper algorithmic
context. Detailed prose in the early-decisions note still covers only unimplemented work.

347. Added exact symmetric Z-matrix rejection:
after exact positive definiteness fails, `check_stability()` now rejects a scaled reduced B matrix whose off-diagonal entries are all
nonpositive. The existing sign scan supplies this fact, so the runtime cost is one mask comparison and full Hadeler enumeration is
avoided. Added an end-to-end 4-by-4 reduced-matrix witness and the certificate-specific reason `F_not_copos_z_matrix`; moved the rule
from the remaining-experiments note into the documented current flow.

348. Corrected the documented reduced-matrix scale:
updated the current exact copositivity flow to reflect the already-implemented removal of Bomze's conventional factor two. The
stored matrix is $M=(d\delta/2)S$, not $M=d\delta S$; this remains a positive scaling and changes no mathematical decision.

349. Clarified cached reduced entries:
documented `reduced_entry()` as one integer-scaled reduced-Hessian entry in difference coordinates and explained why the support
Hessian and later Schur-complement blocks share its cache. No executable code changed.

350. Distinguished current and stored reduced dimensions:
renamed the persistent safe-solver member to `reduced_system_dimension_`. The local `reduced_dimension` remains the dimension derived
from the current support, while the member names the reusable reduced-system storage that must match it. No behavior changed.

351. Returned the outside-reply dimension from scaled reduced B construction:
made `build_scaled_reduced_b()` return the already-computed $|K|$ value and use it directly for the low-dimensional stability
dispatch. This replaces the later matrix-row lookup without another popcount or scan; numerical behavior is unchanged.

352. Numbered the low-dimensional exact stability path:
marked the $|K|=1,2,3$ switch as Path 4 and shifted the following stability comments to Paths 5 through 10. No executable code
changed.

353. Consolidated low-dimensional stability reasons:
replaced the six dimension-specific $|K|=1,2,3$ labels with `T_copos_small` and `F_not_copos_small`. The candidate already records
support and extended-support sizes, so the reason no longer duplicates the small dimension. No stored database row used the
temporary dimension-specific labels.

354. Added negative-part diagonal-dominance early acceptance:
the existing exact triangular sign scan now accumulates each row's diagonal plus its negative off-diagonal entries. When every sum
is strictly positive, `check_stability()` accepts immediately with `T_copos_negative_part_diagonal_dominance`, before the all-ones
witness and positive-definiteness factorization. The implementation cites Liqun Qi, *Linear Algebra and its Applications* 439
(2013), Theorem 10, equation (12), and uses one fixed stack buffer without another matrix scan. Added exact passing and equality-
boundary regressions.

355. Planned an integer-only final stability and copositivity path:
documented how to retain the already-integral scaled reduced B matrix in `fmpz_mat` storage, use FLINT's certified positive-
definiteness, determinant, solve, and nullspace operations, replace the nonsingular full inverse with one solve, replace the singular
cofactor construction with an exact nullspace decision, and verify both mathematical equivalence and measured speed before
retention. No source code changed.

356. Installed and verified project-local FLINT 3.6.0:
built the official release against Fedora's GMP and MPFR with native Release optimization, installed shared and static libraries
under `.local/flint-3.6.0`, and kept all source/build scratch under `.local-tmp`. The upstream `fmpz_mat` suite, a direct four-API
smoke test, the complete isolated FracESSA Release build, and all 10 C++/CLI tests passed; both the CLI and Python extension link to
the local FLINT library. Added both local directories to Git ignore and set FLINT 3.6.0 as the planned integer-path minimum. Fedora's
system FLINT 3.4.0 and production source remain unchanged.

357. Finalized the integer stability and copositivity migration plan:
rechecked every production and test caller, Hadeler's determinant/one-solve/nullspace decisions, and FLINT 3.6.0 contracts. The plan
now defines the complete minimal source and deletion scope, reusable scratch objects, true dimension-four and arbitrary-precision
regressions, rational-versus-integer output equivalence, isolated benchmarking, the configure-time version gate, and implementation
order. Rational arithmetic remains at exact input and candidate-output boundaries; no source code changed.

358. Corrected plan mathematics rendering:
replaced fenced `math` blocks in the integer stability and copositivity plan with standard display-math delimiters so the formulas
render consistently in GitHub and Codex. The formulas themselves are unchanged.

359. Added the mathematical proof behind the integer Hadeler algorithm:
traced each of its six steps to Hadeler's definition and Theorem 3, and distinguished the paper's published equivalences from
FracESSA's derived one-solve and nullspace implementations. The proof explains why only principal submatrices are enumerated, why
cardinality order is essential, and why general minors, cofactors, the inverse, and the explicit adjugate can be removed.

360. Cross-checked the integer Hadeler proof against independent literature:
linked the earlier Cottle-Habetler-Lemke criterion, Väliaho's singular almost-strict classification, Kaplan's
principal-submatrix eigenvalue test, the Hiriart-Urruty-Seeger survey, and the 2020 *Matrix Positivity* restatement. The plan now
marks the fixed-right-hand-side solve and explicit-adjugate removal as FracESSA derivations rather than quotations from a source.

361. Separated the KKT and general fraction-free LDLT implementations:
renamed the existing hot-path file and class to `fraction_free_ldlt_kkt.hpp` and `kkt_fraction_free_ldlt_workspace` without changing
its arithmetic. Added `fraction_free_ldlt_factorization`, which factors a symmetric integer matrix once, exposes its exact signed
determinant, and solves multiple later right-hand sides from the retained lower triangle. Both factorization and solve retain the
FLINT-style immediate-integer path and switch to arbitrary-precision `fmpz` arithmetic when necessary. Focused tests cover ordinary
and arbitrary-precision solves, singular determinants, zero-diagonal coordinate operations, and 280 deterministic symmetric
matrices checked against FLINT. The new general factorization is not yet connected to production stability checking.

362. Moved the complete final stability matrix path from denominator-one rationals to exact integers:
`build_scaled_reduced_b()` now stores its already-integral entries directly in `matrix_int`; the small-dimensional criteria, shared
sign scan, positive-definiteness decision, and Hadeler checker all operate on `fmpz` values. The general fraction-free $LDL^T$
factorization now records exact inertia, supplies the positive-definiteness certificate, provides determinants for symmetric
principal matrices, and solves all identity columns for the unchanged nonsingular adjugate calculation. Singular off-diagonal
cofactors retain the same expensive algorithm but use FLINT integer determinants because those minors are not symmetric. No path
order, reason string, support order, candidate arithmetic, public rational output, or Hadeler decision changed. Genuine
four-dimensional positive, negative-determinant, singular, and arbitrary-precision tests were added. All ten C++/CLI tests pass with
FLINT 3.6 and system FLINT 3.4, all 66 Python tests and 49 subtests pass, and 883 old/new safe-mode outputs match byte-for-byte after
excluding only timing. Pinned end-to-end new/old median ratios were 0.292, 0.986, 1.003, and 1.006 on representative dimensions 14,
19 circular, 23, and 24, so the migration has no material measured regression and greatly improves the copositivity-heavy case.

363. Added a dedicated exact reduced-B test corpus:
created `testdata/Copos_testdata.sqlite3` with 1,069 permutation-inequivalent scaled reduced B matrices. Every row stores the
symmetric integer upper triangle, strict-copositivity result, and source matrix/candidate/support provenance. Extraction replayed only
exact or exact-fallback stored candidates through the current exact solver, then collapsed 6,691 repeated occurrences and 1,411
distinct orientations related by simultaneous row-and-column permutations. Unverified fast-only baselines and three affine image
rows that production does not solve separately were excluded. No benchmark schema or maintained extraction framework was added.

364. Replaced explicit Hadeler adjugates with exact minimal witnesses:
for negative principal determinants, reused the retained fraction-free $LDL^T$ factorization to solve only
$Cy=-\mathbf 1$ and reject exactly when $y>0$. For zero determinants, replaced all cofactor determinants with one exact FLINT
nullspace calculation and the one-sign basis test. Removed the now-unused adjugate helpers and added passing, rejecting, singular,
higher-nullity, and arbitrary-precision regressions. All 1,069 exact reduced-B corpus decisions and seven complete safe candidate
outputs matched. All FLINT 3.6 CTest/Python checks and the focused FLINT 3.4 test passed. Nine alternating, uncontended CPU-2 runs
improved the corpus wall-time median from 4.660 s to 0.995 s (4.68x, 78.65% less time); process CPU time independently produced the
same ratio. The implementation has no material end-to-end regression on the larger source games.

365. Simplified and extended the exact copositivity corpus:
replaced all duplicated candidate/support provenance in `Copos_testdata.sqlite3` with one nullable `fracessa_matrix_id`. Existing
matrix IDs, exact upper triangles, classifications, dates, and FracESSA matrix links remain unchanged. Added nine
permutation-inequivalent references at IDs 9157-9165: published matrices M1-M7 from Bras, Eichfelder, and Judice plus strict and
non-copositive sides of their C5 construction. The former adjugate checker and current solve/nullspace checker both matched all nine
expected strict-copositivity results. The migrated 1,078-row database passes SQLite integrity and upper-triangle shape checks.

366. Removed the retired exact-rational LU and rational positive-definiteness implementations after the integer stability migration
left both without production callers. Deleted the LU test executable and CMake registration, the rational PD unit case, and the
fraction/matrix operations used only by those implementations and their self-focused tests. Rational matrices now remain only at
the parsed-game and candidate-output boundary.

367. Replaced only the final Hadeler principal-submatrix fallback with an exact adaptive simplicial-cone checker while preserving
every preceding small-matrix, sign, diagonal-dominance, all-ones, positive-definiteness, and Z-matrix decision. The independently
written FLINT implementation follows the pair-splitting idea in Mathieu Dutour Sikiric's Polyhedral Common, stores only
$B=VAV^T$, and updates one exact row and column per new generator. Both production and isolated implementations matched all 1,078
stored copositivity matrices, and the complete C++ suite passed. A full-corpus one-shot screen found 21 final-branch matrices: the
cone's median Hadeler/cone speed ratio was 1.328x overall, 1.155x for 14 rejected matrices, and 167.38x for seven strictly copositive
matrices. The previous repeated benchmark independently showed the same qualitative result: optimized Hadeler remains competitive
for quick witnesses but becomes impractical on strictly copositive matrices because it must traverse principal subsets.

368. Removed the final-matrix positive-definiteness and symmetric Z-matrix shortcuts after a complete ablation:
all 1,078 copositivity matrices completed correctly in both variants with no timeout. Of those, 1,057 exited through earlier cheap
routes and 21 reached the compared stage; neither positive definiteness nor the Z-matrix condition decided any of the 21. Skipping
both checks reduced the relevant matrices' median time by 41.82% (1.719x speedup), with median reductions of 12.19% for the seven
strict matrices and 49.45% for the fourteen rejected matrices. The exact cone now handles every unresolved matrix directly. Removed
the checker-owned factorization and positive-sign masks used only by the Z-matrix shortcut; retained the negative-neighbor masks for
the planned connected-component reduction. Archived the now-unreferenced general fraction-free LDLT implementation under `archive/`.
The separate KKT-specialized fraction-free LDLT remains active for candidate solving and reduced-Hessian inertia.

369. Added the exact negative-entry connected-component reduction before the final cone check. The existing sign scan now supplies
the complete adjacency masks without another matrix pass. A connected graph enters the cone checker without copying; a disconnected
graph is checked through its non-singleton principal component matrices, while positive singleton diagonals need no further work.
All 30 focused copositivity tests, all nine C++ test executables, and all 66 Python tests passed. Complete safe output matched the
pre-change binary on all 768 stored matrices through dimension 25 with sub-second calibrations. One CPU-2 measurement of each
81-matrix quick case was neutral at the median: -0.24% elapsed time in fast and -0.07% in safe; this screen primarily measures the
connected or earlier-exit overhead. A focused CPU-2 benchmark screened all 1,078 reduced-B matrices and timed only the 21 that reach
the final checker, using a 0.5-second native-timing target per matrix and method. The seven strictly copositive matrices improved by
15.66% on average and 17.06% at the median (1.206x median speedup). The fourteen non-copositive matrices were neutral at the median
and 5.91% slower on average because their direct cone checks usually find a witness in only a few microseconds.
