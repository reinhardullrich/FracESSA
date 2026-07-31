# C++ Review

Last verified: 2026-07-31

Scope: active C++ analyzer core, CLI, shared parser, CMake, C++/CTest coverage,
and the release workflow. The native Python binding is reviewed separately in
`PYBIND_REVIEW.md`. Frozen experiment source copies are excluded except where a
dated result is cited as evidence.

Correctness is ranked before speed. This file records unresolved findings first,
then explicit reassessment decisions so closed or conditional findings do not
silently reappear in later reviews. Remove an unresolved finding after its fix
and regression coverage are complete.

## Speed

### P2: Partial copositivity retains every reduction matrix

`check_stability()` allocates full histories for masks, sizes, and matrices at
`cpp/src/checkstab.cpp:123`. Each iteration needs only the previous and current
matrices plus current scalar mask state.

Required outcome: use two rolling matrix buffers and scalar current-state
variables; each step can be logged before the buffers are exchanged.

### P2: Exact LU uses a dense rational permutation matrix

`LU_Factorization` builds `m_P` as an `n*n` exact identity matrix at
`cpp/include/linalg/lu_factor_fraction.hpp:30`. Every solve then applies it with
a dense rational matrix-vector multiply at
`cpp/include/linalg/lu_factor_fraction.hpp:109`, performing mostly zero
multiply-adds.

Required outcome: store one permutation-index array and apply it by direct
indexing.

## Build And Release

### P2: C++ tests cannot be disabled

`cpp/CMakeLists.txt:45` declares all four FetchContent projects unconditionally,
`cpp/CMakeLists.txt:74` fetches them together, and `cpp/CMakeLists.txt:210` always
adds tests. `BUILD_TESTING=OFF` is not wired, so every build fetches and builds
GoogleTest and the test targets.

Required outcome: use standard `include(CTest)`/`BUILD_TESTING` and gate
GoogleTest plus `add_subdirectory(tests)`. Add a separate Python-module option
only if a CLI-only build becomes an actual supported workflow.

### P2: Debug configurations are forcibly optimized with assertions disabled

Global compile options at `cpp/CMakeLists.txt:21` and
`cpp/CMakeLists.txt:28` apply `/O2` or `-O3` and `NDEBUG` regardless of the
selected configuration. A nominal Debug build therefore still disables
assertions and compiles production optimization flags.

Required outcome: let CMake's standard build-type flags provide optimization
and `NDEBUG`; scope only FracESSA-specific throughput flags to Release builds.

## Ponytail Simplification Findings

`cpp/include/fracessa/supports.hpp:236-371`: delete: the 136-line experimental
`CircularSupportGeneratorV2` has no caller or test. Delete it and its private
`cpp/include/fracessa/bitset64.hpp:75-82` `rot_left()` helper; nothing replaces
them.

`cpp/include/fracessa/bitset64.hpp:148-165`: delete:
`is_smallest_representation()` is used only by its own test at
`cpp/tests/test_bitset64.cpp:161-168`. Delete both; nothing replaces them.

`cpp/include/linalg/matrix_double.hpp:26` and
`cpp/include/linalg/matrix_fraction.hpp:51`: delete: the mutable matrix `data()`
accessors have no caller. Nothing replaces them.

`cpp/include/fracessa/bitset64.hpp:105-112`: delete:
`find_pos_next_set_bit()` has no production caller after copositivity switched
to one fixed index array. Delete its self-focused tests at
`cpp/tests/test_bitset64.cpp:44-108` and `:222-229`; production set iteration
already uses `extract_set_indices()`.

Net: approximately 250 production and self-testing lines can be deleted.

## Reassessed Non-Findings

- Keeping `find_candidate_exact`, `find_candidate_verified`, and
  `find_candidate_unsafe` as three concrete classes is justified. Each owns
  different reusable scratch state for a distinct algorithm; another shared
  interface or a return to one large `fracessa` implementation would add
  coupling rather than remove it.
- The proof-helper declarations in `find_candidate_verified.hpp` are justified
  by focused tests of the rounding, factorization, residual, error-bound, and
  rejection steps. Splitting the 727-line proof kernel into more classes or
  files would not make the mathematical flow smaller.
- Moving MatrixParser's unchanged implementation into
  `cpp/include/fracessa/matrix_parser.hpp` is behavior-preserving. Its symbols
  have the required inline linkage, and the parser and CLI tests pass.
- Linux and macOS release executables still dynamically depend on FLINT, MPFR,
  and GMP, but `README.md` now documents those dependencies and explicitly says
  the binaries are not universal standalone artifacts. Packaging is therefore
  not a current defect. Reopen it only if portable download-and-run binaries
  become an actual release goal.
- A complete SQLite verification sweep in every ordinary CI run is not required.
  The focused C++/CLI suite is the intended fast gate; long database validation
  remains an explicit manual check.
- The repeated copositivity index scan is closed. Each subset is now extracted
  once into a fixed stack array. Three alternating CPU-2 Release/LTO runs, each
  using the median of 21 hot samples on the strictly copositive, non-positive-
  definite matrix with diagonal 1 and off-diagonal 2, changed by -0.11%,
  +0.19%, +0.37%, and +0.24% at dimensions 4, 6, 8, and 10, respectively. That
  is neutral within measurement noise; the shorter loop is retained without
  claiming a speedup.

## Current Validation State

- A warning-enabled Release build using `-Wall`, `-Wextra`, `-Wpedantic`,
  `-Wconversion`, and `-Wshadow` passed all 11 C++/CLI tests. It found no new
  FracESSA production warning; the remaining diagnostics are test-only
  conversions/bracing and bundled spdlog/fmt diagnostics.
- Streaming Gosper enumeration, the 62-by-62 immediate-rejection regression,
  the singular Hadeler rejection, and the late-pivot 3-by-3 LU solve/inverse
  regression pass all 11 C++/CLI tests. All 54 wrapper tests also pass.
- Exact candidate search now validates outside-support strategies before dense
  vector construction and skips that construction when neither output nor
  logging needs it; the focused requested/successful-vector regression passes.
- The latest complete verified-mode sweep matched the stored ESS count for all
  87 retained SQLite matrices. That long sweep and ASan/UBSan were not rerun for
  the current local source and regression-test changes.
- Wrapper tests: see `PYBIND_REVIEW.md` and `PYTHON_REVIEW.md`.
- The single parser preserves 18-digit direct values and arbitrary-precision
  values, rejects dimensions outside 1-63, and reports failures through
  `std::invalid_argument`.
- Production DFS/FKM generators retain no complete support frontier or
  cardinality layer. An independent order-insensitive comparison matched all
  mathematical candidate rows and ESS results across the former 52-matrix
  verification corpus.
- The canonical SQLite snapshot stores 49,157 candidate representatives whose
  multipliers recover 86,152 candidates and 83,377 ESS across 87 matrices.
- A fixed-seed audit generated 20,000 exact 4-by-4 integer matrices; all 19,890
  nonsingular cases satisfied `A * inverse(A) == I` exactly.
- ASan/UBSan passed all 11 C++/CLI tests on the combined tree.
- Ordinary pushes and pull requests run the same three-platform build and fast
  test matrix as tags; artifact packaging and publication remain tag-only.
- One wrapper integration regression exercises database ID 46 through verified
  and unsafe modes; no complete SQLite matrix-verification runner is wired into
  CTest or CI.
