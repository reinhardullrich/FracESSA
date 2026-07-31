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

## Correctness And Resource Safety

### P1: Copositivity materializes all subsets before early rejection

`CopositivityCheckerV3::is_strictly_copositive()` builds every nonempty subset
at `cpp/include/linalg/copositive_fraction.hpp:150` before testing the first
principal submatrix at `cpp/include/linalg/copositive_fraction.hpp:160`.

This is not merely a speed problem. A supported dimension-63 game can produce a
62-by-62 copositivity matrix whose first diagonal entry is already negative. A
streaming implementation would reject its first 1-by-1 principal submatrix
immediately, whereas the current implementation first attempts to reserve and
store up to `2^62 - 1` masks. The resulting allocation failure makes a valid
input impossible to analyze.

The retained `supports_` member also makes a checker unsafe to reuse directly:
`resize(n)` does not clear surviving inner vectors, so a later call can repeat
old masks or retain masks containing indices outside a smaller matrix. The
production wrapper currently constructs a fresh checker for every call, but
streaming removes this latent state problem as well.

Required outcome: generate fixed-cardinality masks one at a time and stop at
the first rejection. Do not add another generic support-generator abstraction;
the small local mask loop is sufficient. Remove `supports_` and the then-unused
`binomial_coefficient()` helper.

## Speed

### P2: The exact candidate path materializes its full vector too early

The exact path fills `result.vector` at `cpp/src/find_candidate_exact.cpp:65`
before outside-support validation beginning at
`cpp/src/find_candidate_exact.cpp:75`. It also does
this when candidate output and logging are both disabled; stability itself does
not consume the vector.

Required outcome: validate first and materialize a full vector only for a
successful candidate that will be output or logged.

### P2: Copositivity rescans set indices for every row

Candidate search now extracts each stage's support partition once and reuses it
for matrix construction and validation. Bee construction uses the later
`extended_support_reduced`, so it cannot share those earlier partitions. Within
`is_copositive_hadeler()`, however, the same subset mask is rescanned from its
first bit for every row at
`cpp/include/linalg/copositive_fraction.hpp:106` and
`cpp/include/linalg/copositive_fraction.hpp:108`. One fixed index array removes
those nested scans. The direct bit-scanning locations are listed in
`../reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`.

Required outcome: extract the current subset once before its nested matrix
loops and reuse the fixed index array. Benchmark before retaining the change.

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

## Test Coverage

### P2: Two proof-critical exact branches lack focused regressions

All durable LU tests in `cpp/tests/test_lu.cpp` use 2-by-2 matrices. None forces
a row swap after an earlier elimination step, which is the only case that must
swap already-computed columns of `L` at
`cpp/include/linalg/lu_factor_fraction.hpp:55`. An independent randomized exact
audit passed, but it is not a repository regression.

The copositivity suite exercises singular acceptance with the all-ones matrix,
but not the singular Hadeler rejection at
`cpp/include/linalg/copositive_fraction.hpp:127` where `det(A) == 0` and
`adj(A) > 0` entrywise.

Required outcome: add the nonsingular matrix `[[1,1,0],[1,1,1],[0,1,1]]` as a
3-by-3 LU solve/inverse case; it pivots after column zero. Add the exact matrix
`[[1,-1],[-1,1]]` as a strict-copositivity rejection. These two small tests
cover the branches directly.

## Mathematical Documentation

### P3: The singular Hadeler explanation contradicts the implemented branch

The comment at `cpp/include/linalg/copositive_fraction.hpp:119` says that under
the theorem's proper-principal-submatrix precondition, a singular matrix cannot
have the positive-adjugate rejecting pattern. The proposed regression matrix

```text
[[1, -1],
 [-1, 1]]
```

is a direct counterexample to that sentence: both proper principal submatrices
are positive, its determinant is zero, and its adjugate is entrywise positive.
The implementation correctly rejects it because the quadratic form vanishes at
`(1,1)`; only the explanation is wrong.

Required outcome: correct the comment when adding the focused singular Hadeler
regression above. Do not change the rejecting branch.

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

`cpp/include/fracessa/bitset64.hpp:140-157`: delete:
`is_smallest_representation()` is used only by its own test at
`cpp/tests/test_bitset64.cpp:161-168`. Delete both; nothing replaces them.

`cpp/include/linalg/matrix_double.hpp:26` and
`cpp/include/linalg/matrix_fraction.hpp:51`: delete: the mutable matrix `data()`
accessors have no caller. Nothing replaces them.

Net: approximately 172 production and self-testing lines can be deleted.

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

## Current Validation State

- A warning-enabled Release build using `-Wall`, `-Wextra`, `-Wpedantic`,
  `-Wconversion`, and `-Wshadow` passed all 11 C++/CLI tests. It found no new
  FracESSA production warning; the remaining diagnostics are test-only
  conversions/bracing and bundled spdlog/fmt diagnostics.
- The latest complete verified-mode sweep matched the stored ESS count for all
  87 retained SQLite matrices. That long sweep and ASan/UBSan were not rerun for
  this reassessment because the current C++ working-tree change only relocates
  MatrixParser into its header.
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
