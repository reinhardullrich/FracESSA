# C++ Review

Last verified: 2026-07-30

Scope: active C++ analyzer core, CLI, shared parser, CMake, C++/CTest coverage,
and the release workflow. The native Python binding is reviewed separately in
`PYBIND_REVIEW.md`. Frozen experiment source copies are excluded except where a
dated result is cited as evidence.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Speed

### P1: The exact candidate path materializes its full vector too early

The exact path fills `candidate_.vector` at `cpp/src/findeq.cpp:227` before the
outside-support validation beginning at `cpp/src/findeq.cpp:238`. It also does
this when candidate output and logging are both disabled; stability itself does
not consume the vector.

Required outcome: validate first and materialize a full vector only for a
successful candidate that will be output or logged.

### P2: Set indices are rescanned at multiple exact stages

The unsafe filter extracts support and complement once. When exact arithmetic is
required, the support is extracted again at
`cpp/include/fracessa/matrix_server.hpp:95` while building the bordered system
and both support and complement are rebuilt at `cpp/src/findeq.cpp:222`.

Bee construction uses the later `extended_support_reduced`, so it cannot share
the initial extraction. Within `is_copositive_hadeler()`, however, the same
subset mask is rescanned from its first bit for every row at
`cpp/include/linalg/copositive_fraction.hpp:106` and
`cpp/include/linalg/copositive_fraction.hpp:108`. One fixed index array removes
those nested scans. The direct bit-scanning locations are listed in
`../reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`.

Required outcome: keep logically different sets as separate extraction stages,
but pass one exact-stage partition through bordered-system construction and
candidate validation. Benchmark the interface change before retaining it.

### P2: Partial copositivity retains every reduction matrix

`check_stability()` allocates full histories for masks, sizes, and matrices at
`cpp/src/checkstab.cpp:98`. Each iteration needs only the previous and current
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

### P2: Copositivity materializes all subsets before early rejection

`CopositivityCheckerV3::is_strictly_copositive()` builds every nonempty subset
at `cpp/include/linalg/copositive_fraction.hpp:150` before testing the first
principal submatrix at `cpp/include/linalg/copositive_fraction.hpp:160`.

Required outcome: generate fixed-cardinality masks on demand and stop at the
first failure. Benchmark any destructive or move-based LU variant separately.

### P2: CLI elapsed time uses a non-monotonic clock

The CLI measures analyzer duration with `std::chrono::high_resolution_clock` at
`cpp/src/main.cpp:56` and `cpp/src/main.cpp:58`. The C++ standard does not require
that clock to be steady; on the current libstdc++ build it aliases a clock with
`is_steady == false`. A wall-clock correction can therefore distort or even
reverse an elapsed interval used for speed baselines.

Required outcome: use `std::chrono::steady_clock` for elapsed CLI timing. This is
a direct standard-library replacement with no new mechanism.

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

## Build And Release

### P1: CI runs only for release tags

`.github/workflows/release.yml:3` has only a `v*` tag trigger. Pull requests and
ordinary pushes to `main` can therefore merge correctness failures without any
GitHub build or test signal.

Required outcome: run build plus fast tests on pushes and pull requests; keep
artifact publication restricted to release tags.

### P2: C++ tests cannot be disabled

`cpp/CMakeLists.txt:45` declares all four FetchContent projects unconditionally,
`cpp/CMakeLists.txt:74` fetches them together, and `cpp/CMakeLists.txt:212` always
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

### P2: Linux and macOS release executables are not self-contained

CMake prefers `.so`/`.dylib` over static archives at `cpp/CMakeLists.txt:85` and
`cpp/CMakeLists.txt:87`, and the workflow uploads only the executable at
`.github/workflows/release.yml:130`. The runner installs FLINT, MPFR, and GMP,
but an end user still needs ABI-compatible libraries; the macOS binary also
records a Homebrew library path.

Required outcome: either publish a documented system-dependency package or
bundle/link the mathematical runtime consistently. Do not describe the current
Linux and macOS artifacts as portable standalone executables.

## Current Validation State

- Current local Release build: successful.
- Core/CLI CTests: 11/11 passed.
- Wrapper tests: 26/26 passed.
- Full matrix CTests: 55/55 passed, including candidate-rejector-double IDs 45-47.
- The single parser preserves 18-digit direct values and arbitrary-precision
  values, and the former signed-overflow input passes ASan/UBSan exactly.
- The exact-rejection oracle passed all 53 permitted matrices; IDs 33-34 were
  excluded from that oracle run as required. ASan/UBSan passed all 11 core/CLI
  tests and focused verification IDs 45-47.
- The current standard library reports `high_resolution_clock::is_steady ==
  false` and `steady_clock::is_steady == true`.
- Production DFS/FKM generators retain no complete support frontier or
  cardinality layer. An independent order-insensitive comparison matched all
  mathematical candidate rows and ESS results across the original 52-matrix
  generator suite; the current 55-matrix suite also passes.
- A fixed-seed audit generated 20,000 exact 4-by-4 integer matrices; all 19,890
  nonsingular cases satisfied `A * inverse(A) == I` exactly.
- The most recent full sanitizer suite passed its existing tests.
- The current implementation is uncommitted on `certified-filter`;
  no GitHub Actions run exists for this diff.
