# C++ Review

Last verified: 2026-07-29

Scope: active C++ analyzer core, CLI, shared parser, CMake, C++/CTest coverage,
and the release workflow. The native Python binding is reviewed separately in
`PYBIND_REVIEW.md`. Frozen experiment source copies are excluded except where a
dated result is cited as evidence.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Correctness

### P0: The double candidate filter drops valid ESS supports

The temporary default filter exactly normalizes the game and sends small-pivot,
non-finite, and danger-veto cases to exact arithmetic. However,
`is_suspicious()` still uses `minimum_pivot * margin` only as a heuristic proxy
for forward solution error at `cpp/src/findeq.cpp:153`. The negative-probability
and outside-gain decisions at `cpp/src/findeq.cpp:163` and
`cpp/src/findeq.cpp:183` can therefore reject a valid support before exact
arithmetic sees it.

Preserved verification IDs 45-47 at reference commit `2be0207` each have one
exact full-support ESS, while this unsafe filter returns zero candidates and
zero ESS. ID 45 is especially important: its entries are approximately
`-139.66` to `80.15` and its exact minimum support probability is `1e-6`, yet
the normalized bordered system has condition number about `5.34e11`; the
minimum-pivot danger test accepts the wrong negative probability as decisive.

Required outcome: make the later rigorously one-sided Choice 1 filter the
default and retain this heuristic only as explicit `--unsafe` behavior. The
implementation plan is in `../architecture/CHOICE_ONE_CANDIDATE_FILTER.md`.

## Speed

### P1: Normal search materializes the complete exponential support frontier

`Supports::initialize()` stores every nonempty support in cardinality buckets at
`cpp/include/fracessa/supports.hpp:61` before the first normal support is solved
at `cpp/src/fracessa.cpp:67`. For a non-circular game this is `2^n - 1` masks,
in addition to the unavoidable search time. At `n=30`, the masks alone require
almost 8 GiB before vector overhead and before pruning removes any supersets.

A direct current-code measurement at `n=23` stored 8,388,607 masks and reached
about 69 MiB maximum RSS. The circular path stores fewer coprime-cardinality
representatives, but still scans every mask and tests up to `n` rotations at
`cpp/include/fracessa/supports.hpp:72`.

The isolated experiment in `../experiments/SUPPORT_FRONTIER_2026-07-29.md`
preserved byte-identical output on 44 official and 35 generated games. A
byte-indexed fixed-cardinality Gosper generator improved strongly pruned cases
by approximately 30-74% and cut matrix-34 peak RSS roughly in half, at the cost
of approximately 2% median slowdown on weak random games. Recursive DFS was not
predictable from dimension and the larger adaptive implementation was rejected
as unnecessary complexity.

A follow-up on the initially omitted prime circular dimension-19 games found a
52.07% regression on ID 28. Its 1,843 rotated size-7/8 candidates repeatedly
mark approximately 6.73 million supersets, although circular reduction stores
only one representative per 19-mask orbit. ID 29 regressed only 2.30%, proving
that neither dimension nor primality alone is a sufficient gate.

The complete 44-matrix pass found all five major 30-73% wins at composite
circular dimensions 18, 20, 21, 22, and 24. Prime circular IDs 22, 26, 28, and
33 instead regressed by 2.8-52.4%. Candidate output remained byte-identical in
every case.

Required outcome: generate supports by cardinality on demand while preserving
the existing numeric order and exact-equilibrium superset pruning. Use the
byte-indexed Gosper experiment only as a starting point, eliminate its repeated
direct marking for circular representatives without adding a large framework,
and rebenchmark every official dimension 18-24 matrix before retention.

### P1: The exact candidate path materializes its full vector too early

The exact path fills `candidate_.vector` at `cpp/src/findeq.cpp:211` before the
outside-support validation beginning at `cpp/src/findeq.cpp:222`. It also does
this when candidate output and logging are both disabled; stability itself does
not consume the vector.

Required outcome: validate first and materialize a full vector only for a
successful candidate that will be output or logged.

### P2: Set indices are rescanned at multiple exact stages

The unsafe filter extracts support and complement once. When exact arithmetic is
required, the support is extracted again at
`cpp/include/fracessa/matrix_server.hpp:91` while building the bordered system
and both support and complement are rebuilt at `cpp/src/findeq.cpp:206`.

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
`cpp/src/main.cpp:63` and `cpp/src/main.cpp:65`. The C++ standard does not require
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

`cpp/CMakeLists.txt:44` declares all four FetchContent projects unconditionally,
`cpp/CMakeLists.txt:73` fetches them together, and `cpp/CMakeLists.txt:181` always
adds tests. `BUILD_TESTING=OFF` is not wired, so every build fetches and builds
GoogleTest and the test targets.

Required outcome: use standard `include(CTest)`/`BUILD_TESTING` and gate
GoogleTest plus `add_subdirectory(tests)`. Add a separate Python-module option
only if a CLI-only build becomes an actual supported workflow.

### P2: Debug configurations are forcibly optimized with assertions disabled

Global compile options at `cpp/CMakeLists.txt:20` and
`cpp/CMakeLists.txt:27` apply `/O2` or `-O3` and `NDEBUG` regardless of the
selected configuration. A nominal Debug build therefore still disables
assertions and compiles production optimization flags.

Required outcome: let CMake's standard build-type flags provide optimization
and `NDEBUG`; scope only FracESSA-specific throughput flags to Release builds.

### P2: Linux and macOS release executables are not self-contained

CMake prefers `.so`/`.dylib` over static archives at `cpp/CMakeLists.txt:82` and
`cpp/CMakeLists.txt:85`, and the workflow uploads only the executable at
`.github/workflows/release.yml:125`. The runner installs FLINT, MPFR, and GMP,
but an end user still needs ABI-compatible libraries; the macOS binary also
records a Homebrew library path.

Required outcome: either publish a documented system-dependency package or
bundle/link the mathematical runtime consistently. Do not describe the current
Linux and macOS artifacts as portable standalone executables.

## Current Validation State

- Current local Release build: successful.
- Core/CLI CTests: 10/10 passed.
- Wrapper tests: 23/23 passed.
- Full matrix CTests: 44/44 passed, including normalization IDs 38-39.
- The single parser preserves 18-digit direct values and arbitrary-precision
  values, and the former signed-overflow input passes ASan/UBSan exactly.
- All active fast and full verification matrices pass. Reference IDs 45-47
  remain outside active fixtures because unsafe mode is known to fail them.
- The current standard library reports `high_resolution_clock::is_steady ==
  false` and `steady_clock::is_steady == true`.
- A temporary `n=23` support-generation check stored 8,388,607 masks and used
  about 69 MiB maximum RSS before any analyzer work.
- A fixed-seed audit generated 20,000 exact 4-by-4 integer matrices; all 19,890
  nonsingular cases satisfied `A * inverse(A) == I` exactly.
- The most recent full sanitizer suite passed its existing tests.
- The current implementation is uncommitted on `choice-one-candidate-filter`,
  based on `main` at `32f6167`; no GitHub Actions run exists for this diff.
