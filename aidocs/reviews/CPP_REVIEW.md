# C++ Review

Last verified: 2026-07-27

Scope: active C++ analyzer core, CLI, shared parser, CMake, C++/CTest coverage,
and the release workflow. The native Python binding is reviewed separately in
`PYBIND_REVIEW.md`. Frozen experiment source copies are excluded except where a
dated result is cited as evidence.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Correctness

### P0: The double candidate filter drops valid ESS supports

`solve_linear_dbl()` rejects pivots using fixed absolute `1e-12` thresholds at
`cpp/include/linalg/linear_solver.hpp:47` and
`cpp/include/linalg/linear_solver.hpp:63`. `find_candidate_dbl()` treats every
such numerical rejection as conclusive at `cpp/src/findeq.cpp:42`, so exact
arithmetic never sees the support.

The double solver can also reject a rounded support probability below `-1e-10`
at `cpp/include/linalg/linear_solver.hpp:77`. This is another conclusive
floating sign decision on the exact candidate path.

The outside-support test is also conclusive and uses the fixed absolute margin
`1e-4 * dimension` at `cpp/src/findeq.cpp:63`, with rejection at
`cpp/src/findeq.cpp:72`. It therefore has the same scale-dependence problem even
when the linear solve succeeds.

Verification IDs 38 and 39 both have the mixed ESS `(1/2,1/2)`. The default
path returns zero ESS, while `--exact` returns one. Scaling or translating the
same game therefore changes the program's answer. Both cases exercise the
double solve: ID 38 falls below the absolute pivot threshold, while ID 39 loses
the exact `10^20` versus `10^20+1` distinction during conversion to double.
They do not exercise the outside-support loop because their valid candidate has
full support; that unsafe comparison does not yet have a dedicated regression.

Required outcome: remove the double candidate filter from the decision path and
benchmark the exact-only path. Falling through to exact arithmetic after every
double rejection saves no exact work, so keeping the current filter would only
add work. Reintroduce a numerical filter only if a rigorously one-sided method
can reject supports without false negatives.

Design note: `../correctness/CERTIFIED_CANDIDATE_FILTER.md` records exact affine
normalization and a possible rigorous one-sided interval/error-bound filter.

## Speed

### P1: Candidate paths materialize full vectors before they are needed

The double path fills `full_solution[64]` at `cpp/src/findeq.cpp:50`, then reads
only support entries that already exist as `solution(j_pos, 0)`. The entire
array can be deleted.

The exact path fills `candidate_.vector` at `cpp/src/findeq.cpp:97` before
outside-support validation and even when candidate output and logging are both
disabled. Stability does not consume the vector.

Required outcome: use the compact double solution directly; in the exact path,
validate first and materialize a full vector only for a successful candidate
that will be observed.

### P2: Solvers allocate output for every nonsingular back-substitution

`solve_linear_dbl()` and `solve_linear_frc()` assign a freshly constructed
output matrix at `cpp/include/linalg/linear_solver.hpp:66` and
`cpp/include/linalg/linear_solver.hpp:139`. This happens after elimination, so
early singular rejections do not allocate, but every solve that reaches back
substitution does. The exact allocation also constructs a vector of FLINT
objects.

The solved values are copied before the next support and do not escape by
reference. Required outcome: benchmark pre-sized scratch outputs, then make the
solvers overwrite them if the allocation cost is measurable.

### P2: Set indices are rescanned at multiple stages

The original support is extracted independently while building the double and
exact bordered systems and again in each candidate check. Extract it and its
complement once for that support and pass the fixed stack buffers through those
operations.

Bee construction uses the later `extended_support_reduced`, so it cannot share
the initial extraction. Extract that derived set only when it becomes known.
Within `is_copositive_hadeler()`, the same subset mask is rescanned for every
row at `cpp/include/linalg/copositive_fraction.hpp:111`; one fixed index array
would remove the nested bit scans. The direct bit-scanning locations are listed
in `../reference/FIND_POS_FIRST_SET_BIT_CALL_CHAIN.md`.

Required outcome: keep these as separate extraction stages and benchmark the
interface changes before retaining them.

### P2: Partial copositivity retains every reduction matrix

`check_stability()` allocates `bee_vee(r + 1)` at
`cpp/src/checkstab.cpp:97` and stores the complete reduction history. Each
iteration needs only the previous and current matrices plus scalar mask state.

Required outcome: use two rolling matrix buffers.

### P2: Exact LU uses a dense rational permutation matrix

`LU_Factorization` builds `m_P` as an `n*n` identity matrix at
`cpp/include/linalg/lu_factor_fraction.hpp:42`. Every solve then performs dense
rational multiply-add operations, mostly by zero, at
`cpp/include/linalg/lu_factor_fraction.hpp:118`.

Required outcome: store permutation indices and apply them by direct indexing.

### P2: Copositivity materializes all subsets before early rejection

`CopositivityCheckerV3::is_strictly_copositive()` builds every nonempty subset
at `cpp/include/linalg/copositive_fraction.hpp:159` before testing the first
principal submatrix.

Required outcome: generate fixed-cardinality masks on demand and stop at the
first failure. Benchmark any destructive/move LU variant separately.

## Build And Release

### P1: CI runs only for release tags

`.github/workflows/release.yml:3` has only a `v*` tag trigger. Pull requests and
ordinary pushes to `main` can therefore merge correctness failures without any
GitHub build or test signal.

Required outcome: run build plus fast tests on pushes and pull requests; keep
artifact publication restricted to release tags.

### P2: C++ tests cannot be disabled

`cpp/CMakeLists.txt:44` declares all four FetchContent projects unconditionally,
and `cpp/CMakeLists.txt:181` always adds tests. `BUILD_TESTING=OFF` is not wired,
so every build fetches and builds GoogleTest and the test targets.

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

CMake prefers `.so`/`.dylib` over static archives, and the local executable
currently resolves FLINT, MPFR, and GMP dynamically. The release workflow
installs those libraries on the runner but publishes only the executable, so an
end user still needs ABI-compatible libraries; the macOS binary also records a
Homebrew library path.

Required outcome: either publish a documented system-dependency package or
bundle/link the mathematical runtime consistently. Do not describe the current
artifacts as portable standalone executables.

## Current Validation State

- Current local Release build: successful.
- Core/CLI CTests: 10/10 passed.
- Full matrix CTests: 42/44 passed; only IDs 38-39 fail, as described above.
- Direct default and `--exact` CLI checks return zero ESS for IDs 36-37 after
  removal of the floating positive-definiteness certificate.
- Safe input `2#1/0,0,1` now exits normally with a zero-denominator parse error;
  safe-parser and CLI black-box regressions pass.
- The safe parser accepts 63 and rejects 0/64; the unchecked parser has no new
  branch and its boundary test now uses the maximum search dimension 63.
- All 40 expected-green fast verification matrices pass; long-run IDs 32/34
  and the separate known candidate-filter regressions 38/39 were excluded.
- Bit-position scanning now has an explicit nonzero precondition; lowest-bit
  mask extraction remains branch-free and defined for zero.
- GitHub access: working; local `HEAD` and `origin/main` both point to
  `78042fd` before the current uncommitted exact-PD fix.
- Latest `v0.33` Actions run: Linux and macOS built successfully but failed the
  old top-bit test; Windows and publication were then cancelled.
- The top-bit fix is committed; a new release remains blocked by verification
  IDs 38-39.
