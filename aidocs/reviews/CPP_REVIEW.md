# C++ Review

Last verified: 2026-07-26

Scope: active C++ production code, pybind boundary, CMake, C++/CTest
coverage, and the release workflow. Frozen experiment source copies are excluded
except where a dated result is cited as evidence.

Correctness is ranked before speed. This file contains unresolved findings only;
remove a finding after its fix and regression coverage are complete.

## Correctness

### P0: Double positive-definiteness creates false ESS results

`fracessa::check_stability()` accepts a successful double Cholesky test as a
final mathematical certificate at `cpp/src/checkstab.cpp:50`. Rounding can make
an exactly singular or indefinite Bee matrix appear positive definite. The
`--exact` flag still executes this floating stability decision, despite the CLI
describing the mode as exact-only.

The same path can accept non-finite arithmetic: sufficiently large exact input
can convert to infinity, Bee construction can produce NaN through operations
such as `inf - inf`, and the current `diagonal_value <= tolerance` comparison
does not reject NaN.

Verification IDs 36 and 37 reproduce the defect. Both default and `--exact`
currently return two ESS, while exact positive-definiteness and copositivity
return zero.

The smallest correct fix is to delete the double Bee/PD decision and enter the
exact rational PD check directly. The existing exact-PD experiment measured
about 10% overhead for repeated small matrices and about 1% across the 35
historical benchmark matrices; see
`../experiments/speed_comparison_2026-07-26/MICROBENCHMARK_COMPARISON.md`.

Required outcome: no floating-point result may be a final positive certificate,
and `--exact` must not execute a floating stability decision.

### P0: The double candidate filter drops valid ESS supports

`solve_linear_dbl()` rejects pivots using fixed absolute `1e-12` thresholds at
`cpp/include/linalg/linear_solver.hpp:47` and
`cpp/include/linalg/linear_solver.hpp:63`. `find_candidate_dbl()` treats every
such numerical rejection as conclusive at `cpp/src/findeq.cpp:42`, so exact
arithmetic never sees the support.

The outside-support test is also conclusive and uses the fixed absolute margin
`1e-4 * dimension` at `cpp/src/findeq.cpp:63`, with rejection at
`cpp/src/findeq.cpp:72`. It therefore has the same scale-dependence problem even
when the linear solve succeeds.

Verification IDs 38 and 39 both have the mixed ESS `(1/2,1/2)`. The default
path returns zero ESS, while `--exact` returns one. Scaling or translating the
same game therefore changes the program's answer.

Required outcome: remove the double candidate filter from the decision path and
benchmark the exact-only path. Falling through to exact arithmetic after every
double rejection saves no exact work, so keeping the current filter would only
add work. Reintroduce a numerical filter only if a rigorously one-sided method
can reject supports without false negatives.

### P1: A zero denominator crashes the safe parser

The safe path passes a token directly from `parse_fraction_token()` at
`cpp/src/matrix_parser.cpp:13` to the string constructor at
`cpp/include/linalg/fraction.hpp:52`. FLINT/GMP accepts the textual shape
`1/0`, leaving an invalid rational that later crashes during numeric use.

Reproduced against the current Release build:

```text
./cpp/build/fracessa '2#1/0,0,1'
Segmentation fault (exit 139)
```

Required outcome: reject a zero denominator in the shared safe string
constructor before canonicalization or conversion, return a normal parse error,
and add both fraction-level and CLI black-box regressions. This validation is at
the input boundary and does not add a hot-path branch.

### P2: FLINT-allocated strings are freed with the wrong API

`fraction::to_string()` and `operator<<` call `free()` at
`cpp/include/linalg/fraction.hpp:248` and
`cpp/include/linalg/fraction.hpp:259`. FLINT's `fmpq_get_str(nullptr, ...)`
allocates through `flint_malloc`, so its matching deallocator is `flint_free`.
Plain `free` happens to work with the current default allocator but is not the
FLINT contract and can fail with another FLINT memory manager or CRT boundary.

Required outcome: use `flint_free` at both sites and remove the incorrect
`<cstdlib>`/`free()` comments.

### P2: The supported dimension contract is inconsistent at 64

The bitset header advertises `n <= 64`, and the unsafe parser test explicitly
accepts dimension 64. The search uses `1ULL << dimension` in
`two_to_the_power_of()` at `cpp/include/fracessa/bitset64.hpp:62`; shifting by
64 is undefined, and enumerating all 64-bit masks cannot use a one-past-end
`uint64_t` sentinel.

Required outcome: distinguish 64 storage bits from the supported search range
of dimensions 1 through 63. The user-required unsafe path need not add a runtime
check. Keep the parser-only 64-dimensional storage test, but document that an
analyzer call requires `dimension < 64`.

### P2: The safe parser accepts partial dimensions and pybind guesses errors

The safe parser does not check complete `std::stoull` consumption at
`cpp/src/matrix_parser.cpp:40`, so input such as `2x#0,1,0` currently succeeds.
Pybind then reparses failed input in
`infer_safe_parse_status()` at `cpp/src/pybind_module.cpp:59`, which can report
an invalid rational as `INVALID_VALUE_COUNT`.

Required outcome: require `std::stoull` to consume the complete dimension. The
distinct status codes are part of the public wrapper API, so return one parser
status to both CLI and pybind callers instead of maintaining a second parser to
classify the same failed input.

### P3: Pybind narrows the ESS count

The analyzer stores `ess_count_` as `size_t`, but `NativeResult::ess_count` is
`int` at `cpp/src/pybind_module.cpp:41` and receives an explicit narrowing cast
at `cpp/src/pybind_module.cpp:114`. The search space can exceed `INT_MAX`, while
Python integers do not require this truncation.

Required outcome: preserve `size_t` or `uint64_t` through the native result.

### P3: The CLI candidate header has an extra empty column

`candidate::header()` ends with a semicolon at
`cpp/include/fracessa/candidate.hpp:56`, while candidate rows do not. A standard
CSV reader therefore sees one unnamed header column with a missing row value.

Required outcome: remove the trailing delimiter and add a black-box assertion
that header and row field counts match.

### P3: Zero-input bit scanning has an invalid test contract

On GCC/Clang, `ctz64(0)` reaches `__builtin_ctzll(0)` at
`cpp/include/fracessa/bitset64.hpp:44`, which is undefined. Tests nevertheless
assert a sentinel result for `find_pos_first_set_bit(0)` and a zero result for
`lowest_set_bit_as_bit(0)`. Current production callers provide nonzero masks,
so this is a helper-contract and test defect rather than a reproduced
production failure.

Required outcome: keep count-trailing-zeros explicitly unchecked for nonzero
input and remove the test that demands a sentinel from it. Implement
`lowest_set_bit_as_bit()` with the direct unsigned lowest-bit expression, which
naturally returns zero for zero without a branch.

## Speed

### P1: `--fullsupport` builds all supports before checking one

`fracessa` calls `supports_.initialize()` at `cpp/src/fracessa.cpp:50` before
the full-support branch. At dimension 24 this generates 16,777,215 masks before
checking the single requested full support.

Required outcome: construct and check the full-support mask first; initialize
fallback supports only if that check fails.

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
- Full matrix CTests: 40/44 passed; only IDs 36-39 failed, as described above.
- Wrapper unittests that exercise the native module: 23/23 passed.
- Direct CLI probes reconfirmed the safe-parser `1/0` segmentation fault and
  the default/exact result split for IDs 36-39.
- UBSan reconfirmed both `ctz64(0)` and `1ULL << 64` as undefined behavior.
- GitHub access: working; local `HEAD` and `origin/main` both point to
  `81a6849` (`v0.33`) before the current uncommitted review/fix work.
- Latest `v0.33` Actions run: Linux and macOS built successfully but failed the
  old top-bit test; Windows and publication were then cancelled.
- The current worktree fixes that top-bit shift, but a new release remains
  blocked by verification IDs 36-39.
