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

## Allocation And Reallocation Audit

This audit inspected the explicit matrix, vector, string, and output-container
allocation sites in the active C++ source. Temporary variants were built under
`/tmp`; none of the variants below changed the production source tree.

The decision criterion is speed, not memory consumption. Allocation calls and
requested bytes were recorded only to locate repeated work. Saving kilobytes or
megabytes is not useful by itself. A change is retained only if native C++
timings show a repeatable end-to-end improvement without making the code harder
to read.

Release/LTO builds with native architecture enabled were pinned to CPU 2. Small
matrices were repeated to an approximately one-second target and reduced to the
median of the native `elapsed_ns` values returned by C++. A matrix already
slower than one second ran once. The balanced panels included dimensions 5, 10,
19, 20, 22, and 23, circular and non-circular games, and verified, unsafe, and
exact modes. Long one-sample results are reported as such; they are not treated
as precise percentage measurements.

### Result

No production allocation change is currently justified. Several variants
removed large numbers of allocator calls, but none produced a repeatable
end-to-end speed improvement. Two superficially attractive buffer-reuse changes
were materially slower. Keep the current implementation unless a later profile
shows a new runtime bottleneck.

### Exact solution scratch: neutral, do not change

Current code creates the `n`-entry local `matrix_frc solution` at
`cpp/src/find_candidate_exact.cpp:105` for every nonsingular support.

The readable prototype added `matrix_frc::resize()`, kept one `solution_` member
in `find_candidate_exact`, and resized only when its row count changed. A second
prototype called `resize()` unconditionally. That shorter-looking version was
bad: the internal no-op size test still ran on every exact solve and slowed the
balanced exact panel by 6.06% to 53.51%.

Native exact-mode results for the final row-count-check prototype were:

| Matrix | n | Type | Current ns | Reused solution ns | Delta | Unconditional resize ns |
|---:|---:|:---|---:|---:|---:|---:|
| 8 | 5 | circular | 2,583 | 2,500 | -3.213% | 3,833 |
| 9 | 5 | non-circular | 5,917 | 5,708 | -3.532% | 9,083 |
| 15 | 10 | circular | 67,875 | 67,833 | -0.062% | 84,750 |
| 75 | 10 | non-circular | 714,125 | 695,084 | -2.666% | 805,959 |
| 29 | 19 | circular | 115,896,082 | 114,286,208 | -1.389% | 128,316,769 |
| 51 | 20 | non-circular | 87,417 | 86,959 | -0.524% | 125,583 |
| 30 | 20 | circular | 360,488,913 | 360,282,163 | -0.057% | 383,316,372 |
| 32 | 22 | circular | 1,584,130,402 | 1,580,731,152 | -0.215% | 1,680,176,109 |
| 53 | 23 | non-circular | 1,632,023,610 | 1,609,975,193 | -1.351% | 1,736,225,192 |
| 33 | 23 | circular | 41,655,904,817 | 41,731,712,059 | +0.182% | not run |

The apparent small gains above did not repeat as a stable advantage. Three
matched runs on the two established comparison matrices were:

| Matrix | Current native medians (ns) | Reused-solution native medians (ns) | Median delta |
|---:|:---|:---|---:|
| 29 | 108,760,124; 108,793,728; 108,624,520 | 108,908,228; 108,850,562; 108,683,540 | +0.083% |
| 32 | 1,505,384,652; 1,514,806,902; 1,510,372,736 | 1,503,470,861; 1,505,999,819; 1,514,597,318 | -0.290% |

Verified and unsafe mode were likewise mixed rather than faster. The two
current columns are separate reference runs; the prototype column is the final
row-count-check variant.

| Matrix | n/type | Verified current 1 | Verified prototype | Verified current 2 | Unsafe current 1 | Unsafe prototype | Unsafe current 2 |
|---:|:---|---:|---:|---:|---:|---:|---:|
| 8 | 5/C | 9,791 | 9,791 | 9,292 | 3,583 | 3,542 | 3,542 |
| 9 | 5/N | 13,666 | 13,375 | 12,833 | 6,125 | 6,083 | 6,125 |
| 15 | 10/C | 80,250 | 79,708 | 77,833 | 26,208 | 26,250 | 26,250 |
| 75 | 10/N | 359,250 | 354,875 | 350,792 | 51,083 | 51,333 | 50,875 |
| 51 | 20/N | 212,542 | 213,375 | 209,083 | 119,917 | 119,625 | 119,750 |
| 30 | 20/C | 97,837,100 | 97,347,666 | 97,890,184 | 27,744,672 | 28,083,500 | 27,824,796 |
| 32 | 22/C | 547,325,552 | 543,382,412 | 540,481,072 | 211,055,078 | 214,249,415 | 213,172,829 |
| 53 | 23/N | 335,056,058 | 337,500,080 | 338,701,600 | 41,530,590 | 42,034,062 | 41,603,424 |
| 33 | 23/C | 1,215,178,418 | 1,216,855,614 | 1,224,713,586 | 375,861,023 | 374,682,788 | 374,800,648 |

Allocator instrumentation confirms that this was a real allocation hotspot,
but the absence of a speed gain controls the decision:

| Matrix | n/type | Current calls | Prototype calls | Current requested bytes | Prototype requested bytes |
|---:|:---|---:|---:|---:|---:|
| 51 | 20/N | 169 | 150 | 140,569 | 139,961 |
| 53 | 23/N | 123,752 | 16,489 | 16,035,391 | 864,191 |
| 33 | 23/C | 181,631 | 2,366 | 36,074,116 | 580,772 |

Candidate hashes stayed exactly equal for matrices 29 and 32:

- matrix 29: `9fe3562df089aa0db8b0e25c809791bbf98377838703d958e105d5709895f01b`
- matrix 32: `a372a12c795ec77206363bb892a6c302050a5b3817201326cbbcae7b8e38ea1f`

Decision: do not add the `resize()` API, result member, or solver branch for a
memory-only improvement.

### Reusing the next partial-copositivity matrix: slower, reject

Current code constructs `next_matrix` inside every reduction step at
`cpp/src/checkstab.cpp:178`. The prototype hoisted it outside the loop, resized
it, and swapped it with `previous_matrix` after each step. This is readable, but
it saved only 76 allocation calls and 186,240 requested bytes on matrix 32:

| Matrix | Current calls/bytes | Two-buffer calls/bytes |
|---:|:---|:---|
| 32 | 89,924 / 14,660,353 | 89,848 / 14,474,113 |

The three-run native timing comparison was materially worse:

| Matrix | Current median ns | Two-buffer median ns | Delta |
|---:|---:|---:|---:|
| 29 | 108,829,395 | 122,171,833 | +12.260% |
| 32 | 1,506,584,361 | 1,601,847,276 | +6.323% |

The exact candidate hashes above were unchanged. Decision: reject the change;
the current per-step construction is faster.

### Mutating the reusable Bee matrix directly: slower, reject

After positive definiteness fails, current code copies `Bee` into
`previous_matrix` at `cpp/src/checkstab.cpp:131`. A one-line prototype changed
that declaration to `auto& previous_matrix = Bee`, because the original Bee
values are not read again during that stability check. The reference is simple,
but reducing `bee_matrix_` destroys the full-size scratch capacity that later
candidates can reuse.

| Matrix | Current median ns | Bee-reference median ns | Delta |
|---:|---:|---:|---:|
| 29 | 108,829,395 | 115,688,249 | +6.302% |
| 32 | 1,506,584,361 | 1,613,763,776 | +7.114% |

The candidate hashes remained unchanged. Decision: preserve the current copy;
it is faster across successive candidates.

### Copositivity submatrix/LU/inverse reuse: kernel win, no application win

`is_copositive_hadeler()` constructs a principal matrix and exact LU for every
subset at `cpp/include/linalg/copositive_fraction.hpp:102-110`.
`LU_Factorization::inverse()` also constructs `b`, `bp`, `y`, and `x` for every
inverse column at `cpp/include/linalg/lu_factor_fraction.hpp:95-118`.

The prototype reused the checker submatrix and LU, allocated inverse `y` and
`x` once, solved unit columns directly through the stored permutation, and
scaled the adjugate in place. On the synthetic strictly-copositive,
non-positive-definite matrix with diagonal 1 and off-diagonal 2, it looked very
good:

| n | Current calls/bytes | Prototype calls/bytes | Current median ns | Prototype median ns | Time delta |
|---:|:---|:---|---:|---:|---:|
| 4 | 122 / 7,712 | 32 / 2,656 | 9,584 | 8,042 | -16.09% |
| 6 | 674 / 65,712 | 111 / 13,616 | 85,083 | 73,792 | -13.27% |
| 8 | 3,290 / 450,112 | 403 / 63,360 | 636,375 | 576,250 | -9.45% |
| 10 | 15,314 / 2,743,760 | 1,560 / 325,520 | 4,323,917 | 4,123,708 | -4.63% |

All 16 focused copositivity tests and all 7 exact-LU tests passed. The balanced
real exact-mode panel did not reproduce the kernel result:

| Matrix | n/type | Current ns | Copositivity prototype ns | Delta |
|---:|:---|---:|---:|---:|
| 8 | 5/C | 2,584 | 2,750 | +6.424% |
| 9 | 5/N | 5,916 | 6,375 | +7.759% |
| 15 | 10/C | 67,375 | 69,042 | +2.474% |
| 75 | 10/N | 703,584 | 701,917 | -0.237% |
| 29 | 19/C | 113,000,929 | 113,324,437 | +0.286% |
| 51 | 20/N | 87,375 | 88,750 | +1.574% |
| 30 | 20/C | 351,214,164 | 351,245,602 | +0.009% |
| 32 | 22/C | 1,556,901,218 | 1,571,947,938 | +0.966% |
| 53 | 23/N | 1,604,377,265 | 1,663,324,537 | +3.674% |
| 33 | 23/C | 41,633,209,989 | 61,915,798,177 | +48.717% outlier |

The last result did not repeat: a standalone prototype run took
42,498,045,495 ns, about 2.08% above that current reference. A later matched
three-run comparison was neutral: matrix 29 was +0.337% and matrix 32 was
-0.203%. The same-order real allocation panel also showed effectively no
change: matrix 51 was 169 versus 169 calls, matrix 53 was 123,752 versus 123,745,
and matrix 33 was 181,631 versus 181,631.

Decision: reject the three-file mathematical rewrite. The synthetic kernel is
not representative of the application workload, and the real program is not
faster.

### Candidate-vector reserve: no consistent speed gain

When candidate output is enabled, `cpp/src/fracessa.cpp:41-42` reserves
`250 * dimension_` candidate rows. The retained database has 31 circular games
with an average of 33.5 stored representatives and a maximum of 345; none
exceeds the reserve. The 56 non-circular games average 859.3 rows, reach 9,060,
and three exceed the reserve. Therefore the only plausible production change
was to skip the eager reserve for circular games while preserving it for
non-circular games.

A focused vector test showed the tradeoff:

| Shape | Eager reserve | No reserve | Median time delta |
|:---|:---|:---|---:|
| n=22, 156 rows | 157 calls, 934,912 bytes, 23,583 ns | 165 calls, 136,672 bytes, 22,000 ns | -6.712% |
| n=23, 9,060 rows | 9,062 calls, 6,094,080 bytes, 1,422,959 ns | 9,075 calls, 8,576,800 bytes, 1,526,126 ns | +7.250% |

Real verified-mode circular runs with candidate output enabled were:

| Matrix | n | Current reference ns | No-reserve ns | Delta |
|---:|---:|---:|---:|---:|
| 8 | 5 | 9,625 | 9,417 | -2.161% |
| 15 | 10 | 79,938 | 78,875 | -1.329% |
| 30 | 20 | 94,048,478 | 94,249,374 | +0.214% |
| 32 | 22 | 524,553,350 | 524,742,036 | +0.036% |
| 33 | 23 | 1,169,497,614 | 1,172,110,988 | +0.223% |

Unsafe-mode deltas for the same matrices were 0.000%, +0.627%, +0.260%,
-0.251%, and +0.214%, respectively. Decision: keep the current reserve. The
small verified cases improve slightly, but the result is not consistent across
modes or large matrices, and memory savings alone have no value here.

### LDLT diagonal storage: neutral, do not change

`matrix_frc::is_positive_definite()` currently allocates a matrix `L` and a
separate vector `d` at `cpp/include/linalg/matrix_fraction.hpp:99-102`. The
prototype stored `D` on the otherwise-unused diagonal of the same matrix and
treated the unit diagonal of mathematical `L` implicitly. This removes the
`d` allocation and the full identity initialization without adding a workspace
object.

The focused kernel changed from 2 allocation calls to 1. Requested bytes were
320 to 256 at n=4, 672 to 576 at n=6, 1,152 to 1,024 at n=8, and 1,760 to 1,600
at n=10. All 5 focused matrix tests passed, including exact positive
definiteness, and the exact matrix-29 and matrix-32 candidate hashes remained
unchanged.

The final timings were neutral:

| Workload | Current ns | LDLT-diagonal ns | Delta |
|:---|---:|---:|---:|
| matrix 29 exact | 108,829,395 | 108,931,248 | +0.094% |
| matrix 32 exact | 1,506,584,361 | 1,506,959,860 | +0.025% |
| matrix 61 unsafe, n=20 non-circular, 1,458 candidates | 364,755,330 | 367,593,371 | +0.778% |
| matrix 90 unsafe, n=23 non-circular, 9,060 candidates | 6,401,719,543 | 6,368,000,314 | -0.527% |

Decision: do not alter the clear existing LDLT representation for mixed noise
around zero.

### Allocation sites reviewed without a retained prototype

- Verified and unsafe candidate search already keep their double matrices and
  fixed-size index arrays as reusable members or stack storage. Their remaining
  MPFR work is one-time exact-to-binary conversion required by the proof.
- `linear_system_` and `bee_matrix_` already reuse their storage while the
  matrix size stays unchanged. The experiments above show that more aggressive
  Bee reuse can be slower.
- Support generators stream masks and retain only exact forbidden supports
  needed for pruning. Reserving those vectors speculatively has no measured
  speed basis.
- Parser storage is one-time input work. A prior compact-first reserve experiment
  gained about 2% on circular parsing but regressed symmetric parsing by 4-23%,
  so it remains rejected.
- Stored candidate vectors, Python strings/objects, and logging strings are the
  requested output. Avoiding them would change the API or logging behavior, not
  merely remove accidental allocation.
- The one input game-matrix copy is not in an enumerated hot loop. Changing
  constructor ownership to save it would broaden the API for no measured speed
  benefit.

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
- The dense exact LU permutation matrix is replaced by a source-row index
  vector. A fixed-seed audit found the same 1,984 nonsingular systems and exact
  old/new result hash across 2,000 matrices; every solve satisfied `A*x == b`
  and every inverse satisfied `A*inverse(A) == I`. A separate exhaustive audit
  covered all 46,233 row permutations through dimension 8, verified 2,856,811
  exact inverse entries, and reused one factorization while shrinking from
  dimension 8 to 2. Three alternating CPU-2 Release/LTO runs on the
  diagonal-1/off-diagonal-2 matrix, each using the median of 21 hot samples,
  improved determinant-only factorization by 3.80-5.86%, solve by
  10.76-17.97%, and inverse by 26.77-28.78% at dimensions 4, 8, and 12.
- The partial-copositivity reduction history is replaced by scalar mask/size
  state and previous/next matrix buffers. Each reduction is still logged before
  the buffers are exchanged. For the observed size-15-to-6, nine-step case,
  retained reduction storage falls from 1,185 to at most 421 exact fraction
  cells. Full exact candidate hashes for database IDs 29 and 32 were unchanged.
  Three CPU-2 Release runs per matrix, using native medians and an approximately
  one-second target, changed the median by -0.21% and +0.26%, respectively;
  runtime is therefore neutral within measurement noise.
- C++ tests now use CMake's standard `CTest`/`BUILD_TESTING` switch. With testing
  disabled, GoogleTest is neither declared nor fetched and the test subdirectory
  is not added. The always-supported CLI and Python module remain unchanged; no
  speculative CLI-only build option was added. A fresh disconnected
  `BUILD_TESTING=OFF` configure exposed no GoogleTest or test targets, built the
  CLI and Python module, and reported zero CTests. The default configuration
  still passed all 11 C++/CLI and 54 Python tests.
- Build configurations now use CMake's standard optimization and `NDEBUG`
  flags. Fresh GNU compile databases showed that every Debug FracESSA source
  had `-g` and none had `-O3`, `NDEBUG`, loop unrolling, `-march=native`, or
  LTO. RelWithDebInfo had only CMake's standard `-O2 -g -DNDEBUG`; Release had
  CMake's `-O3 -DNDEBUG` plus the project throughput flags and Release-only IPO.
  The verified proof source retained strict floating-point flags and no IPO.
  The fresh Debug CLI built and ran, while the Release suite passed all 11
  C++/CLI and 54 Python tests and ASan/UBSan passed all 11 C++/CLI tests.

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
  87 retained SQLite matrices. That long sweep was not rerun for the current
  local source and regression-test changes.
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
