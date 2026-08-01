# Project Knowledge

Last verified: 2026-08-01

## Source-Code Approval Gate

1. Do not modify any C++ or Python source file (`*.cpp`, `*.hpp`, or `*.py`)
   without Reinhard's explicit approval for the exact files and intended changes.
2. Before editing, work read-only, name every source file that would change, and
   describe the smallest necessary patch. A general request to investigate or
   fix an issue is not permission to modify additional or unlisted source files.
3. Approval is limited to the stated files and scope. If implementation requires
   another source file or a broader change, stop and request approval again.
4. Make only the approved minimal changes. Do not include opportunistic cleanup,
   refactoring, formatting, or adjacent fixes.
5. After editing, present the actual source changes and verification results for
   Reinhard's direct review and acceptance before continuing with further code
   changes.

## Line Width

1. Use 120 columns as a soft line-width target for C++, Python, comments, and Markdown.
2. Do not wrap lines merely to satisfy the traditional 80-column limit.
3. Exceed 120 columns when splitting a formula, command, URL, matrix, or readable expression would make it harder to understand.

## Priorities

1. Correctness is absolute; speed is second. Other concerns are secondary.
2. FracESSA evaluates a `2^n` support space and performs exact fractional
   operations millions or billions of times. It is optimized for many
   operations on small matrices, not a few operations on large matrices.
3. Keep validation at input boundaries. Do not add checks, allocations,
   abstractions, or branches to proven hot paths without a demonstrated
   correctness need and a benchmark.
4. Intentionally unchecked bitset operations exist for raw speed. Matrix input
   has one validating parser at the input boundary; values up to 18 digits use
   direct integer construction and larger values use exact FLINT text conversion.
5. Use the Ponytail skill for code work: understand the complete path, then use
   the smallest correct implementation.
6. Treat allocation counts and ordinary kilobyte/megabyte memory reductions as
   diagnostics, not as performance goals. Retain an allocation optimization
   only when real end-to-end benchmarks show a repeatable speed improvement;
   lower memory use does not justify slower or more complicated code.
7. Performance samples must cover the affected methods, small and large games,
   and both circular and non-circular matrices. Include dimensions around 20
   and 23 when feasible. A synthetic kernel result alone is not sufficient.
8. Keep dimension-2 matrices in the canonical test data and correctness
   verification, but exclude them from performance benchmark runs, tables, and
   aggregate performance statistics. Dimension 3 remains benchmarked for now.
9. Keep numerical code human-readable. Prefer a small reuse or two-buffer
   change when it is measurably faster, but do not add custom allocators, pools,
   generic workspace frameworks, or extensive plumbing merely to reduce
   allocations.
10. FracESSA algorithm and orchestration code under `cpp/include/fracessa/` and `cpp/src/` uses C++ numerical types only. Raw
    FLINT types and `fmpz_*`/`fmpq_*` calls are confined to the thin `linalg` wrappers and specialized numerical kernels.

## Repository

- Root: `/home/reinhard/projects/fracessa`
- Remote: `git@github.com:reinhardullrich/FracESSA.git`
- Main implementation: `cpp/`
- Python API: `python/`
- Canonical test data: `testdata/fracessa_testdata.sqlite3`
- Research papers: `research/papers/`
- Historical benchmark material: `experiments/`
- Inactive historical tooling: `archive/`
- Agent documentation: `aidocs/`
- Public GitHub introduction: `README.md`
- `AGENTS.md` must remain a pointer only.

Generated or local-only paths include `cpp/build*/` and experiment `builds/`,
`sources/`, and `logs/` directories.
`archive/callgrind/` preserves the four former JSON-fed profiling scripts
unchanged for reference. They are not active tooling and do not run against the
current SQLite matrix store without adaptation. The tracked
`cpp/callgrind/callgrind.out.1` through `.35` files are historical Callgrind
3.15.0 profiles from the former x86-64 Linux build; newly generated
`callgrind.out.*` files remain ignored unless added deliberately.
`zzz_legacy/` is the tracked collection of preserved REF/EFR predecessors.
Its six top-level folders are `EFR`, `REF_2016-10-06`, `REF_2016-11-16`,
`REF_2016-11-20-Werner`, `REF_2019-09-20`, and `REF_R`. They are preserved historical
material, not active project source. `EFR/` contains the four selected C# timeline snapshots
`EFR_2016-04`, `EFR_2016-09`, `EFR_2018-03`, and `EFR_2019-08`, plus the sole
`NewRational.2.1` version as `NewRational`. These project directories were moved
directly from their former locations without copying. The other two
byte-identical `NewRational.2.1` copies are in the desktop Trash.
The duplicate April and September 2016 timeline trees and both copies of the
skipped intermediate September 2016 port are in the desktop Trash.
`REF_2016-10-06` is the source-only GMP/hybrid milestone, `REF_2016-11-16` is
the last source-only version before circular-symmetry optimization,
`REF_2016-11-20-Werner` is the Werner circular-symmetry version, and `REF_2019-09-20`
is the substantial 2019 rewrite with its minor December 2025 timing and build
changes retained. `REF_R` contains one copy of each preserved research script
and the newest EFR test driver, dated October 21, 2016. The former mixed-history
tree and its duplicate, intermediate, generated, profiling, IDE, and unrelated
material are in the desktop Trash.
The backed-up private Werner download page was last updated on November 19,
2016, matching the email that directed Werner to the new circular-symmetry
release and his confirmation that it worked. The downloaded source archive's
algorithm files match `REF_2016-11-20-Werner`; the 2019 Werner correspondence
contains no later executable or download.
The collection is intentionally storage-cleaned rather than independently
buildable in every folder: `REF_2019-09-20/dependencies/` retains the verified
Boost 1.71 tarball through Git LFS, and the canonical Boost 1.62 tree remains under
`REF_2016-11-20-Werner/include/boost`. Redundant extracted Boost trees and
generated `build*`/`obj` directories are not retained.
The legacy collection intentionally contains no nested Git-control metadata;
its former `.git` directories were removed from the collection.

## Product Surface

FracESSA is a C++17 ESS (evolutionarily stable strategy) analyzer for symmetric
payoff matrices. It has two entry surfaces:

- CLI: `cpp/src/main.cpp`, built as `cpp/build/fracessa`.
- PyFracESSA package: `python/fracessa/`, backed by the native extension in
  `cpp/src/pybind_module.cpp` named `fracessa_core`.

CLI matrix format is `dimension#values`. Values are either the upper triangle
of a symmetric matrix (`n*(n+1)/2` entries) or the compact circular-symmetric
form (`floor(n/2)` entries).

The parser accepts dimensions 1 through 63, rejects textual fractions with a
zero denominator, and validates value syntax during the same scan that builds
the exact fractions. Support masks have 64 storage bits, but complete
enumeration requires the exclusive `2^n` bound, so dimension 64 is not
supported.

Every entry surface requires one of three search methods before the matrix; there is no default. `safe` uses exact arithmetic for
every candidate decision. `fast` uses the unnormalized raw-double heuristic only after an exact integer precision-span check and
then confirms every surviving support exactly. A matrix with $P\geq10^9$ switches entirely to `safe`, and a support with a pivot
below $10^{-12}$ reaches exact checking. Its remaining probability and outside-payoff rejections are heuristic. Experimental
`test` is an independent source copy of `fast`, not a wrapper around it; it currently has identical behavior and can be changed
without changing production fast search.

## Computation Flow

```text
CLI or pybind input
  -> matrix_parser
  -> fracessa constructor
  -> support generation and pruning
  -> optional find_candidate_fast or find_candidate_test heuristic
  -> find_candidate_safe
  -> check_stability
  -> ESS/candidate output
```

Important implementation points:

- Supports are represented by `uint64_t` masks.
- Fixed stack buffers hold extracted support indices; clear-lowest-set-bit
  iteration avoids bit-test branches in inner matrix loops.
- `--fullsupport` constructs and checks the full mask directly; on failure its
  fallback solves only cardinalities below `n`.
- Production generates one mask at a time through a header-defined templated
  callback and never materializes the complete frontier or a cardinality layer.
  The generator owns the cardinality sweep and calls the analyzer with both the
  mask and its size. `NonCircularSupportGenerator` uses fixed-cardinality binary DFS;
  `CircularSupportGenerator` uses fixed-content FKM necklace recursion plus
  reflection reduction. Both prune partial branches against earlier exact
  candidate supports bucketed by lowest bit. Circular rules expand every
  distinct rotation/reflection only as compact forbidden masks. The analyzer
  stores one solved bracelet representative with their count as `multiplier`.
  See `plans/SUPPORT_GENERATORS.md`.
- Newly found exact candidates are pending until the next cardinality, keeping
  each generator layer's pruning rules stable. Stability is irrelevant to this
  pruning rule: every exact equilibrium support forbids later strict supersets.
- `find_candidate_fast` owns the direct, unnormalized binary64 conversion and bordered solve. One exact matrix-wide precision-span
  test reads the safe solver's already prepared integer game and common denominator. A matrix with $P\geq10^9$ switches to safe
  search before double allocation or conversion. Every pivot below $10^{-12}$ is inconclusive and reaches exact checking. Its
  probability and outside-payoff branches remain heuristic, so `fast` is not a correctness certificate.
- `find_candidate_test` intentionally duplicates fast double storage and the candidate kernel without sharing its implementation
  or double state. It currently has the same precision-span, pivot, probability, and outside-payoff decisions as fast.
- `find_candidate_safe` clears the game's rational denominators once, eliminates
  the normalization/payoff border, and constructs the integer-scaled symmetric
  reduced system $dH y=dr$ in reusable FLINT storage. One fraction-free
  $LDL^T$-style factorization solves the candidate, proves singularity, and
  records the exact inertia needed by stability. Rational values are constructed
  only for successful public candidate output.
- `fracessa` owns the rational game used by stability. `find_candidate_fast` and `find_candidate_test` refer to it for their double
  matrices, while `find_candidate_safe` owns the one integer-scaled copy used by all exact candidate solves and both one-time
  precision-span decisions.
- Exact candidate factorization and validation use FLINT `fmpz_t` integers;
  header-only `linalg::integer` and `linalg::matrix_int` wrappers own their storage and expose inline exact operations to
  FracESSA. Public rational results and the retained Bomze stability fallback use FLINT `fmpq_t` through `linalg::fraction`.
- Stability reuses the exact reduced-Hessian inertia. A non-negative-definite
  support Hessian rejects ESS immediately; a negative-definite Hessian proves
  ESS immediately when extended support equals support. Only the rare
  negative-definite case with outside best replies constructs Bee and enters
  the retained Bomze reduction/copositivity path. A binary64 result is never
  accepted as a final mathematical certificate.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md` documents the concrete failures and
  proves why tolerance tuning cannot recover an exact PD certificate.
- `correctness/FAST_CANDIDATE_FALSE_REJECTION.md` gives exact ESS counterexamples for all three former fast per-support rejection
  rules. Current fast recovers the cutoff example through per-support pivot fallback and recovers the probability and
  outside-payoff examples because their precision spans select matrix-wide safe search. These fallbacks are heuristics, not a
  general correctness proof.
- `safe` does not initialize or allocate any double candidate-search state, and all final stability decisions use exact rational
  arithmetic.
- The raw-double algorithm at revision `32f61679da64` used six one-time input checks and rejected small pivots. Current fast instead
  uses the exact precision-span gate and treats small pivots as inconclusive. The retired normalized heuristic fixed IDs 38-39 but
  introduced misses on IDs 45-47 and is not a production method.

Key files:

- `cpp/include/fracessa/bitset64.hpp`: support-mask primitives.
- `cpp/include/fracessa/supports.hpp`: support generation and pruning.
- `cpp/include/linalg/integer.hpp` and `cpp/include/linalg/matrix_integer.hpp`: header-only owning C++ wrappers around FLINT exact
  integers and integer matrices.
- `cpp/include/linalg/fraction.hpp`: FLINT rational wrapper.
- `cpp/include/linalg/copositive_fraction.hpp`: exact copositivity checks.
- `cpp/include/fracessa/find_candidate_fast.hpp` and
  `cpp/src/find_candidate_fast.cpp`: production raw-double class, exact precision-span gate, small-pivot fallback, heuristic
  inequalities, and reusable scratch.
- `cpp/include/fracessa/find_candidate_test.hpp` and
  `cpp/src/find_candidate_test.cpp`: independent experimental copy of fast search.
- `cpp/include/fracessa/find_candidate_safe.hpp` and
  `cpp/src/find_candidate_safe.cpp`: exact class, border elimination,
  integer candidate validation, and candidate construction.
- `cpp/include/linalg/flint_style_fraction_free_ldlt.hpp`: reusable in-place
  fraction-free symmetric solve, exact inertia, and zero-diagonal coordinate
  handling.
- `cpp/include/fracessa/fracessa.hpp` and `cpp/src/fracessa.cpp`: exact game
  ownership, method coordination, support search, and candidate lifecycle.
- `cpp/src/checkstab.cpp`: stability classification.

## Build And Dependencies

From repository root:

```bash
./build.sh
```

Equivalent manual build:

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build -j"$(nproc)"
```

Required system dependencies are a C++17 compiler, CMake 3.18 or newer, Python
3.14 or newer with development headers, GMP, MPFR, and FLINT. CMake
FetchContent downloads:

- `spdlog`: optional rotating diagnostic logs.
- `argparse`: the cross-platform CLI parser.
- `googletest`: C++ unit-test executables only; it is fetched only when
  `BUILD_TESTING=ON` and is not linked into the production executable.
- `pybind11` v3.0.4: the native Python module.

`BUILD_TESTING` uses CMake's standard `CTest` option and defaults to `ON`.
Configure with `-DBUILD_TESTING=OFF` to skip GoogleTest and every C++/CLI test
target. The other three FetchContent dependencies remain part of production, so
a clean configure still needs network access unless their sources are cached.

Local non-MSVC Release builds default to `FRACESSA_NATIVE_ARCH=ON`
(`-march=native`); Release CI sets it to `OFF`. Debug and other configurations
use CMake's standard flags without FracESSA's throughput options. IPO/LTO is
enabled only for Release and only when CMake confirms support.
The fast heuristic uses the normal Release throughput settings of its historical implementation. MPFR remains an explicit link
dependency because static FLINT linkage requires the `FLINT -> MPFR -> GMP` order; production source no longer calls MPFR
directly.
When sandboxing blocks the normal ccache directory, rerun the build with
escalated filesystem access rather than disabling or redirecting ccache.

`cmake --install cpp/build` installs the CLI target only. It does not install
GMP, MPFR, or FLINT.

## Tests And Test Data

```bash
./test.sh # build, core/CLI tests, and wrapper tests
```

The non-matrix CTest suite consists of nine GoogleTest executables plus one CLI
black-box parser test. Wrapper tests use Python `unittest`. Matrix correctness
is no longer wired as one CTest per matrix.

`testdata/fracessa_testdata.sqlite3` is the canonical matrix, expected-result,
and timing store. Its strict schema is in `testdata/schema.sql`; the current
snapshot has 90 matrices and 49,161 stored candidate representatives. Nullable
multipliers recover weighted totals of 86,156 candidates and 83,381 ESS:
circular rows store one smallest dihedral representative and its orbit count,
while non-circular rows store null. Candidate IDs and row order remain
reproducibility checks; complete weighted candidate sets and ESS
classifications are the mathematical contracts.

Every dimension from 2 through 25 has at least one circular and one
non-circular matrix. IDs 67-79 fill the previously missing combinations with
deterministic random integers; exact and the then-current normalized heuristic agreed on their complete rational candidate
contracts before insertion.

IDs 45-47 preserve the verified-search regressions against the retired normalized heuristic: the dimension-20 heuristic
counterexample, the LU-boundary fallback case, and the failed-proof exact-fallback case. IDs 48-55 cover non-circular dimensions 15-24 through
Hilbert, Hadamard, Paley conference, MINIJ, Fiedler, deterministic random
families, and a dense weighted-Laplacian game with one full-support ESS.

IDs 91-93 are exact non-circular false-rejection regressions for the former fast rules covering all three per-support candidate
conditions. ID 91 has a
nonsingular full-support ESS but produces a $7.5\times10^{-13}$ double pivot below the $10^{-12}$ cutoff. ID 92 has an exact
positive probability $10^{-10}$ that the double solve computes as negative. ID 93 has an outside payoff exactly $10^{-4}$ below
the equilibrium payoff that the double solve computes above its rejection margin. Current fast and test recover them through the
small-pivot or matrix-wide precision-span fallback; their fast/test/safe ESS counts are `1/1/1`, `1/1/1`, and `2/2/2`.

Every distinct matrix in Tables 1 and 2 of the Bomze-Schachinger-Ullrich
ESS-growth paper is present exactly once. IDs 18 and 26 are the exact published
Table 1 matrices for dimensions 12 and 17. IDs 80-81 are the two previously
missing Table 2 circular base matrices, and IDs 82-90 are its nine constructed
non-circular matrices. Same-property alternatives formerly stored at IDs 12
and 21 were removed; the former contents of IDs 18 and 26 were replaced by the
published vectors.

The timing snapshot has 670 CPU-2 persistent-Pybind median rows with a one-second target. Werner's default and the
preserved pre-mode default are stored as `fast`; Werner's exact run is stored as `safe`, matching their current semantic
equivalents. The later `current-main` three-mode snapshot retains historical `safe`, `unsafe`, and `exact` labels because its
`safe` rows are the removed verified proof rather than today's exact safe method. Build label and revision disambiguate them. The
raw historical build mismatches IDs 38-39, the retired normalized heuristic mismatches IDs 45-47, and the removed verified proof
and exact search match all 87 matrices in that timing snapshot. The historical `fast` build at revision `8697ebaf` has 78
rows covering every matrix with dimension at least 3 that existed when it was measured; all observed ESS counts match. New IDs
91-93 have no stored timings yet. Timing
reports include matrix dimension, circularity, and the derived paper-style
lower bound `gamma_lower_bound = expected_ess ** (1 / dimension)` without
storing it in SQLite.

The first `reduced-hessian-ldlt` benchmark measured 85 matrices; IDs 33-34 were
not included in that run. All 85 ESS counts match. Against `current-main` on
the same matrices, summed exact
medians fall from 1,386.743 to 1,184.045 seconds (14.62%), and the median
per-matrix ratio is 0.6842. Eighty-two matrices are faster. IDs 45 and 47 are
material regressions at 168.882 versus 74.655 seconds and 132.230 versus 72.686
seconds; excluding those two adversarial cases, summed time improves by 28.76%.
ID 46 is 3.8 milliseconds slower (13.4%); its absolute effect is small, but a
repeat would be needed before classifying it as signal or timing noise.

The isolated integer-solver experiment under
`experiments/exact_integer_solver_comparison_2026-07-31/` compares the old
bordered rational Gaussian solver, current reduced-Hessian rational $LDL^T$,
integer bordered FFLU, a complete FFLU-plus-candidate-$LDL^T$ hybrid, and
an FFLU-plus-fraction-free-reduced-Hessian hybrid, as well as rational
fraction-free FFLU. Its CPU-2 one-second sweep covers 82 matrices;
IDs 45, 47, 65, 66, and 90 are excluded because current reduced-Hessian exact
time exceeds two minutes. IDs 33-34 are ordinary included rows. All candidate
contracts match. The fraction-free reduced-Hessian kernel matches current
exact inertia in detailed candidate comparisons and in the 74-row
cross-procedure portion of the sweep. Its summed medians are 83.198 seconds
versus 326.256 seconds for current $LDL^T$ (74.50% lower) and only 0.30% above
candidate-only FFLU. It loses all 26 dimension-2-to-8 rows, wins 26 of 28
dimension-9-to-16 rows, and wins 27 of 28 dimension-17-to-25 rows; the large
exception is ID 51, whose 20 visited supports are all candidates. This supports
the fraction-free reduced-Hessian solver used by production. The original
experiment remains the immutable comparison against the former rational kernel.

The first production comparison against preserved rational revision `799be715`
used one persistent Pybind process per build, native nanosecond medians, a
one-second target, CPU 2, and nine circular/non-circular matrices spanning
dimensions 3 through 23. The fraction-free exact path improved the median
per-matrix time by 78.33% and the arithmetic mean percentage by 56.32%. The
only regression was the dimension-4 circular ID 69, from 0.792 to 0.833
microseconds. Substantive cases improved by 61.48% to 94.84%; ID 51 improved
by 2.42% because its 20 visited supports are all candidates.

The former JSON/CSV verification, baseline-generation, speed-benchmark, and
Callgrind runners were removed. There is no replacement matrix-verification
runner yet. The small `python -m fracessa.timing` tool reads matrices from
SQLite, measures one build and one matrix at a time on a user-selected Linux
CPU, and writes normalized nanosecond samples to the same database. Reusing a
session name groups separately invoked builds. Each row records `source_ref`
(for a moving name such as `main`), its immutable `revision`, the binary hash,
backend, historical `mode` database value, CPU, comment, observed ESS count, target and measured wall time,
iteration count, and median native elapsed time. One Pybind process stays open
for every selected method and matrix in a build. New runs require `--method fast` or `--method safe`; adapters map those choices
onto old Pybind and CLI interfaces when benchmarking historical builds. The first returned native
duration chooses `ceil(target / duration)` total samples and remains in the
sample; a duration at or above the target chooses one run. The stored result is
the median returned `elapsed_ns`. Measured wall time is metadata only and
never chooses or determines the reported timing. The CLI backend remains
available for legacy inspection but starts a child process per sample and must
not be mixed with persistent-Pybind microbenchmarks. Dated material under
`experiments/` and `aidocs/experiments/` remains immutable historical
evidence.

Database IDs 45-47 preserve the retired normalized-heuristic correctness regressions tracked in `reviews/CPP_REVIEW.md`, and IDs
91-93 preserve the three former fast per-support false rejections. The wrapper integration suite checks IDs 38 and 46 through
fast and safe routes, but no complete SQLite matrix suite is currently wired into `./test.sh` or release CI.

## Pybind Boundary

`cpp/src/pybind_module.cpp` exposes the C++ analyzer as the native
`fracessa_core` module. It owns Python/native argument and result conversion,
native status codes, GIL release, and native timing. Binding-specific open
findings are tracked in `reviews/PYBIND_REVIEW.md`, separately from both the
analyzer core and Python orchestration.

The safe parser throws `std::invalid_argument` with a detailed diagnostic and
does not write to `stderr`. The CLI catches and prints that diagnostic; Pybind
catches the same exception and exposes one `PARSE_ERROR` status with the
diagnostic in `error_message`, without reparsing the input.

The analyzer and native binding both store the ESS count as `size_t`; Pybind
converts that value directly to Python's arbitrary-precision integer.

Native analyzer timing uses `std::chrono::steady_clock` and is always returned
as integer nanoseconds in `elapsed_ns`. The CLI `--timing` output uses the same
clock and unit. There is no wrapper timing-suppression option.

Matrix IDs are signed 64-bit values at the CLI, analyzer, Pybind, and file-sink
boundaries. `Matrix` accepts only built-in Python integers in that range and
rejects booleans and coercible float/string values before native execution.

The binding releases the GIL during native execution. Logging-enabled calls are
serialized by one process-wide native mutex because each analyzer writes and
rotates the same log file; non-logging calls remain concurrent.

## PyFracESSA

`python/fracessa/` calls `fracessa_core` in-process and is the maintained API.
It supports sequential execution, process-based parallelism across matrices,
and CSV/JSON/Parquet disk sinks. One matrix is always computed by one worker
process; parallelism is across matrices. File sinks create output paths
exclusively and never overwrite existing run-ID output. Each format writes
summary and candidate data plus a format-specific JSON sidecar for per-matrix
metadata. Empty outputs have stable readable schemas; Parquet buffers 1,024
rows per row group. JSON writers replace non-finite floats with `null` and use
strict encoding, so they never emit the non-standard `NaN` or `Infinity`
literals.

Sink construction and consumption are transactional across each exclusive
output triplet. A caught initialization, computation, write, or finalization
exception closes the result iterator and sink resources, removes only paths
made by that attempt, and re-raises the original exception so callers may retry
the same run ID. Cleanup errors never replace an active operation error.

All maintained wrapper execution paths return one flat dictionary. Candidate
rows are plain dictionaries; there are no Python result-row classes. Circular
rows contain one bracelet representative with an integer `multiplier`;
non-circular rows use `None`. `candidate_count` is the number of returned
representatives, while `ess_count` remains the weighted mathematical total.

No production wrapper or matrix workflow imposes a per-matrix computation
timeout. A matrix may legitimately run for hours. Worker-liveness handling must
not be implemented as a computation timeout.

`RunConfig()` contains only analysis options. `run`, `run_multiprocessing`, and `compute_matrix` require `"fast"`, `"safe"`, or
experimental `"test"` as their first argument; there is no default method. Fast and test can miss candidates and ESS results,
while safe bypasses every floating-point candidate procedure.

`run` and `run_multiprocessing` are the only public execution functions. Both
accept a required method followed by one `Matrix` or an iterable and accept an optional sink. One matrix
returns one dictionary, an iterable returns a result iterator, and passing a
sink eagerly writes the results and returns the matrix count. `compute_matrix`
is the public low-level native adapter. `run_multiprocessing` adds only a final
optional `MPConfig`; its default uses the CPUs available to the Python process.

The private multiprocessing helper uses one shared matrix queue and one shared
result queue, yields completion order, bounds pending matrices to
`min(queue_maxsize, workers * prefetch_per_worker)`, serializes each matrix
before counting it, detects dead workers while waiting, and cancels workers
when iteration stops early. It does not batch multiple matrices into one queue
item. Native logging is rejected before multiprocessing workers are created;
it remains available in sequential wrapper execution.

New API work belongs in `python/fracessa/` and `fracessa_core`.

Every production Python module, class, function, and method has a standard
docstring consumable by `pydoc`, `pdoc`, or Sphinx `autodoc`. Tests are excluded
because they are verification code rather than generated API documentation.

The generic JSON loader accepts a top-level row list or an object containing
the configured matrix key. It requires a list of object rows and rejects a
missing key or malformed row container instead of silently returning no matrices.
It validates integer/string fields without lossy coercion; values-only matrices
require a built-in integer dimension.

## Release Workflow

`.github/workflows/release.yml` builds and checks Ubuntu and macOS with Python
3.14 for pull requests and `main`; feature branches are not built a second time
by the push trigger. Windows is temporarily restricted to pushed `v*` tags
until its dependency installation is fast enough for normal CI. Native
integration tests require the built module, and each platform build installs
PyArrow before the wrapper suite, so binding and Parquet coverage cannot turn
into successful skips. Packaging, artifact upload, write permission, and GitHub
release publication run only for `v*` tags.

The artifacts are architecture-specific and are not uniformly self-contained:

- Linux is x86-64 and currently depends on system FLINT at runtime.
- macOS is ARM64 and currently records a Homebrew FLINT dylib path.
- Windows is x86-64 and currently links its third-party/runtime dependencies
  statically.

GitHub installing dependencies on its runners makes compilation succeed; it
does not install those dependencies on an end user's machine. Published
`v0.22` and `v0.24` macOS binaries used statically linked mathematical
libraries, but the current release configuration uses system FLINT.

## Documentation Policy

- `KNOWLEDGE.md` contains only current facts and durable project policy.
- `CHANGES.md` is the append-only human-readable history. Update it when a
  meaningful project change benefits from a concise historical record; read it
  only when history or drift matters.
- `reviews/CPP_REVIEW.md`, `reviews/PYBIND_REVIEW.md`, and
  `reviews/PYTHON_REVIEW.md` contain only unresolved findings for their
  respective implementation areas.
- Dated experiment reports are immutable snapshots and must state their scope.
- Git remains authoritative for exact diffs and commit history.
- Do not store generated source concatenations, session transcripts, or stale
  profiling output in `aidocs/`.
