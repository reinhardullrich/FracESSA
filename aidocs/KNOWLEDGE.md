# Project Knowledge

Last verified: 2026-07-31

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
7. Performance samples must cover the affected modes, small and large games,
   and both circular and non-circular matrices. Include dimensions around 20
   and 23 when feasible. A synthetic kernel result alone is not sufficient.
8. Keep numerical code human-readable. Prefer a small reuse or two-buffer
   change when it is measurably faster, but do not add custom allocators, pools,
   generic workspace frameworks, or extensive plumbing merely to reduce
   allocations.

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

The numerical modes are verified candidate search, unsafe candidate search, and
exact candidate search. No flag selects the rigorously one-sided verified
procedure; `--unsafe` selects the faster heuristic; `--exact` bypasses both.
`--exact` and `--unsafe` are mutually exclusive. If build-time or runtime
floating-point requirements are unavailable, default mode stops before support
enumeration and requires an explicit `--exact` or `--unsafe` choice. An
inconclusive proof for an individual support still falls back to exact
arithmetic normally.

## Computation Flow

```text
CLI or pybind input
  -> matrix_parser
  -> fracessa constructor
  -> support generation and pruning
  -> find_candidate_verified, find_candidate_unsafe, or exact bypass
  -> find_candidate_exact
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
- `find_candidate_verified` lazily translates and positively scales the exact
  game into `[-1,1]`, encloses every entry in binary64, and reuses one bordered
  LU scratch matrix. It rejects only after a rigorous solution-error bound
  proves a nonpositive support probability or a profitable outside strategy.
- `find_candidate_unsafe` owns its normalized binary64 matrix and bordered
  scratch, then applies the existing danger veto. Constant games and unusable
  conversions select exact candidate solving.
- `find_candidate_exact` owns the reusable exact bordered system and constructs
  the authoritative candidate in `candidate_`.
- `fracessa` owns one exact game. All three concrete candidate classes store a
  reference to it, so no game matrix is copied between procedures.
- Exact arithmetic uses FLINT `fmpq_t` through `linalg::fraction`.
- Stability uses exact rational positive-definiteness; a binary64 result is not
  accepted as a final mathematical certificate.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md` documents the concrete failures and
  proves why tolerance tuning cannot recover an exact PD certificate.
- `--exact` does not initialize or allocate either double candidate-search
  state, and all final stability decisions use exact rational arithmetic.
- `--unsafe` uses a cheap pivot/margin danger veto, not a proof. IDs 45-47 are
  active regressions for verified search and remain counterexamples
  to the unsafe heuristic.

Key files:

- `cpp/include/fracessa/bitset64.hpp`: support-mask primitives.
- `cpp/include/fracessa/supports.hpp`: support generation and pruning.
- `cpp/include/linalg/fraction.hpp`: FLINT rational wrapper.
- `cpp/include/linalg/linear_solver.hpp`: exact bordered-system solver.
- `cpp/include/linalg/copositive_fraction.hpp`: exact copositivity checks.
- `cpp/include/fracessa/find_candidate_verified.hpp` and
  `cpp/src/find_candidate_verified.cpp`: verified class, strict one-sided proof,
  availability check, and focused proof-helper contracts.
- `cpp/include/fracessa/find_candidate_unsafe.hpp` and
  `cpp/src/find_candidate_unsafe.cpp`: unsafe class, normalized game, heuristic
  solve, and reusable scratch.
- `cpp/include/fracessa/find_candidate_exact.hpp` and
  `cpp/src/find_candidate_exact.cpp`: exact class, bordered system, and candidate
  construction.
- `cpp/include/fracessa/fracessa.hpp` and `cpp/src/fracessa.cpp`: exact game
  ownership, mode coordination, support search, and candidate lifecycle.
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
The `find_candidate_verified` object target overrides normal throughput flags with
strict floating-point semantics, contraction disabled, and IPO/LTO disabled.
One centralized availability function combines compiler support with the
runtime binary64, round-to-nearest, and subnormal checks. An unsupported build
still provides exact and unsafe modes, but refuses verified search.
Do not run verification IDs 33 or 34 with `--exact` without Reinhard's
separate approval.
When sandboxing blocks the normal ccache directory, rerun the build with
escalated filesystem access rather than disabling or redirecting ccache.

`cmake --install cpp/build` installs the CLI target only. It does not install
GMP, MPFR, or FLINT.

## Tests And Test Data

```bash
./test.sh # build, core/CLI tests, and wrapper tests
```

The non-matrix CTest suite consists of ten GoogleTest executables plus one CLI
black-box parser test. Wrapper tests use Python `unittest`. Matrix correctness
is no longer wired as one CTest per matrix.

`testdata/fracessa_testdata.sqlite3` is the canonical matrix, expected-result,
and timing store. Its strict schema is in `testdata/schema.sql`; the current
snapshot has 87 matrices and 49,157 stored candidate representatives. Nullable
multipliers recover weighted totals of 86,152 candidates and 83,377 ESS:
circular rows store one smallest dihedral representative and its orbit count,
while non-circular rows store null. Candidate IDs and row order remain
reproducibility checks; complete weighted candidate sets and ESS
classifications are the mathematical contracts.

Every dimension from 2 through 25 has at least one circular and one
non-circular matrix. IDs 67-79 fill the previously missing combinations with
deterministic random integers; exact and unsafe runs agreed on their complete
rational candidate contracts before insertion.

IDs 45-47 preserve the verified-search regressions: the dimension-20 unsafe
heuristic counterexample, the LU-boundary fallback case, and the failed-proof
exact-fallback case. IDs 48-55 cover non-circular dimensions 15-24 through
Hilbert, Hadamard, Paley conference, MINIJ, Fiedler, deterministic random
families, and a dense weighted-Laplacian game with one full-support ESS.

Every distinct matrix in Tables 1 and 2 of the Bomze-Schachinger-Ullrich
ESS-growth paper is present exactly once. IDs 18 and 26 are the exact published
Table 1 matrices for dimensions 12 and 17. IDs 80-81 are the two previously
missing Table 2 circular base matrices, and IDs 82-90 are its nine constructed
non-circular matrices. Same-property alternatives formerly stored at IDs 12
and 21 were removed; the former contents of IDs 18 and 26 were replaced by the
published vectors.

The timing snapshot has one CPU-2 session with a one-second target and 348
persistent-Pybind median rows. Current unsafe, verified, and exact use one
Release/native/LTO binary at algorithm revision `34e003168607`; historical
default, very unsafe uses raw-double algorithm revision `32f61679da64` with a
temporary normal-parser/nanosecond Pybind adapter. Every mode covers all 87
matrices. Historical very unsafe mismatches IDs 38-39, current unsafe mismatches
regression IDs 45-47, and current verified and exact match all 87. Timing
reports include matrix dimension, circularity, and the derived paper-style
lower bound `gamma_lower_bound = expected_ess ** (1 / dimension)` without
storing it in SQLite.

The former JSON/CSV verification, baseline-generation, speed-benchmark, and
Callgrind runners were removed. There is no replacement matrix-verification
runner yet. The small `python -m fracessa.timing` tool reads matrices from
SQLite, measures one build and one matrix at a time on a user-selected Linux
CPU, and writes normalized nanosecond samples to the same database. Reusing a
session name groups separately invoked builds. Each row records `source_ref`
(for a moving name such as `main`), its immutable `revision`, the binary hash,
backend, mode, CPU, comment, observed ESS count, target and measured wall time,
iteration count, and median native elapsed time. One Pybind process stays open
for every selected mode and matrix in a build. The first returned native
duration chooses `ceil(target / duration)` total samples and remains in the
sample; a duration at or above the target chooses one run. The stored result is
the median returned `elapsed_ns`. Measured wall time is metadata only and
never chooses or determines the reported timing. The CLI backend remains
available for legacy inspection but starts a child process per sample and must
not be mixed with persistent-Pybind microbenchmarks. Dated material under
`experiments/` and `aidocs/experiments/` remains immutable historical
evidence.

Database IDs 45-47 preserve the known unsafe-filter correctness regressions
tracked in `reviews/CPP_REVIEW.md`. The wrapper integration suite checks ID 46
through both verified and unsafe routes, but no complete SQLite matrix suite is
currently wired into `./test.sh` or release CI.

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

`RunConfig()` selects verified candidate search by default. `unsafe=True`
selects the heuristic that can miss candidates and ESS results, while
`exact=True` bypasses both numerical procedures. Exact and unsafe together are
rejected by the native boundary.

`run` and `run_multiprocessing` are the only public execution functions. Both
accept one `Matrix` or an iterable and accept an optional sink. One matrix
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
