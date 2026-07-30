# Project Knowledge

Last verified: 2026-07-30

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

## Priorities

1. Correctness is absolute; speed is second. Other concerns are secondary.
2. FracESSA evaluates a `2^n` support space and performs exact fractional
   operations millions or billions of times. It is optimized for many
   operations on small matrices, not a few operations on large matrices.
3. Keep validation at input boundaries. Do not add checks, allocations,
   abstractions, or branches to proven hot paths without a demonstrated
   correctness need and a benchmark.
4. Intentionally unchecked bitset operations exist for raw speed. Matrix input
   instead has one validating parser at the input boundary; values up to 18
   digits use direct integer construction and larger values use exact FLINT text
   conversion.
5. Use the Ponytail skill for code work: understand the complete path, then use
   the smallest correct implementation.

## Repository

- Root: `/home/reinhard/projects/fracessa`
- Remote: `git@github.com:reinhardullrich/FracESSA.git`
- Main implementation: `cpp/`
- Python API and verification: `python/`
- Reproducible benchmark material: `experiments/`
- Agent documentation: `aidocs/`
- Public GitHub introduction: `README.md`
- `AGENTS.md` must remain a pointer only.

Generated or local-only paths include `cpp/build*/`, `python/results/`, raw
Callgrind output, and experiment `builds/`, `sources/`, and `logs/` directories.

## Product Surface

FracESSA is a C++17 ESS (evolutionarily stable strategy) analyzer for symmetric
payoff matrices. It has two entry surfaces:

- CLI: `cpp/src/main.cpp`, built as `cpp/build/fracessa`.
- Native Python module: `cpp/src/pybind_module.cpp`, exposed through
  `python/wrapper_v1/` as `fracessa_core`.

CLI matrix format is `dimension#values`. Values are either the upper triangle
of a symmetric matrix (`n*(n+1)/2` entries) or the compact circular-symmetric
form (`floor(n/2)` entries).

The parser accepts dimensions 1 through 63, rejects textual fractions with a
zero denominator, and validates value syntax during the same scan that builds
the exact fractions. Support masks have 64 storage bits, but complete
enumeration requires the exclusive `2^n` bound, so dimension 64 is not
supported.

The temporary numerical modes are unsafe filtering and exact solving. No flag
and `--unsafe` select the same uncertified filter; `--exact` bypasses it.
`--exact` and `--unsafe` are mutually exclusive.

## Computation Flow

```text
CLI or pybind input
  -> matrix_parser
  -> fracessa constructor
  -> support generation and pruning
  -> find_candidate_dbl
  -> find_candidate_frc
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
- In non-exact mode, `MatrixServer` exactly translates and scales the game into
  `[-1,1]`, converts it once to binary64, and reuses one bordered scratch matrix.
  Constant games and unusable zero/subnormal conversions permanently fall back
  to exact candidate solving before enumeration.
- Exact arithmetic uses FLINT `fmpq_t` through `linalg::fraction`.
- Stability uses exact rational positive-definiteness; a binary64 result is not
  accepted as a final mathematical certificate.
- `correctness/DOUBLE_PD_FALSE_POSITIVES.md` documents the concrete failures and
  proves why tolerance tuning cannot recover an exact PD certificate.
- `--exact` does not initialize or allocate the double candidate-filter state,
  and all final stability decisions use exact rational arithmetic.
- The default unsafe filter uses a cheap pivot/margin danger veto, not a proof.
  It can still reject valid exact candidates; see `reviews/CPP_REVIEW.md`.

Key files:

- `cpp/include/fracessa/bitset64.hpp`: support-mask primitives.
- `cpp/include/fracessa/supports.hpp`: support generation and pruning.
- `cpp/include/fracessa/matrix_server.hpp`: reusable matrix construction.
- `cpp/include/linalg/fraction.hpp`: FLINT rational wrapper.
- `cpp/include/linalg/linear_solver.hpp`: exact bordered-system solver.
- `cpp/include/linalg/copositive_fraction.hpp`: exact copositivity checks.
- `cpp/src/findeq.cpp`: candidate construction.
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
3 development headers, GMP, MPFR, and FLINT. CMake FetchContent downloads:

- `spdlog`: optional rotating diagnostic logs.
- `argparse`: the cross-platform CLI parser.
- `googletest`: nine C++ unit-test executables only; it is not linked into the
  production executable.
- `pybind11`: the native Python module.

These four dependencies are currently fetched unconditionally, so a clean
configure needs network access unless the FetchContent sources are cached.
CMake also builds the tests unconditionally; `BUILD_TESTING=OFF` is not wired
up yet.

Local non-MSVC builds default to `FRACESSA_NATIVE_ARCH=ON` (`-march=native`).
Release CI sets it to `OFF`. IPO/LTO is enabled only when CMake confirms support.
When sandboxing blocks the normal ccache directory, rerun the build with
escalated filesystem access rather than disabling or redirecting ccache.

`cmake --install cpp/build` installs the CLI target only. It does not install
GMP, MPFR, or FLINT.

## Tests And Benchmarks

```bash
./test.sh          # core/CLI, wrapper, and fast matrix correctness
./test.sh --full   # core/CLI, wrapper, and all matrix correctness
./python.sh        # fast speed benchmark
./python.sh --full # all-matrix speed benchmark
```

The non-matrix CTest suite consists of nine GoogleTest executables plus one CLI
black-box parser test. Wrapper tests use Python `unittest`. Matrix correctness
is one CTest per matrix so CTest can run matrices in parallel.

Fast mode is a static policy: all verification matrices except IDs 32 and 34.
The speed script records timings; correctness belongs to CTest.

All 52 active verification matrices (IDs 1-44 and 48-55) pass in the temporary
default unsafe mode; IDs 38-39 specifically cover the corrected scale and
translation behavior. Preserved reference IDs 45-47 are intentionally not
active yet: they prove the unsafe heuristic can still miss an exact ESS and
belong to the later certified Choice 1 phase. IDs 48-55 add non-circular
dimensions 15-24 through Hilbert, Hadamard, Paley conference, MINIJ, Fiedler,
deterministic random matrix families, and a dense weighted-Laplacian game with
one full-support ESS.

## Verification Data

`python/verification/` is the only active verification-data source. There is no
second fixture copy under `cpp/tests/`.

- `verification_matrices.json`: 52 matrices, active IDs 1-44 and 48-55,
  maximum dimension 24.
- `baseline_candidates.csv`: 1,196 candidate representative rows for all 52
  active matrices, representing 38,044 candidates after multiplication.
- `baseline_result.json`: historical timing results for IDs 1-35.
- `ctest_verify_matrix.py`: matrix correctness comparison.
- `matrix_selection.py`: fast/full static selection.
- `create_baselines.py`: baseline regeneration using the archived executable.

`testdata/fracessa_testdata.sqlite3` is a staged SQLite migration snapshot, not
an active test input yet. It currently contains 63 matrices and 29,114 stored
candidate representatives. Their multipliers represent 65,962 candidates and
63,369 ESS. Its `matrices` rows use stable IDs rather than names and include
dimension, size class, circular symmetry, exact input, weighted candidate/ESS
counts, weighted support-size structures, and required provenance/purpose text
in `origin`. Qualitative categories live in the `tags` JSON array; quantitative
facts are not duplicated as tags. IDs 56-66 are staged non-circular
complete-multipartite matrices with many ESS for support-frontier benchmarks.
The schema is in `testdata/schema.sql`; benchmark-run storage remains deferred.
Do not switch Python or CTest consumers from the existing JSON/CSV files without
separate approval.

There is no `in_use` field. Candidate IDs and row order are deterministic and
remain useful correctness checks while enumeration is unchanged. Circular rows
store the smallest dihedral support representative and a non-null `multiplier`;
non-circular rows store null. Weighted candidate and ESS totals are recovered by
summing `multiplier`, treating null as one. The DFS/FKM generator adoption and
representative-output conversion were checked across all 52 active matrices and
regenerated the affected circular fixtures in both the CSV baseline and SQLite
snapshot.
`T_pd_dbl` and `T_pd_frc` are normalized as the same `T_pd_*` classification for
baseline comparison because floating-point branch selection can change while
the mathematical result does not.

The archived baseline executable is x86-64. Do not run it on this ARM64 Linux
machine: it invokes the system `binfmt` dispatcher without a configured x86-64
emulator. Use native CTest correctness checks for current development.

## Pybind Boundary

`cpp/src/pybind_module.cpp` exposes the C++ analyzer as the native
`fracessa_core` module. It owns Python/native argument and result conversion,
native status codes, GIL release, and native timing. Binding-specific open
findings are tracked in `reviews/PYBIND_REVIEW.md`, separately from both the
analyzer core and Python orchestration.

## Python Wrapper

`python/wrapper_v1/` calls `fracessa_core` in-process and is the maintained API.
It supports sequential execution, process-based parallelism across matrices,
and stream/CSV/JSON/Arrow/Parquet sinks. One matrix is always computed by one
worker process; parallelism is across matrices.

For circular matrices, `CandidateRow` contains one bracelet representative with
an integer `multiplier`; non-circular rows use `None`. `SummaryRow.candidate_count`
is the number of returned representative rows, while `ess_count` remains the
weighted mathematical total.

No production wrapper or matrix workflow imposes a per-matrix computation
timeout. A matrix may legitimately run for hours. Worker-liveness handling must
not be implemented as a computation timeout.

`run_jobs_mp` uses bounded streaming submission. The lower-level
`MPQueueRunner` still has backpressure, dead-worker, and shutdown problems; see
`reviews/PYTHON_REVIEW.md` before using it for unbounded generated input.

`python/fracessa_py.py` is legacy subprocess code still used by the speed script.
New API work belongs in `python/wrapper_v1/` and `fracessa_core`.

## Release Workflow

`.github/workflows/release.yml` runs only for pushed `v*` tags. It builds and
checks Ubuntu, macOS, and Windows before publishing three executables.

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
