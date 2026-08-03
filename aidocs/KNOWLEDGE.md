# Project Knowledge

Last verified: 2026-08-03

## Worktree Policy

1. Work in the main worktree at `/home/reinhard/projects/fracessa` unless Reinhard explicitly approves another worktree first.
2. Do not switch to, run project commands in, or modify another worktree without that prior approval.

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

`research/papers/` contains the retained PDFs and complete audited Markdown transcriptions of all five papers. The 2014, 2015,
and 2018 Bomze-Schachinger-Ullrich papers also retain the exact LaTeX sources that generated their PDFs and seven rendered figure
assets used by the Markdown versions.

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
every candidate decision. `fast` removes the game's common denominator, applies an exact integer precision-span check, converts
and equilibrates one symmetric triangle of the complete binary64 game once, mirrors it, and solves each border-eliminated system with
Bunch-Kaufman $LDL^T$. A matrix with $P\geq10^9$ switches entirely to `safe`, and a support with an inconclusive pivot below
$10^{-12}$ reaches exact checking. Its remaining probability and outside-payoff rejections are heuristic. Experimental `test` is
an independent source copy of `fast`, not a wrapper around it; it currently has identical behavior and can be changed without
changing production fast search.

Every result exposes `safe_fallback`: null means the selected fast/test method reached its double search, while `precision_span`,
`equilibration_invalid`, or `equilibration_non_convergence` identifies the whole-matrix preparation step that switched the run to
safe search. A per-support exact retry does not set it. CLI timing output prints this as line 3; Pybind and all Python sinks use the
same nullable field. Historical timing rows may retain the legacy generic value `equilibration`.

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
  production `CircularSupportGeneratorV3` uses the direct binary fixed-density
  bracelet recursion of Karim-Alamgir-Husnine. Both prune partial branches against earlier exact
  candidate supports bucketed by lowest bit. Circular rules expand every
  distinct rotation/reflection only as compact forbidden masks. The analyzer
  stores one solved bracelet representative with their count as `multiplier`.
  The former FKM-plus-reflection implementation remains available as `CircularSupportGenerator` (V1), and compact
  bit-parallel `CircularSupportGeneratorV2` also remains test-only.
  See `plans/SUPPORT_GENERATORS.md`.
- Newly found exact candidates are pending until the next cardinality, keeping
  each generator layer's pruning rules stable. Stability is irrelevant to this
  pruning rule: every exact equilibrium support forbids later strict supersets.
- `find_candidate_fast` reuses the safe solver's denominator-cleared integer game. It switches the complete matrix to safe search
  when the remaining integer entries or distinct gaps give $P\geq10^9$; otherwise one common power-of-two normalization precedes
  one LAPACK-derived symmetric BIN equilibration $A\mapsto DAD$ of the complete game. Conversion and final equilibration scaling
  each process one triangle and mirror it. Each support eliminates the
  payoff/normalization border, forms the transformed reduced Hessian, and solves it with the adapted lower-triangle Bunch-Kaufman
  $LDL^T$ factorization and solve. The transformed normalization vector and row scales recover probabilities and payoffs in the
  original game coordinates. Exact all-zero rows and columns retain scale one and are excluded individually from the active BIN
  iteration; every failure or non-convergence in the remaining nonzero block, every small or non-finite pivot, and every non-finite
  completed probability or payoff accumulation are inconclusive and reach exact checking. Outside strategies are enumerated
  directly from the complement mask only after factorization and probability checks pass. Probability and outside-payoff
  rejections remain heuristic, so `fast` is not a correctness certificate.
- `find_candidate_test` intentionally duplicates fast double storage, equilibration, and reduced-system kernel without sharing its
  implementation or double state. Its source is currently identical to fast apart from class and header names.
- `find_candidate_safe` clears the game's rational denominators once, eliminates
  the normalization/payoff border, and constructs the integer-scaled symmetric
  reduced system $dH y=dr$ in reusable FLINT storage. One fraction-free
  $LDL^T$-style factorization solves the candidate, proves singularity, and
  records the exact inertia needed by stability. Rational values are constructed
  only for successful public candidate output.
- `fracessa` owns the rational game used by stability. `find_candidate_safe` owns the one integer-scaled copy used by exact
  candidate solves and both one-time precision-span decisions; fast and test each own the converted and equilibrated double copy
  prepared from that integer game.
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
  `cpp/src/find_candidate_fast.cpp`: production exact precision-span gate, whole-game double equilibration, reduced symmetric
  Bunch-Kaufman solve, small-pivot fallback, heuristic inequalities, and reusable scratch.
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

### Canonical Performance Build

Unless the user explicitly requests a compiler experiment, every stored performance benchmark and every build compared with it
must use this configuration:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S cpp -B cpp/build-benchmark -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DFRACESSA_NATIVE_ARCH=ON
cmake --build cpp/build-benchmark --target fracessa fracessa_core --parallel
```

The CMake configuration is the source of truth; do not reproduce or amend its flags by hand. On the current Fedora/GCC toolchain
it produces `-O3 -DNDEBUG -funroll-loops -march=native` plus IPO/LTO (`-flto=auto -fno-fat-lto-objects`). A configure where CMake's
IPO check fails is not a canonical benchmark build. Do not add floating-point, aliasing, sanitizer, profiling, debug, or alternate
optimization flags. In particular, fast, safe, and test are runtime methods in the same `fracessa_core` binary and must never be
compiled with method-specific settings.

Both sides of a comparison must be clean revisions built afresh with the same compiler executable and version, CMake version,
generator, dependency versions, and command above. The current reference toolchain is GCC 16.1.1, CMake 4.3.0, and Ninja 1.13.2.
A toolchain upgrade starts a new benchmark lineage and must be recorded in the timing comment; it must not be presented as a
source-only comparison with older rows. Canonical timing uses CPU 2, one persistent Pybind process, matrix-owned fast/safe
calibrations, a 0.5-second default target, and the median native nanosecond duration. An explicitly requested compiler experiment
must be labelled as an experiment, record every deviation, and must not replace canonical rows.

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
and timing store. Its strict schema is in `testdata/schema.sql`; the current snapshot has 1,064 distinct strategically normalized
matrices. The 713 analyzed rows have 65,800 stored candidate representatives; nullable multipliers recover weighted
totals of 104,098 candidates and 91,134 ESS:
circular rows store one smallest dihedral representative and its orbit count,
while non-circular rows store null. Candidate IDs and row order remain
reproducibility checks; complete weighted candidate sets and ESS
classifications are the mathematical contracts. Of the 351 catalog-only rows, 206 exceeded the two-minute safe retry cutoff and
145 generator-catalogue rows exceeded their initial one-second safe cutoff. They keep all four baseline-summary fields null as one
group; null never means zero candidates or zero ESS. SQLite enforces stored-input uniqueness on `(dimension, matrix)`;
`matrix` alone cannot be unique because compact input omits its dimension. Import audits additionally reject matrices whose stored
value vectors have the same dimension and circular-storage flag and differ only by a positive nonzero rational multiplier; they
retain the lowest matrix ID. Negative multipliers remain distinct because they reverse payoff comparisons.

Corpus sampling retains exactly three diagonal matrices: dimension-one ID 314, compact all-zero dimension-50 ID 2203, and
nonzero dimension-60 ID 2180. The other 27 diagonal rows were removed to avoid overrepresenting this structurally trivial family.
This is a benchmark-corpus sampling choice, not a claim that distinct diagonal games are mathematically equivalent.

Matrix size is represented only by `size_class`: `small` covers dimensions 1-8, `medium` 9-16, `large` 17-25, and
`super_large` 26-63. Do not duplicate this fact in `tags`.

For optimization-problem collections, import only explicitly stored symmetric
objective matrices. Do not add matrices whose only role is as a constraint. If
such a matrix already exists from an independent source, retain the existing
row rather than treating the constraint occurrence as grounds for deletion.

Matrix rows also carry human and machine-readable catalog metadata: name,
family, subfamily, description, source URL, original format and ID, and a
creation/first-use date. `origin` remains the established prose provenance
field, while `source_url` stores an external website or DOI. Legacy REF
fixtures use `family = "historical"`; when their exact day is unknown, January
1 of the known year is a documented placeholder rather than an asserted
historical day.

IDs 91-117 originally held the 27 distinct SuiteSparse Matrix Collection imports selected
from the official 2,904-entry index by square shape, exact numerical symmetry,
real values, and dimension at most 63. Matrix Market pattern, integer, and
finite printed real tokens are converted to exact rational values in FracESSA format. This
preserves the downloaded decimal tokens but does not certify an unprinted
physical value as symbolically rational. Exact dense rational duplicates are
skipped; SuiteSparse ID 2758 duplicates matrix ID 1. Imported dimensions 26-63
use `size_class = "super_large"`. Diagonal sampling later removed ID 105 (`HB/bcsstm01`), leaving 26 SuiteSparse rows.

The finite NIST Matrix Market search across all symmetric categories and
dimensions 1-63 produced nine downloadable exact-symmetric files. Eight exactly
duplicate SuiteSparse rows; `BFW62B` is the same source and sparsity pattern as
SuiteSparse `Bai/bfwb62` with eight last-digit differences no larger than
`1e-21`; and NIST withdrew its incorrect `LAP 25` Matrix Market assembly. The
earlier `FIDAP005` check likewise identifies a rounded `FIDAP/ex5` source copy.
NIST pages remain alternate provenance, not separate matrix rows.

The QAPLIB audit selected 152 exact distinct matrices from the official `qapdata` archive at IDs 118-269. Of 136 fixed instances,
109 have dimensions at most 63; their A and B integer matrices are tested independently, yielding 182 symmetric occurrences and
30 internal duplicates. The later zero-game consolidation removed circular zero matrix ID 158, leaving 151 QAPLIB rows. They cite
DOI `10.7488/ds/3428`; redistributed data is attributed to Burkard, Cela, Karisch, Rendl, Anjos, and Hahn under CC BY 4.0.

IDs 270-280 are the 11 exact symmetric integer matrices from the official
TSPLIB95 symmetric-TSP archive that declare `EDGE_WEIGHT_TYPE: EXPLICIT` and
have dimension at most 63. The three retained explicit layouts were expanded
exactly and checked edge-for-edge against the independent official TSPLIB95 XML
files. The six other explicit symmetric TSPs exceed dimension 63; coordinate-
derived instances, tours, and asymmetric problem categories are excluded. The
11 imports are non-circular, globally distinct, and catalog-only.

IDs 281-313 are the 33 in-range exact integer matrices from the official Biq
Mac Library archive: 10 Beasley and 13 Glover-Kochenberger-Alidaee
binary-quadratic matrices plus 10 Rudy Max-Cut graphs. The full archive has 468
files representing 343 logical instances; all 125 dense/sparse pairs match
exactly, and 310 logical instances exceed dimension 63. Sparse Q matrices use
the library's symmetry contract, while Max-Cut edge lists are stored as
symmetric zero-diagonal adjacency matrices. The retained rows are non-circular,
globally distinct, and catalog-only.

IDs 314-319 are the six exact symmetric in-range representatives from Magma's
Hadamard database: representative 1 at degrees 1, 2, 4, 8, 16, and 32. The
compact binary decoding matches the exact degree-16 example printed in Magma's
handbook, and every one of the 4,474 ordinary representatives through degree 63
was checked for the Hadamard identity and exact symmetry. The separate 638-row
skew-Hadamard database was also audited; all rows have degrees 36, 44, or 52
and none is symmetric. The six source matrices are globally distinct and catalog-only. Dimension-one ID 314 is mathematically
circulant but remains in full non-circular storage because the compact representation would contain no values; IDs 315-319 are
non-circular.

Thirty of the 35 in-range QPLIB problems have an explicitly stored quadratic
objective. The database retains only these 30 objective Hessians; quadratic-
constraint matrices and linear objective vectors are excluded. Each stored
lower triangle is mirrored exactly, and finite printed coefficients are parsed
as exact fractions. The objectives are pairwise distinct, duplicate no earlier
catalog row, and retain their original noncontiguous IDs between 320 and 1925.
An independent PyQPLIB 0.1.7 parse matches all coefficients and roles. The
QPLIB imports are non-circular, catalog-only, and attributed under CC BY 4.0
with DOI `10.1007/s12532-018-0147-4`.

IDs 1926-1982 are 57 exact distinct matrices from every locally maintained
problem family in J.E. Beasley's current OR-Library index plus its still-hosted
urban-transit page: 23 binary-quadratic Q matrices, 23 capacitated minimum-
spanning-tree cost matrices, six aircraft-
separation matrices, two CAB hub-location matrices, one portfolio correlation
matrix, and two urban-transit demand matrices. Eligibility requires an
explicitly stored rational square matrix, exact symmetry, and dimension at most
63; coordinate-derived, shortest-path-derived, rectangular, and externally
linked data are excluded. The OR-Library Q matrices use the source's
maximization sign and are therefore negatives of, not duplicates of, the Biq
Mac minimization copies. Exact deduplication removes ten repeated `capmst`
occurrences and the repeated `portreb1` correlation matrix. Two in-range
`capmstnew` matrices, two aircraft tables, and six corporate tax tables fail
symmetry. The archived official `td1` and `td2` demand files are retained;
transit time files contain nonnumeric absent-link markers. Every import is
non-circular, globally new, catalog-only, and covered by OR-Library's MIT
license.

COMPl_e_ib 1.1 was audited as a control-system benchmark library rather than
treated as a generic bag of arrays. It defines 168 state matrices `A`: 57 exceed
the dimension-63 limit, and none of the 111 in-range matrices is exactly
symmetric. Input/output channels and synthesized identity, zero, and weighting
arrays are not independent benchmark matrices. The audit therefore adds no
catalog row.

The SLICOT model-reduction collection has 18 linear-system benchmarks. Seventeen
exceed dimension 63; the sole in-range state matrix is the order-48 building
model, whose `A` matrix is asymmetric. Its `build.mat` file is byte-identical
to COMPl_e_ib's `lah.mat`, so this audit also adds no catalog row.

IDs 1983-1990 are eight exact adjacency matrices from KONECT's 23 downloadable
unipartite networks through dimension 63. All nine undirected files are
symmetric; directed `moreno_taro` is also exactly symmetric because its arcs
occur in reciprocal pairs, while the other 13 directed matrices are
asymmetric. Seven small bipartite files are excluded because their native
tables are rectangular. Exact dense deduplication maps KONECT Dolphins and
Zachary karate club to existing SuiteSparse IDs 114-115 as alternate
provenance. Every new KONECT row is non-circular and catalog-only. KONECT does
not state one collection-wide dataset license, so project metadata makes no
broader licensing claim.

IDs 1991-2137 began as a deterministic 147-graph stratified sample from all 23,988 House of Graphs entries with dimensions
1-63 on August 2, 2026. Five dimension bands are crossed with acyclic, connected bipartite cyclic, connected planar
non-bipartite, connected nonplanar, disconnected cyclic, regular, dense, vertex-transitive, asymmetric, and unrestricted
control categories. Up to three graphs per populated stratum are chosen by SHA-256 rank with seed `20260802`; one stratum
has only one eligible graph and one graph occurs in two selected strata. Every canonical graph6 value matched the API
adjacency list, and no imported matrix duplicated an existing dense matrix. The later zero-game consolidation removed empty graphs
1992, 2080, and 2127, leaving 144 rows, including five compact circular matrices. No collection-wide House of Graphs data license
was found, so metadata makes no broader claim.

IDs 2138-2176 are the 39 exact, globally new Network Data Repository matrices selected from 1,241 index rows reporting
dimensions 1-63. The deterministic sample crosses the repository category with dimension bands `1-8`, `9-16`, `17-25`,
`26-44`, and `45-63`, ranks with SHA-256 seed `20260802`, and retains at most three directly represented symmetric matrices
per populated stratum after deduplication against the 558-row pre-import database. The rows comprise 15 animal-social, 15
cheminformatics, six protein, two DIMACS, and one biological matrix. Exact Matrix Market expansion and undirected or
reciprocal edge-list reconstruction are accepted; rectangular or unsymmetric matrices, temporal streams, bipartite tables
without a square adjacency matrix, malformed files, and forced symmetrization are rejected. All 39 source archives
round-trip exactly, one row is circular, and all are catalog-only. The repository states a Creative Commons
Attribution-ShareAlike license without identifying a version.

IDs 2177-2206 originally held the 30 exact SDPLIB `F0` objective coefficient matrices whose complete block-diagonal dimensions are at
most 63. The 92-problem SDPLIB 1.2 mirror has exactly 30 such eligible problems. Only matrix number zero is retained from
each; all 1,799 `F1...Fm` constraint matrices and separate source blocks are excluded. Exact decimal and scientific tokens
become reduced fractions. The source objectives are pairwise distinct, duplicate no earlier source matrix, and remain catalog-only;
15 use `size_class = "super_large"`. Diagonal sampling later removed control IDs 2177-2179 and truss IDs 2204-2206, leaving 24
SDPLIB rows: one control, 15 H-infinity, two infeasible-dual, two infeasible-primal, three quadratic-assignment, and one theta
matrix. `theta1` is circulant with source diagonal `1`; retained ID 2203 stores the strategically equivalent compact zero-diagonal
matrix obtained by subtracting `1` from every entry, and its description records that normalization. The current GitHub mirror
declares GPL-3.0; project metadata makes no broader licensing claim.

IDs 2210-2687 contain 450 retained matrices from the combined Anymatrix, TypedMatrices.jl, and Matrix Depot exact-generator
audit. Treat the three repositories as one overlapping catalogue. The audited pool contains 2,678 distinct exact symmetric
matrices from 51 eligible families and 34 property classifications. SHA-256 seed `20260802` selects up to three representatives
per populated property-by-band stratum over `1-8`, `9-16`, `17-25`, `26-44`, and `45-63`, then adds any family still lacking a
representative. The final selection has 484 distinct matrices in 166 populated strata; six exact matches at existing IDs 48, 49,
314, 1995, 2001, and 2155 are tagged rather than duplicated. The raw import added 478 rows. Exact circular normalization of 53
full matrices exposed ID 2215 as an exact duplicate, and the later positive-scale audit removed IDs 2212-2214, 2220, and 2254.
The dimension-one consolidation then removed IDs 2210 and 2211 in favor of older ID 314. Diagonal sampling removed all 20 Strakos
rows. The retained new rows include 78 compact circular matrices and no dimension-one or Strakos row. One-second safe calibration
supplied complete baselines for 305 retained rows; the other 145 retain null baselines. The exhaustive source, exclusion, sampling, normalization, and
validation record is `aidocs/reference/MATRIX_GENERATOR_CATALOGUE_AUDIT.md`.

Every dimension from 2 through 25 has at least one circular and one
non-circular matrix. IDs 67-79 fill the previously missing combinations with
deterministic random integers; exact and the then-current normalized heuristic agreed on their complete rational candidate
contracts before insertion.

IDs 45-47 preserve the verified-search regressions against the retired normalized heuristic: the dimension-20 heuristic
counterexample, the LU-boundary fallback case, and the failed-proof exact-fallback case. IDs 48-55 cover non-circular dimensions 15-24 through
Hilbert, Hadamard, Paley conference, MINIJ, Fiedler, deterministic random
families, and a dense weighted-Laplacian game with one full-support ESS.

IDs 2207-2209 are exact non-circular false-rejection regressions for the former fast rules covering all three per-support candidate
conditions. ID 2207 has a
nonsingular full-support ESS but produces a $7.5\times10^{-13}$ double pivot below the $10^{-12}$ cutoff. ID 2208 has an exact
positive probability $10^{-10}$ that the double solve computes as negative. ID 2209 has an outside payoff exactly $10^{-4}$ below
the equilibrium payoff that the double solve computes above its rejection margin. Current fast and test recover them through the
small-pivot or matrix-wide precision-span fallback; their fast/test/safe ESS counts are `1/1/1`, `1/1/1`, and `2/2/2`.

Every distinct matrix in Tables 1 and 2 of the Bomze-Schachinger-Ullrich
ESS-growth paper is present exactly once. IDs 18 and 26 are the exact published
Table 1 matrices for dimensions 12 and 17. IDs 80-81 are the two previously
missing Table 2 circular base matrices, and IDs 82-90 are its nine constructed
non-circular matrices. Same-property alternatives formerly stored at IDs 12
and 21 were removed; the former contents of IDs 18 and 26 were replaced by the
published vectors.

IDs 2688-2695 store the eight previously missing payoff games derived from the 2014 and 2015 Bomze-Schachinger-Ullrich
copositive constructions at dimensions 6, 7, 8, 9, 11, 12, 13, and 14. They use the primitive integer transformation
`A = (dJ - S) / g`, where `d` is the common diagonal and `g` is the positive gcd; on the simplex this turns each isolated
copositive zero into a global strict maximizer and hence an ESS. Exact safe analysis gives ESS counts `8, 14, 20, 30, 33, 60,
108, 192`, respectively. The papers state 18 zeroes at dimension 8 and 27 at dimension 9; the complete ESS search finds two and
three additional nonglobal strict local maxima. The order-15 `S15` construction was already present as ID 24 with 360 ESS.

The timing snapshot has 710 CPU-2 persistent-Pybind rows. It retains Werner fast and safe, the July 27 pre-refactor GitHub build
renamed `classic`, the paired safe builds immediately before and after the C++ FLINT-wrapper extraction, and four 81-matrix fast
panels for `equilibrated-fast`, rejected FP-S01, retained FP-S02, and combined FP-S02+FP-S03. The experimental test panels use
historical timing mode `fast` because the strict schema does not admit `test`; their build labels and comments identify the actual
native method. Classic is
revision `32f61679`; it predates both the support-generator redesign and exact $LDL^T$, and its raw-double result mismatches IDs
38-39 in the historical corpus. The retained Werner panels now have 81 fast rows and 66 safe rows. The paired builds each have 77 dimension-3-or-larger rows with one-second native
medians; IDs 65-66 and 90 were not attempted because retained evidence did not establish a sub-30-second run, while pre-wrapper
ID 47 was measured at 50.394 seconds and then removed. All paired ESS counts match. Build label, revision, and binary hash identify
every stored build. All four newer panels also match their expected ESS counts. The generated matrix column
`gamma_lower_bound = ess_count ** (1 / dimension)` is the paper's lower bound for $\gamma$ implied by the matrix and is null when
the exact ESS count is unknown. `matrix_overview` places it immediately after `ess_structure`; timing reports print the same value.

Exact circular-storage normalization uses `A - dJ`, never `A - dI`: subtracting the common diagonal value `d` from every entry
preserves candidates and ESS on the simplex and shifts every payoff by `-d`. `scripts/normalize_circular_matrices.py` performs the
audit with exact fractions and records `d` in each description. Retained IDs 1 and 2203 use compact normalized storage with
`d = 0, 1`, respectively. Normalized IDs 39 and 41 duplicated older ID 1 and were removed. The generator catalogue added 53
further `A - dJ` conversions; normalized ID 2215 exactly duplicated the later-removed ID 44. Positive-scale deduplication then
removed IDs 38 and 44 in favor of ID 1. Zero-game consolidation removed ID 43 in favor of the retained dimension-50 ID 2203, and
dimension-one consolidation left ID 314 as the sole full-storage exception. Eighteen stale `classic`/Werner rows for the old
matrices were replaced by a current-build fast/safe panel; after both cleanups it has four rows for IDs 1 and 2203.

Fast calibration covers all 1,072 matrices with 774 positive durations and 298 `-1` timeouts; safe calibration has 718 positive
durations and 354 timeouts. No calibration remains null.
The reusable
`scripts/calibrate_matrices.py` processes only null calibration fields by default;
either method also stores missing candidate data when it finishes within the cutoff. Its `fast|safe --retry-timeouts` form
processes only that method's `-1` rows, runs each matrix once, and leaves unsuccessful rows at `-1`; the retained two-minute pass
uses `--cutoff-seconds 120`. Repeating `--cpu ID` processes matrices concurrently with one child pinned to each selected core
while the main process serializes SQLite writes. A completed fast result is exact when `safe_fallback` is non-null; otherwise it
is explicitly an unverified fast result while `safe_calibration_ns` remains `-1`. These calibrations choose future iteration
counts and are not timing-history rows.

`scripts/calibrate_matrices.py audit` is the resumable full consistency-calibration pass. It dispatches matrices in
`dimension ASC, matrix_id ASC` order, runs fast before safe, and uses each positive stored calibration to request
`ceil(1 second / calibration)` native samples. Each method independently uses
`max(120 seconds, 1.2 * its stored calibration)`; null and `-1` use 120 seconds. A fast timeout skips safe and sets both
calibrations to `-1`; a non-null fast `safe_fallback` supplies both exact
calibrations without a second safe run. Complete candidate output is compared with stored data, missing data is filled, and
conflicting stored data is preserved. `calibration_timestamp` stores the latest audit time, while `calibration_comment` is an
human-readable, indented, append-only JSON history containing the assigned CPU, actual cutoffs, and all outcomes. Repeated
`--cpu ID` options run matrices concurrently while keeping each matrix's fast and safe calls on the same CPU. Each matrix commits
independently. The default
selects only rows without an audit timestamp; `--refresh-all` deliberately starts a new audit over every row.

The August 2 fast retry attempted all 319 previous timeouts with a 120-second per-matrix cutoff. Two rows completed before CPU 2
was reserved; the remaining 317 attempts used performance CPUs 3 through 9. Twenty-one matrices completed, adding 683
representative rows for 841 weighted candidates and 236 ESS. Two completed through exact fallback and 19 are unverified fast
results. The remaining 298 timeouts comprise 41 exact-fallback classifications and 257 unverified fast rows.

A historical `reduced-hessian-ldlt` benchmark measured 85 matrices; IDs 33-34
were not included in that run. Those rows are no longer in the canonical timing table. All 85 ESS counts matched. Against
`current-main` on the same matrices, summed exact
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

The paired wrapper benchmark compares clean adjacent revisions: raw FLINT-owning storage at `3547df5d` and the C++ `integer` and
`matrix_int` wrappers at `29799de8`. Both use identical Release/native/LTO flags, one persistent Pybind process per build, CPU 2,
and one-second native medians. Across 77 matrices, the wrapper build wins 59, ties 4, and loses 14. Its arithmetic mean and median
per-matrix changes are $-7.01\%$ and $-8.50\%$; summed medians fall from 73.083 to 68.499 seconds, a $6.27\%$ reduction.
Non-circular and circular median changes are $-11.34\%$ and $-6.52\%$ respectively.

This is a compiler/code-generation gain, not a mathematical or storage-copy change. The wrapper entry references are non-owning,
the small methods inline to the same FLINT operations, and the algorithm is unchanged. A reverse-order ten-second matrix-60 check
measured 400.201 milliseconds before and 332.415 milliseconds after. Perf recorded approximately $9.3\%$ fewer instructions,
$16.8\%$ fewer cycles, and $16.2\%$ fewer branches per solve in the wrapped Pybind build. Rebuilding both revisions with
`-fno-strict-aliasing` left the wrapper result unchanged at 331.810 milliseconds but improved the raw-pointer build to 369.745
milliseconds, reducing the gap from $16.9\%$ to $10.3\%$. This isolates alias optimization as one contributor; the remaining
difference is consistent with the changed inline/LTO code shape and layout, but the diagnostic does not attribute it to one source
expression. Catalog-only imports have no timing rows, and the timing runner excludes rows without baselines by default.

The former JSON/CSV verification, baseline-generation, speed-benchmark, and
Callgrind runners were removed. There is no replacement matrix-verification
runner yet. The small `python -m fracessa.timing` tool reads matrices from
SQLite, measures one build and one matrix at a time on a user-selected Linux
CPU, and writes normalized nanosecond samples to the same database. Reusing a
session name groups separately invoked builds. Each row records `source_ref`
(for a moving name such as `main`), its immutable `revision`, the binary hash,
backend, historical `mode` database value, nullable whole-matrix `safe_fallback`, CPU, comment, observed ESS count, target and measured wall time,
iteration count, and median native elapsed time. One Pybind process stays open
for every selected method and matrix in a build. New runs require `--method fast` or `--method safe`; adapters map those choices
onto old Pybind and CLI interfaces when benchmarking historical builds. Each matrix row owns nullable `fast_calibration_ns` and
`safe_calibration_ns` values. Positive values choose `ceil(target / calibration)` samples, while `-1` records a calibration killed
at its cutoff and chooses one sample; a missing value is a hard error. `scripts/calibrate_matrices.py fast|safe` fills only missing
values with a one-second default cutoff; it classifies and stores `matrices.safe_fallback` before the timed process so a timeout does
not lose the reason, and either method stores missing candidate results. Use
`scripts/calibrate_matrices.py fast|safe --retry-timeouts --cutoff-seconds 120` to retry only that method's `-1` rows once. Clear a field manually
before an intentional refresh. Calibrations are not timing-history rows. The default target is 0.5 seconds, and the stored result is the
median returned `elapsed_ns`. Measured wall time is metadata only and never chooses or determines the reported timing. The CLI
backend remains
available for legacy inspection but starts a child process per sample and must
not be mixed with persistent-Pybind microbenchmarks. Timing reports retain fallback rows but exclude non-null fast fallbacks from
speed-ratio summaries. Dated material under
`experiments/` and `aidocs/experiments/` remains immutable historical
evidence.

Database IDs 45-47 preserve the retired normalized-heuristic correctness regressions tracked in `reviews/CPP_REVIEW.md`, and IDs
2207-2209 preserve the three former fast per-support false rejections. The wrapper integration suite checks IDs 46 and 2208 through
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
clock and unit, followed by `safe_fallback` on line 3. There is no wrapper timing-suppression option.

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
