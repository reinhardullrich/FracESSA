# FracESSA Project Overview

Last verified: 2026-08-08

## Repository Map

- Root: `/home/reinhard/projects/fracessa`
- Remote: `git@github.com:reinhardullrich/FracESSA.git`
- C++ implementation and CLI: `cpp/`
- Python package: `python/pyfracessa/`
- Python tests: `python/tests/`
- Canonical test data: `testdata/fracessa_testdata.sqlite3`
- Agent documentation: `aidocs/`
- Public project page: `README.md` and `site/`; Sphinx sources for the combined Python/C++ API site are under `docs/`.
- Preserved predecessor source: `zzz_legacy/`; it is historical, not active source.
- Local-only material: `research/`, `experiments/`, and database-maintenance `scripts/` are ignored by Git.
- Generated/local paths include `cpp/build*/`, `.local/`, and `.local-tmp/`.

## Product Contract

FracESSA is a C++17 analyzer for evolutionarily stable strategies in symmetric payoff games. It has two maintained entry surfaces:

- CLI: `cpp/src/main.cpp`, normally built as `cpp/build/fracessa`.
- PyFracESSA: `python/pyfracessa/`, backed by the native `fracessa_core` module in `cpp/src/pybind_module.cpp`.

Input uses `dimension#values`. Values are either the upper triangle of a symmetric matrix or the compact circular-symmetric form.
The validating parser accepts every positive dimension whose dense matrix size is representable and exact integer or rational values.
Dimensions 1 through 64 use one `uint64_t`; dimensions 65 and above use fixed-width multiword support masks. Runtime remains
exponential even when a dimension is representable, so the practical limit is normally much smaller than the storage limit.

Every entry surface requires an explicit method; there is no default:

- `safe`: exact candidate and stability decisions.
- `fast`: binary64 candidate search with exact matrix-wide and per-support fallbacks; remaining probability and outside-payoff
  rejections are heuristic, so fast is not a correctness certificate.
- `test`: an independent experimental copy of fast, currently behaviorally identical but intentionally separable.

`safe_fallback` is null when fast/test used the double search. A non-null preparation reason means the whole matrix switched to safe
search. A per-support exact retry does not set it. All final stability decisions use exact arithmetic.

Compact circular input automatically uses exact cyclic symmetry reduction. Output retains one dihedral representative and a
`multiplier`; weighted candidate and ESS counts include that multiplier.

## Computation Flow

```text
CLI or Pybind input
  -> validating matrix parser
  -> support generation and exact-equilibrium superset pruning
  -> fast/test candidate attempt when selected
  -> exact candidate solver when selected or required as fallback
  -> exact stability check
  -> candidate and ESS output
```

Core implementation facts:

- Supports through dimension 64 are raw `uint64_t` masks; larger supports use fixed-width word vectors selected once before search.
  Production generates one support at a time and never materializes the complete support frontier.
- Non-circular supports use fixed-cardinality binary DFS. Circular supports use direct fixed-density bracelet generation.
- Every exact equilibrium support forbids later strict supersets; ESS status is irrelevant to that pruning rule.
- Fast/test convert and equilibrate the complete matrix once, solve reduced symmetric systems in binary64, and fall back to exact
  checking when preparation or a support solve is inconclusive.
- Safe clears denominators once, solves the border-eliminated symmetric candidate system with fraction-free exact integer arithmetic,
  and constructs rational probabilities and payoff only for successful public output.
- Stability reuses exact reduced-Hessian inertia. The rare case with unresolved outside best replies constructs an integer-scaled
  reduced Bomze matrix through an exact Schur complement and decides strict copositivity exactly.
- The final strict-copositivity path applies low-dimensional and sign checks, splits the negative-entry graph into connected
  components, and sends unresolved components to exact Hadeler enumeration.
- Hadeler checks principal submatrices by increasing cardinality. Each unresolved submatrix uses an exact fraction-free determinant,
  one retained solve when nonsingular, or one exact nullspace when singular. A binary64 value is never accepted as a stability
  certificate.
- Rational matrices remain at parsed-input and public-output boundaries. Exact candidate and stability kernels use
  `linalg::integer` and `linalg::matrix_int`.

## Build And Dependencies

Required dependencies are a C++17 compiler, CMake 3.18 or newer, Python with development headers when building Pybind, GMP, MPFR,
and FLINT. CMake FetchContent obtains `spdlog`, `argparse`, `pybind11`, and, when `BUILD_TESTING=ON`, GoogleTest. A clean configure
therefore needs network access unless those sources are cached.

Local non-MSVC Release builds default to `FRACESSA_NATIVE_ARCH=ON`; release CI disables native-CPU code generation. IPO/LTO is used
only when CMake confirms support. The optimized local FLINT 3.6 installation is `.local/flint-3.6.0`.

`cmake --install cpp/build` installs only the CLI. It does not install GMP, MPFR, or FLINT.

## Verification And Test Data

The maintained automated checks are the C++/CLI CTest suite and the Python unittest suite. Matrix correctness is not currently
registered as one CTest per database row.

`testdata/fracessa_testdata.sqlite3` is the canonical game, result, and timing store. Its schema, provenance, inventory, calibration
rules, and maintenance commands are in `testdata/README.md`.

Canonical performance comparisons use CPU 2, one persistent Pybind process, stored per-matrix calibrations, a 0.5-second default
target, and median native nanoseconds. CLI process timing and persistent-Pybind timing are not comparable.

## Python API

`python/pyfracessa/` is the maintained API. Public execution uses `run()` or `run_multiprocessing()` with an explicit `fast`, `safe`,
or `test` method. `compute_matrix()` is the public low-level native adapter. One matrix is always computed by one worker process;
multiprocessing is across matrices.

No production CLI, native, or Python workflow imposes a per-matrix computation timeout. A valid matrix may run for hours.

CLI, JSON, and CSV expose support masks as arbitrary-width decimal integers. Python returns ordinary arbitrary-precision `int`
objects. The current Parquet candidate schema remains `uint64` and rejects a candidate support above that range explicitly.

The public API, result schema, multiprocessing details, timing contract, examples, and generated-documentation commands are in
`aidocs/pyfracessa/README.md`.

## Public Documentation

GitHub Pages serves the landing page from `site/` and builds the combined API documentation into `site/docs/`. Sphinx reads Python
docstrings directly and Breathe renders Doxygen XML generated from `cpp/include/`; generated HTML and XML are not committed.

## Release

Releases are manual GitHub Actions runs from `main`. Calendar versioning is defined in `cpp/CMakeLists.txt`; successful workflow runs
create the tag and GitHub release and publish CLI binaries, wheels, and the source distribution. Failed builds create no release.

The complete release procedure is in `aidocs/RELEASING.md`.
