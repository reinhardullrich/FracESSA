<p align="center">
  <img src="logo.jpeg" width="900" alt="FracESSA logo" />
</p>

# FracESSA
Fast ESS analysis for symmetric games with exact rational verification.

FracESSA is a C++17 command-line tool for Evolutionary Stable Strategy (ESS) search on symmetric payoff matrices.
It is designed for raw speed in large support-space scans, with an exact-only
safe method, a faster heuristic method, and an isolated experimental copy of that heuristic.

## Why FracESSA
- Required method choice: `fast` uses binary64 rejection before exact checks;
  `safe` uses exact FLINT arithmetic from the start.
- Experimental `test` independently copies `fast` so proposed numerical changes can be measured without changing production fast search.
- Bitset-based support enumeration over a `2^n` search space.
- Optimized for many repeated operations on small/medium matrix dimensions.
- Circular-symmetric and general symmetric matrix input support.

## Quick Start

### 1. Build
```bash
./build.sh
```

### 2. Run Tests
```bash
./test.sh
```

### 3. Run One Matrix
```bash
./cpp/build/fracessa safe "3#4,13/2,1/2,5,11/2,3"
```

## Input Format
Matrix input is:
```text
dimension#values
```

Supported encodings:
- Symmetric upper-triangular: `n*(n+1)/2` values.
- Circular-symmetric compact: `floor(n/2)` values.

The parser accepts dimensions 1 through 63 and validates the complete matrix
syntax. The 64-bit mask is storage, not support for a dimension-64 search.

## CLI Flags

The first positional argument is required and must be `fast`, `safe`, or experimental `test`; there
is no default method.

- `-c, --candidates` include candidate rows in output.
- `-l, --log` write detailed log output.
- `-f, --fullsupport` evaluate full support first.
- `-t, --timing` print analyzer timing in nanoseconds and the whole-matrix safe fallback.
- `-m, --matrixid` optional signed 64-bit matrix ID for logging/verification runs.

For every circular-symmetric matrix, the cyclic symmetry filter checks exact index multipliers that preserve the complete matrix,
then solves only one bracelet from each detected affine orbit. If it finds no extra multiplier, it disables itself after this
one-time check. Candidate output remains universal: every stored row represents only its rotations and reflections, so its
`multiplier` is at most twice the dimension. Affine-equivalent dihedral rows are reconstructed by exact permutation without
solving their systems again. The filter does not run on non-circular input or change the selected candidate method.

`fast` removes the game's common denominator and switches the whole matrix to safe search when the remaining exact integer
precision span satisfies $P\geq10^9$. Otherwise it normalizes and equilibrates the complete binary64 game once, eliminates the
candidate border, and solves each reduced symmetric system with Bunch-Kaufman $LDL^T$. An inconclusive pivot below $10^{-12}$
sends that support to exact checking. The remaining probability and outside-payoff rejections are heuristic, so `fast` is not a
mathematical correctness certificate. `safe` bypasses floating-point candidate rejection and uses exact arithmetic for every
support.

`test` is an independent source copy of `fast` for numerical experiments. It currently has the same behavior as `fast`.

Output format:
- line 1: ESS count
- optional line 2 (`-t`): runtime in nanoseconds
- optional line 3 (`-t`): `null`, `precision_span`, `equilibration_invalid`, or `equilibration_non_convergence`
- optional extra lines (`-c`): candidate CSV rows

Candidate CSV and log rows print `payoff_dbl` with enough significant digits to reconstruct the original binary64 value.

## Build Dependencies
- C++17 compiler
- CMake >= 3.18
- GMP, MPFR, FLINT
- Python 3.14 or newer with development headers (for the pybind module)

Notes:
- Third-party C++ dependencies (`spdlog`, `argparse`, `googletest`, `pybind11`) are pulled by CMake `FetchContent`.
- First configure in a fresh `cpp/build/` directory requires internet access unless dependencies are cached.
- Current Linux and macOS release binaries dynamically link GMP/MPFR/FLINT; they are not universal standalone binaries.

## Versioning

FracESSA uses calendar versions in the form `YEAR.MONTH.DAY.RELEASE_OF_DAY`, for example `2026.8.3.1`. The final number starts at
1 each day and increments for additional releases on that date. CMake owns the version compiled into the CLI; a release tag must
contain the same value with a `v` prefix, such as `v2026.8.3.1`, or the release workflow stops before publishing artifacts.

## Repository Layout
- `cpp/`: core engine, tests, and local build
- `python/`: maintained Python wrapper
- `testdata/`: canonical SQLite test matrices and expected candidates
- `aidocs/`: current project knowledge, open issues, wrapper docs, benchmark reports, and references
- `experiments/`: historical benchmark code, metadata, and results
- `.github/workflows/`: build/test CI and tagged-release automation

## Project Focus
FracESSA is explicitly tuned for high-throughput ESS workflows where millions of support-level operations are common.
Micro-overheads in hot paths matter and are optimized aggressively by design.
