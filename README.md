<p align="center">
  <img src="logo.jpeg" width="900" alt="FracESSA logo" />
</p>

# FracESSA
Fast ESS analysis for symmetric games with exact rational verification.

FracESSA is a C++17 command-line tool for Evolutionary Stable Strategy (ESS) search on symmetric payoff matrices.
It is designed for raw speed in large support-space scans, with an exact-only
mode when correctness must not depend on floating-point filtering.

## Why FracESSA
- Default two-stage pipeline: an uncertified floating-point rejection filter,
  then exact FLINT fraction checks for surviving or suspicious supports.
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
./cpp/build/fracessa "3#4,13/2,1/2,5,11/2,3"
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
- `-c, --candidates` include candidate rows in output.
- `-l, --log` write detailed log output.
- `-e, --exact` disable floating-point candidate rejection and use exact solving.
- `-u, --unsafe` explicitly select the default uncertified numerical filter.
- `-f, --fullsupport` evaluate full support first.
- `-t, --timing` print analyzer timing in nanoseconds.
- `-m, --matrixid` optional signed 64-bit matrix ID for logging/verification runs.

`--exact` and `--unsafe` are mutually exclusive. With no numerical flag,
FracESSA temporarily uses unsafe numerical filtering. It can miss exact
candidates and ESS results; suspicious or unusable floating cases fall back to
exact arithmetic and may therefore run as slowly as exact mode.

Output format:
- line 1: ESS count
- optional line 2 (`-t`): runtime in nanoseconds
- optional extra lines (`-c`): candidate CSV rows

## Build Dependencies
- C++17 compiler
- CMake >= 3.18
- GMP, MPFR, FLINT
- Python 3.14 or newer with development headers (for the pybind module)

Notes:
- Third-party C++ dependencies (`spdlog`, `argparse`, `googletest`, `pybind11`) are pulled by CMake `FetchContent`.
- First configure in a fresh `cpp/build/` directory requires internet access unless dependencies are cached.
- Current Linux and macOS release binaries dynamically link GMP/MPFR/FLINT; they are not universal standalone binaries.

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
