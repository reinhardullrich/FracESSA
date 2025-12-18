# FracESSA

**FracESSA - Fractional ESS Analyzer**

A high-performance C++ solver for finding Evolutionary Stable Strategies (ESS) in Standard Quadratic Problems using exact rational arithmetic.

## Overview

FracESSA solves the **Standard Quadratic Problem (SQP)**, which involves finding all **Evolutionary Stable Strategies (ESS)** for a given symmetric payoff matrix in evolutionary game theory. 

In the SQP, we seek strategies (probability distributions over pure strategies) that are:
- **Nash equilibria**: Best responses to themselves
- **Evolutionarily stable**: Resistant to invasion by alternative strategies

The problem is formulated as finding all local maximizers of a quadratic form over the standard simplex, where the payoff matrix defines the interactions between strategies. FracESSA uses exact rational arithmetic to guarantee correctness, avoiding floating-point errors that can occur with ill-conditioned matrices.

FracESSA computes all ESS using:
- **Exact rational arithmetic** via FLINT (Fast Library for Number Theory)
- **Optimized algorithms** including Hadeler's copositivity criterion and Bomze's stability checks
- **High-performance implementation** with aggressive optimizations (LTO, native architecture, static linking)

## Features

- ✅ Exact rational arithmetic (no floating-point errors)
- ✅ Support for both general symmetric and circular symmetric matrices
- ✅ Fast double-precision filtering before exact computation
- ✅ Python bindings for easy integration
- ✅ Highly optimized C++ implementation
- ✅ Cross-platform (Linux, macOS, Windows)

## Pre-built Binaries

Pre-built executables are available for all platforms in the [Releases](https://github.com/reinhardullrich/FracESSA/releases) section. Simply download the binary for your operating system:
- **Linux**: `fracessa-linux`
- **macOS**: `fracessa-macos`
- **Windows**: `fracessa-windows.exe`

Make the binary executable (Linux/macOS: `chmod +x fracessa-linux`) and run it directly—no compilation needed!

## Building

### Prerequisites

**Linux/macOS:**
- CMake 3.14+
- C++17 compiler (GCC/Clang)
- GMP and MPFR development libraries
- FLINT is built automatically from source

```bash
# Ubuntu/Debian
sudo apt-get install libgmp-dev libmpfr-dev cmake build-essential

# macOS
brew install gmp mpfr cmake
```

**Windows:**
- CMake 3.14+
- Visual Studio 2019+ with C++ support
- vcpkg for dependencies

```powershell
vcpkg install flint:x64-windows-static
```

### Build Instructions

```bash
cd cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
```

The executable will be at `build/fracessa` (or `build/Release/fracessa.exe` on Windows).

## Usage

### Command Line

FracESSA accepts matrices in the format `dimension#values`, where values are comma-separated rational numbers.

**Basic usage:**
```bash
./fracessa "2#0,1,0"
```

**With options:**
```bash
./fracessa -c -t "3#4,13/2,1/2,5,11/2,3"
```

**Options:**
- `-c, --candidates`: Include detailed candidate information in output
- `-t, --timing`: Output computation time
- `-e, --exact`: Use only exact arithmetic (slower, but handles extreme values)
- `-f, --fullsupport`: Search full support directly after size 1
- `-l, --log`: Enable detailed logging to `fracessa.log`
- `-m, --matrixid ID`: Optional matrix ID for logging
- `-u, --unsafe`: Use unsafe parsing (no validation, maximum performance)

**Matrix Formats:**
- **General symmetric**: Upper triangular format with `n*(n+1)/2` elements
  - Example: `"3#4,13/2,1/2,5,11/2,3"` for a 3×3 matrix
- **Circular symmetric**: Compact format with `floor(n/2)` elements
  - Example: `"5#1,3"` for a 5×5 circular symmetric matrix

### Python API

```python
from fracessa_py import Fracessa, Matrix

# Initialize
fracessa = Fracessa()

# Create a matrix
matrix = Matrix("0,1,0", dimension=2, is_circular=False)

# Compute ESS
result = fracessa.compute_ess(
    matrix=matrix,
    include_candidates=True,
    exact_arithmetic=False
)

print(f"Found {result.ess_count} ESS")
for candidate in result.candidates:
    if candidate.is_ess:
        print(f"ESS {candidate.candidate_id}: payoff={candidate.payoff_double:.6f}")
```

## Project Structure

```
fracessa/
├── cpp/
│   ├── include/
│   │   ├── fracessa/          # Main application headers
│   │   └── rational_linalg/   # Rational linear algebra library
│   └── src/                    # Implementation files
├── python/
│   ├── fracessa_py.py         # Python bindings
│   ├── run_matrices.py        # Batch processing script
│   └── verification/          # Test matrices and baselines
└── .github/workflows/         # CI/CD for releases
```

## Algorithm

FracESSA implements:
1. **Support enumeration** with efficient pruning
2. **Double-precision filtering** for fast rejection
3. **Exact rational computation** for candidates
4. **Stability checking** using:
   - Positive definiteness tests
   - Partial copositivity (Bomze 1992)
   - Full copositivity via Hadeler's criterion with memoization

## Performance

- Optimized for speed with aggressive compiler flags
- Static linking for maximum portability
- Efficient matrix operations and caching
- Thread-safe (designed for external multithreading)

## License

See LICENSE file for details.

## Citation

If you use FracESSA in your research, please cite appropriately.
