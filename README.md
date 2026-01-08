<p align="right">
  <img src="logo.png" width="150" alt="FracESSA Logo">
</p>

# FracESSA
### Fractional Evolutionary Stable Strategy Analyzer

FracESSA is a high-performance computational engine for identifying **Evolutionary Stable Strategies (ESS)** within the framework of **Standard Quadratic Problems (SQP)**. The software provides a robust implementation of stability analysis in evolutionary game theory, utilizing exact rational arithmetic to ensure mathematical precision.

By employing exact algebraic methods, FracESSA eliminates numerical instabilities inherent in floating-point computations, enabling reliable analysis of complex payoff matrices where stability boundaries are extremely sensitive to precision.

---

## Technical Overview

The identification of ESS represents a fundamental problem in evolutionary biology and game theory, requiring determination of local maximizers of quadratic forms over the standard simplex. FracESSA implements a combinatorial approach optimized for high-dimensional matrices, combining exact algebraic methods with modern computational techniques.

### Algorithmic Foundation
The core methodology systematically explores supports of the strategy simplex, verifying necessary and sufficient conditions for evolutionary stability through combinatorial enumeration and exact verification.

---

## Core Features and Implementation

### Exact Rational Computation
At the heart of FracESSA lies the **Fast Library for Number Theory (FLINT)**. Every critical computation—from determinant evaluation to linear system solving—is performed using exact fractions. This approach ensures mathematical rigor essential for research applications where stability region boundaries demand absolute precision.

### Performance Optimization
To maintain computational efficiency while preserving exactness, FracESSA employs a multi-stage filtering pipeline:
1. **Double-Precision Screening**: Potential supports are initially evaluated using standard double-precision arithmetic, rapidly eliminating non-viable candidates.
2. **Exact Algebraic Verification**: Surviving candidates undergo full verification using exact rational arithmetic to confirm both Nash equilibrium conditions and stability properties.

### Support for Complex Symmetries
While supporting general symmetric payoff matrices, the analyzer includes optimized routines for **Circular Symmetric Matrices**. This feature significantly reduces computational complexity for large-scale periodic models prevalent in spatial ecology and evolutionary dynamics.

---

## Usage and Integration

### Command Line Interface
The primary interface for FracESSA is a standalone CLI tool. Matrices are specified in a compact triangular format: `dimension#values`.

**Standard Analysis:**
```bash
# Analyze a 3-dimensional symmetric matrix
./fracessa "3#4,13/2,1/2,5,11/2,3"
```

**Advanced Parameters:**
The application supports several configuration flags:
- `--candidates (-c)`: Outputs detailed reports of examined supports, including payoffs and stability status.
- `--timing (-t)`: Provides high-resolution execution benchmarks.
- `--exact (-e)`: Forces exclusive use of exact arithmetic, bypassing double-precision filtering.
- `--log (-l)`: Generates comprehensive execution profiles in `fracessa.log`.

### Python API
For integration into research workflows, a native Python interface is provided via `fracessa_py`.

```python
from fracessa_py import Fracessa, Matrix

# Initialize the analyzer
analyzer = Fracessa()

# Define payoff matrix structure
matrix = Matrix("0,1,0", dimension=2)

# Compute complete set of ESS
result = analyzer.compute_ess(matrix, include_candidates=True)

print(f"Total ESS Identified: {result.ess_count}")
for candidate in result.candidates:
    if candidate.is_ess:
        print(f"Refined Support: {candidate.support_bits}")
```

---

## Building and Environment

### System Dependencies
- **Linux/macOS**: Requires development headers for `GMP` and `MPFR`. `FLINT` is integrated into the build process.
- **Windows**: Distribution via `vcpkg` is supported using the `flint:x64-windows-static-release` triplet.

### Compilation
The project uses a standard CMake build system:
```bash
cd cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
```

---

<p align="center">High-Performance Computational Tool for Evolutionary Game Theory • Developed by Reinhard Ullrich</p>
