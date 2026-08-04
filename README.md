<p align="center">
  <img src="logo.png" width="600" alt="FracESSA logo" />
</p>

# FracESSA

FracESSA finds evolutionarily stable strategies (ESS) for symmetric rational payoff matrices. It provides an exact search for
reliable results and a faster floating-point-assisted search for exploratory work.

## The problem

Let $A=A^\mathsf{T}$ be a symmetric payoff matrix and let $x$ be a probability vector. FracESSA studies the quadratic form

$$
\max_{x\in\Delta_n} x^\mathsf{T} A x,
\qquad
\Delta_n=\{x\in\mathbb{R}^n \mid x_i\geq0,\ \sum_{i=1}^n x_i=1\}.
$$

For a symmetric matrix, the ESS are exactly the strict local maximizers of this problem. FracESSA searches the possible supports,
checks the corresponding equilibrium candidates, and classifies their stability. It returns the number of ESS and can also return
their exact probability vectors, supports, payoffs, and stability results.

The support space has size $2^n$, so running time can grow quickly with the dimension and structure of the matrix. FracESSA is
optimized for repeated operations on many small support systems.

## Capabilities and limitations

FracESSA accepts symmetric payoff matrices whose entries are integers or exact fractions. It can find every equilibrium candidate
and every ESS rather than only one solution. Every reported probability vector and payoff is an exact rational value, regardless of
whether `safe` or `fast` was selected. The search method affects completeness, not the exactness of the reported output. Numerators
and denominators are not restricted to the precision or exponent range of binary64 floating-point numbers.

The main limitations are:

- the matrix must be symmetric and contain rational values;
- the dimension must be smaller than 64;
- the $2^n$ support space makes some larger or difficult matrices computationally expensive;
- one matrix is processed on one CPU core, although the Python API can process several matrices in parallel.

## Choose a search method

FracESSA always requires an explicit method; there is no default.

| Method | Use it when | Guarantee |
|---|---|---|
| `safe` | The result must be complete and mathematically reliable | Uses exact arithmetic for every candidate decision |
| `fast` | Speed matters more than a completeness guarantee | Uses floating-point filtering before exact checks and can miss candidates |

Candidates returned by `fast` are checked exactly, but its early floating-point filtering can make the returned set incomplete.
When `fast` detects certain dangerous whole-matrix conditions, it automatically falls back to `safe` and reports the reason.

The CLI and Python API also expose `test`, an experimental copy of `fast` used for development. It is not intended for normal use.

## Build from source

Required system dependencies:

- a C++17 compiler;
- CMake 3.18 or newer;
- Python 3.14 or newer, including development headers;
- GMP, MPFR, and FLINT;
- Git and internet access for CMake's first download of `argparse`, `pybind11`, `spdlog`, and the test framework.

Build the command-line program and Python extension from the repository root:

```bash
git clone https://github.com/reinhardullrich/FracESSA.git
cd FracESSA
./build.sh
```

The resulting programs are:

- `cpp/build/fracessa`: command-line interface;
- `cpp/build/fracessa_core.*`: native module used by Python.

Run the automated test suite with:

```bash
./test.sh
```

## Command-line use

The general form is:

```bash
./cpp/build/fracessa [OPTIONS] METHOD "DIMENSION#VALUES"
```

For example, the symmetric matrix

$$
A=\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$

is written as `2#0,1,0`. Its exact ESS is $(1/2,1/2)$:

```bash
./cpp/build/fracessa safe --candidates "2#0,1,0"
```

```text
1
candidate_id;vector;support;support_size;extended_support;extended_support_size;multiplier;is_ess;stability;payoff;payoff_dbl
1;1/2,1/2;3;2;3;2;;1;T_pd_frc;1/2;0.5
```

Without `--candidates`, FracESSA prints only the ESS count.

In the candidate table, `vector` and `payoff` retain exact fractions, `support` is a bit mask, and `is_ess` is `1` for an ESS.
For compact circular input, `multiplier` says how many rotations and reflections the displayed representative covers.

Useful options:

| Option | Meaning |
|---|---|
| `-c`, `--candidates` | Print the candidate table after the ESS count |
| `-t`, `--timing` | Print native runtime in nanoseconds and the whole-matrix safe-fallback reason |
| `-f`, `--fullsupport` | Check the full support first |
| `-l`, `--log` | Write detailed diagnostics to `log/fracessa.log` |
| `-m`, `--matrixid ID` | Attach a signed 64-bit matrix ID to the log |

Run `./cpp/build/fracessa --help` for the current command-line reference.

## Matrix input

Every input begins with the dimension, followed by `#` and a comma-separated value list:

```text
dimension#values
```

Values must be integers or exact integer fractions such as `-3/5`. Decimal notation is not accepted; write `1/2` instead of
`0.5`.

### General symmetric matrices

Provide the upper triangle row by row. For

$$
A=\begin{pmatrix}
a_{11}&a_{12}&a_{13}\\
a_{12}&a_{22}&a_{23}\\
a_{13}&a_{23}&a_{33}
\end{pmatrix},
$$

the input is:

```text
3#a11,a12,a13,a22,a23,a33
```

Exactly $n(n+1)/2$ values are required.

### Circular-symmetric matrices

A circular-symmetric matrix with zero diagonal can use a compact list of $\lfloor n/2\rfloor$ values. Entry $c_k$ is the payoff
at circular distance $k$:

```text
n#c1,c2,...,c_floor(n/2)
```

For example, `5#1,3` expands to a 5-by-5 matrix with first row `0,1,3,3,1`. FracESSA recognizes this encoding from the shorter
value count and automatically applies its circular-symmetry reductions.

## Python use

The maintained Python API calls the same C++ engine in-process. Build the project first, then make the package visible from the
repository root:

```bash
PYTHONPATH=python python3 your_script.py
```

A complete single-matrix example is:

```python
from fracessa import Matrix, RunConfig, run

matrix = Matrix(1, "2#0,1,0")
result = run(
    "safe",
    matrix,
    config=RunConfig(include_candidates=True),
    run_id="example",
)

print(result["ess_count"])
print(result["candidates"])
```

The result is a plain dictionary containing the status, candidate and ESS counts, support-size structures, native runtime,
optional safe-fallback reason, optional candidate rows, and the input metadata.

For many independent matrices, use `run_multiprocessing()` to distribute matrices across worker processes. One matrix remains
single-core; parallelism is across matrices. The [Python API guide](aidocs/pyfracessa/README.md) documents sequential execution,
multiprocessing, JSON input, CSV/JSON/Parquet output, and every result field.

## Further documentation

- [Python API](aidocs/pyfracessa/README.md)
- [Test-matrix database](testdata/README.md)
- [Technical documentation index](aidocs/INDEX.md)
- [Research papers](research/papers/)

The root README is intentionally an introduction and usage guide. Mathematical derivations, numerical details, benchmarks, and
historical experiments are kept in the linked documentation instead.
