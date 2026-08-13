<p align="center">
  <img src="https://reinhardullrich.github.io/fracessa/logo.png" width="600" alt="FracESSA logo" />
</p>

# FracESSA

## The problem

We consider the *Standard Quadratic Optimization Problem (StQP)* given by

$$
\max_{\mathbf x\in \Delta^n} \mathbf x^\top{\mathsf A}\mathbf x
$$

where ${\mathsf A}$ is a symmetric $n\times n$ matrix and $\Delta^n$ is the standard simplex

$$
\Delta^n=\lbrace\mathbf x\in\mathbb{R}^{n}:\sum_{i=1}^{n}x_i=1,\ x_i\geq0\quad\text{for }i=1,\ldots,n\rbrace.
$$

When ${\mathsf A}$ is interpreted as the payoff matrix of a symmetric partnership game, $\mathbf x^\top{\mathsf A}\mathbf x$ is the
population's average payoff under the mixed strategy $\mathbf x$. Its ESS are exactly the strict local maximizers of this quadratic
form. FracESSA reports the number and support-size structure of the candidates and ESS, and can also return every exact probability
vector, support, payoff, and stability result. It is available as a self-contained command-line program and as the `pyfracessa`
Python package.

## Capabilities and limitations

FracESSA accepts integers, fractions, finite decimals, and scientific notation as exact rational values. For example, `0.125`, `1/8`,
and `1.25e-1` are identical inputs. Numerators and denominators use arbitrary-precision integers, so their size is limited by memory
rather than a fixed machine type. The reported `vector` and `payoff` fields are exact; `payoff_dbl` is only a convenient floating-point
approximation.

The main limitations are:

- the matrix must be symmetric and contain rational values;
- the search may have to inspect up to $2^n-1$ nonempty supports, so difficult matrices can take a long time;
- one matrix is processed on one CPU core, although the Python API can process several matrices in parallel.

There is no fixed dimension cap, but exponential running time is normally the practical limit.

## Choose a search method

FracESSA always requires an explicit method; there is no default.

| Method | Use it when | Guarantee |
|---|---|---|
| `safe` | The result must be complete and mathematically reliable | Uses exact arithmetic for every candidate decision |
| `fast` | Speed matters more than a completeness guarantee | Uses floating-point filtering before exact checks and can miss candidates |

Candidates returned by `fast` are checked exactly, but its early floating-point filtering can make the returned set incomplete.
When `fast` detects certain dangerous whole-matrix conditions, it automatically falls back to `safe`. Python returns the reason in
`safe_fallback`; the CLI includes the same field in every JSON summary.

## Install

The GitHub [Releases](https://github.com/reinhardullrich/fracessa/releases) page provides self-contained command-line binaries for
Linux x86-64, Linux ARM64, macOS Intel, macOS Apple Silicon, Windows x86-64, and Windows ARM64. Python, GMP, MPFR, and FLINT do not
need to be installed separately. On Linux and macOS, make the downloaded file executable before running it:

```bash
chmod +x fracessa-*
```

Windows users can run the downloaded `.exe` directly. The Linux binaries require glibc 2.28 or newer, and the macOS binaries
require macOS 11 or newer.

Python 3.11 through 3.14 users can install the native extension and Python API directly from PyPI:

```bash
python -m pip install pyfracessa
```

The distribution and import package are both named `pyfracessa`. Parquet output is optional:

```bash
python -m pip install "pyfracessa[parquet]"
```

## Command-line use

Use the path to the downloaded binary, or `./cpp/build/fracessa` after building from source. The general form for inline input is:

```bash
./cpp/build/fracessa [OPTIONS] METHOD "DIMENSION#VALUES"
```

The final argument can instead be the path to a file containing the same compact syntax or a Matrix Market matrix.

For example, the symmetric matrix

$$
A=\begin{pmatrix}
-1 & 0 & 0 \\
0 & -2 & 0 \\
0 & 0 & -3
\end{pmatrix}
$$

is written as `3#-1,0,0,-2,0,-3`. Its exact ESS is

$$
\mathbf{x}=\left(\frac{6}{11},\frac{3}{11},\frac{2}{11}\right).
$$

```bash
./cpp/build/fracessa safe --candidates "3#-1,0,0,-2,0,-3"
```

```text
{"run_id":null,"matrix_id":-1,"status":0,"candidate_count":1,"ess_count":1,"candidate_structure":{"3":1},"ess_structure":{"3":1},"elapsed_ns":33291,"safe_fallback":null,"error_message":""}
candidate_id;vector;support;support_size;extended_support;extended_support_size;multiplier;is_ess;stability;payoff;payoff_dbl
1;6/11,3/11,2/11;7;3;7;3;;1;T_reduced_hessian_nd;-6/11;-0.54545454545454541
```

Every run writes one JSON summary containing the status, candidate and ESS counts, their support-size structures, runtime, any
automatic safe fallback, and any error. `--candidates` appends the candidate CSV shown above. Its `vector` and `payoff` columns are
exact; `payoff_dbl` is only an approximation. For circular input, `multiplier` records how many rotations and reflections the
representative covers. The displayed runtime is only an example and varies between runs.

Useful options:

| Option | Meaning |
|---|---|
| `-c`, `--candidates` | Print the candidate CSV after the JSON summary |
| `-f`, `--fullsupport` | Check the full support first |
| `-l`, `--log` | Write detailed diagnostics to `log/fracessa.log` |
| `-m`, `--matrixid ID` | Set the signed 64-bit matrix ID in the JSON summary and log |

Run `./cpp/build/fracessa --help` for the current command-line reference.

## Matrix input

Inline input begins with the dimension, followed by `#` and a comma-separated value list:

```text
dimension#values
```

Values can be integers, fractions, finite decimals, or scientific numbers, such as `-3/5`, `-0.6`, or `-6e-1`. All forms are parsed
as exact rational numbers. A fraction's denominator must be a positive integer; put an optional sign before the numerator, writing
`-1/2` rather than `1/-2`.

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

Exactly $\frac{n(n+1)}{2}$ values are required.

### Circular-symmetric matrices

For dimensions 2 and larger, a circular-symmetric matrix with zero diagonal can use a compact list of
$\left\lfloor \frac{n}{2} \right\rfloor$ values. Entry $c_k$ is the payoff at circular distance $k$:

```text
n#c1,c2,...,c_floor(n/2)
```

For example, `5#1,3` expands to a 5-by-5 matrix with first row `0,1,3,3,1`. FracESSA recognizes this encoding from the shorter
value count and automatically applies its circular-symmetry reductions.

### Matrix files

The CLI also accepts a file path in place of inline input. The file can contain either `dimension#values` or a symmetric Matrix
Market matrix. Python accepts the same two text formats through `Matrix.matrix`; read a file with `Path.read_text()` when needed.

## Python use

The maintained Python API calls the same C++ engine in-process. Install `pyfracessa` from PyPI as shown above, or build the project
from source and make the source package visible from the repository root:

```bash
PYTHONPATH=python python3 your_script.py
```

A complete single-matrix example is:

```python
from pyfracessa import Matrix, RunConfig, run

matrix = Matrix(1, "3#-1,0,0,-2,0,-3")
result = run(
    "safe",
    matrix,
    config=RunConfig(include_candidates=True),
    run_id="example",
)

print(result["candidate_count"])
print(result["candidate_structure"])
print(result["ess_count"])
print(result["ess_structure"])
print(result["candidates"])
```

The result is a plain dictionary containing the same summary information as the CLI JSON line. Counts and support-size structures are
always present. `RunConfig(include_candidates=True)` also returns the individual candidate rows with exact vectors and payoffs.

### Process several matrices in parallel

`run_multiprocessing()` distributes independent matrices across worker processes. One matrix still uses one CPU core; the
parallelism is across matrices. Results are yielded in completion order, so use `matrix_id` to identify each one.

```python
from pyfracessa import MPConfig, Matrix, RunConfig, run_multiprocessing

matrices = [
    Matrix(1, "3#-1,0,0,-2,0,-3"),
    Matrix(2, "3#4,13/2,1/2,5,11/2,3"),
]

if __name__ == "__main__":
    mp_config = MPConfig(workers=8)

    for result in run_multiprocessing(
        "safe",
        matrices,
        config=RunConfig(include_candidates=False),
        mp_config=mp_config,
    ):
        print(result["matrix_id"], result["ess_count"])
```

`MPConfig()` without arguments uses all CPUs available to the Python process. Set `workers` to choose a smaller number. The
`if __name__ == "__main__":` guard is required by the default portable process-start method.

The [Python API guide](https://reinhardullrich.github.io/fracessa/python-api.html) documents sequential execution, JSON input,
CSV/JSON/Parquet output, sinks, and every result field.

## Build from source

Required system dependencies:

- a C compiler and a C++17 compiler;
- CMake 3.18 or newer;
- Python 3.11 or newer, including development headers;
- GMP, MPFR, and FLINT;
- Git and internet access for CMake's first download of `argparse`, `pybind11`, `spdlog`, and the test framework.

Build the command-line program and Python extension from the repository root:

```bash
git clone --recurse-submodules https://github.com/reinhardullrich/fracessa.git
cd fracessa
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build --parallel
```

The resulting programs are:

- `cpp/build/fracessa`: command-line interface;
- `cpp/build/fracessa_core.*`: native module used by Python.

Run the automated test suite with:

```bash
ctest --test-dir cpp/build --output-on-failure --parallel
PYTHONPATH=python python3 -m unittest discover -s python/tests -p 'test_*.py'
```

## Further documentation

- [Public documentation](https://reinhardullrich.github.io/fracessa/)
- [Python API](https://reinhardullrich.github.io/fracessa/python-api.html)
- [Test-matrix database](testdata/README.md)
- [Release procedure](aidocs/RELEASING.md)
- [Technical documentation index](aidocs/INDEX.md)

## License

FracESSA is free software licensed under [GPL-3.0-or-later](LICENSE). Prebuilt releases statically link third-party libraries; their
licenses and exact source locations are recorded in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
