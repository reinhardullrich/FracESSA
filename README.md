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

## Algorithmic workflow

The potential search space consists of the $2^n-1$ nonempty supports, processed by increasing cardinality without
storing the full frontier. For a general matrix, a bitwise depth-first generator streams the supports; for a
circular-symmetric matrix, a fixed-density bracelet generator keeps one representative from each rotation/reflection
class, and an exact affine-symmetry filter removes further equivalent cases. Following
[Bomze's support criterion](research/papers/bomze_1992.md), every exact equilibrium support rules out all of its
strict supersets, so those branches are pruned before their systems are built. For each remaining support of size $k$,
FracESSA uses the simplex normalization to eliminate one probability and removes the equilibrium payoff from the
bordered system, producing a symmetric $(k-1)\times(k-1)$ reduced Hessian. In `safe`, one exact fraction-free
$LDL^\mathsf{T}$-style factorization then does two jobs: it solves for the candidate and supplies the Hessian inertia
needed for stability. If the inertia shows that the Hessian is not negative definite, the candidate is not an ESS. If it
is negative definite and no unused outside strategy ties the candidate's payoff, the candidate is an ESS; only the
remaining cases continue to Bomze's exact copositivity test. The `fast` method can reject some supports earlier with
floating-point calculations, but every candidate it reports and every final stability decision are exact.

Circular-symmetric matrices are symmetric matrices that remain unchanged when all strategy labels are shifted together
around a cycle. This creates many equivalent supports related by rotations and reflections. Circular symmetry does not
guarantee many ESS, but this matrix class has produced important games with unusually many strict local maxima and hence
many ESS. FracESSA therefore provides specialized support generation and symmetry reductions for these matrices.

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

## Install

The GitHub [Releases](https://github.com/reinhardullrich/fracessa/releases) page provides self-contained command-line archives for
Linux x86-64, Linux ARM64, macOS Intel, macOS Apple Silicon, and Windows x86-64. Extract the archive and run `fracessa`
(`fracessa.exe` on Windows); Python, GMP, MPFR, and FLINT do not need to be installed separately.

Python 3.11 through 3.14 users can install the native extension and Python API directly from PyPI:

```bash
python -m pip install pyfracessa
```

The distribution is named `pyfracessa`; the import remains `fracessa`. Parquet output is optional:

```bash
python -m pip install "pyfracessa[parquet]"
```

## Build from source

Required system dependencies:

- a C++17 compiler;
- CMake 3.18 or newer;
- Python 3.11 or newer, including development headers;
- GMP, MPFR, and FLINT;
- Git and internet access for CMake's first download of `argparse`, `pybind11`, `spdlog`, and the test framework.

Build the command-line program and Python extension from the repository root:

```bash
git clone https://github.com/reinhardullrich/fracessa.git
cd fracessa
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

The maintained Python API calls the same C++ engine in-process. Install `pyfracessa` from PyPI as shown above, or build the project
from source and make the source package visible from the repository root:

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
- [Release procedure](aidocs/RELEASING.md)
- [Technical documentation index](aidocs/INDEX.md)
- [Research papers](research/papers/)

## License

FracESSA is free software licensed under [GPL-3.0-or-later](LICENSE). Prebuilt releases statically link third-party libraries; their
licenses and exact source locations are recorded in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The root README is intentionally an introduction and usage guide. Mathematical derivations, numerical details, benchmarks, and
historical experiments are kept in the linked documentation instead.
