# coposit Integration And Exact-Type Ownership

Status: implemented and verified locally. The coposit foundation, FracESSA representation migration, strict-copositivity backend,
candidate-vector and scalar-output cleanup, obsolete exact-code deletion, and release packaging are complete.

## Decision

The general direction is correct, with two important qualifications:

1. The common input denominator must belong to the parsed input as a whole. It must not be hidden inside `integer` or
   `matrix_integer`.
2. The representation cleanup and the replacement of FracESSA's strict-copositivity algorithm must be separate measured changes.
   Otherwise, a timing or behavior change could not be attributed reliably.

The intended dependency direction is:

```text
FracESSA -> coposit
```

coposit must never depend on FracESSA. A third common repository is unnecessary.

## Current Checkpoint Evidence

FracESSA now uses coposit's parser and exact integer types, owns one parsed integer game plus its denominator, and keeps fast-mode
preparation outside the exact solver. Candidate probabilities use `std::vector<fraction>`; the duplicate parser and `matrix_frc`
storage have been deleted. The final strict-copositivity decision is one call to `coposit::safe`.

Verification on 2026-08-12 passed all 8 retained C++/CLI tests and all 70 Python tests after the obsolete backend tests were deleted.
The paired persistent-Pybind quick-suite representation benchmark on CPU 2 matched all 81 saved ESS results in both `fast` and `safe`.
Excluding dimension 2, the per-matrix median time
changed by `-0.929%` for `fast` and `-0.716%` for `safe`. The separately measured coposit backend changed the median by `-0.125%`
for `fast` and `+0.149%` for `safe`; among baseline runs of at least 1 ms, both mean and median changes stayed within `0.4%`.
Conditioning rational payoff and probability construction on requested candidate output then changed the median by `-0.888%` for
`fast` and `-0.719%` for `safe`. Among baseline runs of at least 1 ms, the corresponding medians were `-0.574%` and `-0.256%`.
Every benchmark row retained the expected ESS count; no stored baseline was rerun.
Slimming the remaining rational type to its output-value contract was neutral against that checkpoint: excluding dimension 2, the
median changes were `-0.002%` for `fast` and `-0.218%` for `safe`; among baseline runs of at least 1 ms, they were `-0.006%` and
`-0.196%`. All 162 matched rows retained their expected ESS count. coposit's complete 65-test suite passes independently at pinned
commit `901cdfc1bb550588977d8255a4483ebe72e3b979`.

Release jobs now check out the pinned submodule, the PyPI sdist contains only the coposit files required to build FracESSA, and the
GitHub release attaches a complete source archive rather than mislabeling GitHub's submodule-less automatic archive. The final
131,613-byte sdist was extracted without Git metadata, rebuilt into a CPython 3.14 ARM64 wheel, installed in a clean environment,
and passed a native safe-mode ESS smoke test with exact candidate output.
The separate complete source archive contains both repositories.

## Final Ownership

| Component | Decision | Owner |
|---|---|---|
| `matrix_dbl` | Keep unchanged | FracESSA |
| `matrix_frc` | Remove | None |
| `fraction` | Keep as a small output-only rational value | FracESSA |
| `integer` | Use coposit's implementation | coposit |
| `matrix_int` | Alias coposit's `matrix_integer` | coposit |
| Generic fraction-free LDLT | Use coposit's implementation | coposit |
| Candidate-specific fraction-free KKT LDLT | Keep | FracESSA |
| Exact matrix parsers | Use coposit's parsers | coposit |
| Strict-copositivity implementation | Call one public coposit function | coposit |

## Keep `matrix_dbl`

`matrix_dbl` is a small row-major dense matrix used by the fast candidate filter for:

- the complete converted double game;
- the reusable reduced candidate system;
- in-place Bunch-Kaufman factorization.

Replacing it with `std::vector<double>` would not eliminate meaningful code. It would merely distribute expressions such as
`row * dimension + column` throughout the numerical solver.

coposit does not need a corresponding dense double factorization matrix. A shared or templated matrix class would therefore add an
abstraction without removing duplication.

The fast-filter constructor now receives only the dimension and prepares its double matrix directly from the shared integer matrix.

## Remove `matrix_frc`

Before removal, `matrix_frc` combined four unrelated responsibilities:

| Former responsibility | Replacement |
|---|---|
| Parsed game storage | coposit `parsed_matrix` |
| Compact circular-matrix expansion | coposit parser |
| Matrix retained for logging and symmetry detection | Shared integer matrix plus denominator |
| Candidate probability vector | `std::vector<fraction>` |

The former input path was:

```text
text
  -> complete matrix_frc
  -> copy retained by FracESSA
  -> temporary FLINT rational matrix
  -> denominator clearing
  -> complete matrix_int
```

The implemented path is:

```text
text
  -> coposit parser
  -> integer matrix Z and common denominator d
  -> move into FracESSA
```

Here, `A = Z / d` and `d > 0`.

This removes:

- one complete rational matrix and its copy;
- the temporary `fmpq_mat` used only for denominator clearing;
- the later rational-to-integer conversion;
- `matrix_int::set_from_fraction_matrix()`;
- FracESSA's duplicate parser;
- `create_symmetric()` and `create_circular_symmetric()`.

### Logging

Logging does not require a rational matrix object. A small FracESSA formatting function can display each entry as `Z_ij / d`.
Pretty diagnostic output is a FracESSA responsibility and must not become a method of coposit's generic integer matrix.

### Circular symmetry

Exact symmetry detection can compare entries of `Z` directly because the entire matrix has one common positive denominator:

$$
A_{ij}=A_{kl}\quad\Longleftrightarrow\quad Z_{ij}=Z_{kl}.
$$

The affine-symmetry detector therefore receives the integer matrix. No rational matrix is required.

## Represent The Parsed Input Explicitly

coposit returns one small aggregate:

```cpp
struct parsed_matrix {
    matrix_integer matrix;       // Z, where A = Z / denominator
    integer denominator{1};      // Always positive
    bool compact_circular = false;
};
```

The denominator must not be stored in `integer`: one integer has no denominator.

It must also not be stored in `matrix_integer`. Most integer matrices in both projects are derived systems, factors, right-hand
sides, Schur complements, or workspaces. They do not all retain the original game's scale. Hidden denominator metadata would make
matrix operations ambiguous and could allow a derived matrix to retain a meaningless input denominator.

The denominator is metadata for the parsed input as a whole. Keeping it in `parsed_matrix` makes that invariant explicit.

coposit uses `parsed.matrix` and ignores the denominator because copositivity is invariant under multiplication by a positive
scalar. FracESSA uses:

- `matrix` for exact candidate and stability calculations;
- `matrix` for fast double conversion;
- `denominator` when constructing exact payoffs;
- `compact_circular` when selecting circular support generation.

Candidate probabilities do not need the game denominator. A common positive scaling of the payoff matrix does not change the
equilibrium vector.

## Store The Complete Integer Game Once

The analyzer now owns the parsed integer game once:

```text
fracessa::basic_analyzer
  owns parsed_matrix { Z, d, compact_circular }

  exact_candidate_solver
    holds const references to Z and d

  fast_candidate_filter
    reads Z once and owns its converted matrix_dbl
```

The parser result is moved into the analyzer. `matrix_integer` has an inexpensive move operation implemented by swapping FLINT
matrix storage.

The analyzer is already non-copyable and non-movable. References retained by its solvers are therefore valid for the analyzer's
complete lifetime.

## Keep And Reduce `fraction`

Exact probabilities and payoffs are public FracESSA results, so the scalar `fraction` class should remain. Replacing it with separate
numerator and denominator fields would recreate sign normalization, greatest-common-divisor reduction, canonical formatting, and
copy/move management already handled correctly by FLINT's `fmpq`.

After the parser migration, `fraction` should retain only operations that output actually uses:

- default construction and destruction;
- copy and move operations;
- `set_zero()`;
- `set_ratio()`;
- `to_string()`;
- `to_dbl()`;
- stream output if candidate CSV serialization continues to use it.

Likely removable operations include:

- parsing constructors;
- arithmetic operators;
- equality and `is_zero()` if no remaining caller uses them;
- static constants such as `zero()`, `one()`, `neg_one()`, and `two()`.

Candidate probabilities should become:

```cpp
std::vector<fracessa::numeric::fraction> vector;
```

That reflects their actual meaning as an `n`-element vector rather than an `n x 1` matrix.

coposit should not acquire this rational output class. Its maintained numerical core remains integer-only.

## Materialize Rational Output Only When Requested — Implemented

The exact solver previously constructed the exact rational payoff and `payoff_dbl` for every successful candidate, even when candidate
output was not requested. Only the dense probability vector was conditional.

The same condition should cover:

- the exact probability vector;
- the exact payoff;
- `payoff_dbl`.

Counting, pruning, and stability already use integer numerators and denominators. They do not require the public rational objects.
The former `materialize_vector` argument is now `materialize_output`.

The measured result is recorded in Current Checkpoint Evidence above.

## Use coposit's Integer Types Through Aliases

The former FracESSA and coposit integer wrappers were nearly identical. The generic string operations now live in coposit's
`integer`.

FracESSA's `fracessa/types.hpp` keeps the small FracESSA-owned double matrix beside zero-cost aliases for coposit's exact types:

```cpp
namespace fracessa::numeric {
using integer = coposit::integer;
using matrix_int = coposit::matrix_integer;
class matrix_dbl; // FracESSA's row-major binary64 scratch matrix.
}
```

These aliases introduce no object, copy, function call, or runtime layer. Retaining them permanently is reasonable. Replacing every
`fracessa::numeric::matrix_int` spelling with `coposit::matrix_integer` would create a large mechanical diff without a speed or
correctness benefit.

## Share The Generic LDLT, Not The Candidate-Specific LDLT

FracESSA's generic `fraction_free_ldlt.hpp` was removed in favor of coposit's generic implementation. Generic rank, inertia, solve,
and nullspace operations belong in coposit.

FracESSA's `fraction_free_ldlt_kkt.hpp` must remain. It is specialized for the hot candidate path:

- one right-hand side;
- border-reduced equilibrium equations;
- reusable workspaces;
- exact reduced-Hessian inertia;
- stability reuse.

Moving that solver into coposit would place ESS-specific candidate-solving behavior in a generic copositivity library and reverse the
intended conceptual dependency.

## Provide An Embeddable coposit Target

coposit's standalone build includes command-line programs, analysis companions, Python modules, and tests. Its embedded build does
not trigger those products.

coposit provides:

```text
coposit::core
    integer and matrix types
    parsers
    generic LDLT
    model-independent preprocessing headers

coposit::safe
    selected exact strict-copositivity implementation
```

It also exposes one small public function equivalent to:

```cpp
bool coposit::safe::is_strictly_copositive(const matrix_integer& matrix);
```

This function owns the selected preprocessing and Dickinson Final implementation. FracESSA must not know which coposit model source
file supplies the decision.

Embedding options disable unrelated products when coposit is a subproject:

```text
COPOSIT_BUILD_APPS=OFF
COPOSIT_BUILD_PYTHON=OFF
COPOSIT_BUILD_TESTS=OFF
```

Those options default to enabled when coposit is the top-level project and disabled when another project embeds it.

## Parser Behavior

coposit's parser already supports:

- integers;
- exact fractions;
- exact decimals;
- scientific notation;
- compact circular input;
- full upper-triangular input;
- Matrix Market input.

It retains the common positive denominator in `parsed_matrix`.

FracESSA uses this broader exact syntax. Decimal and scientific input become exact rational values; no floating-point parsing is
involved.

For the CLI, the simplest shared rule is the existing coposit behavior:

- an argument beginning like `dimension#...` is inline matrix input;
- another argument is treated as a file path;
- Matrix Market files are selected by their banner.

Python does not initially need another native file-reading API. Python can read a file with `Path.read_text()` and pass the contents
to the same parser. A convenience `Matrix.from_file()` should be added only if real use requires it.

## Move Fast Preparation Out Of The Exact Solver

Before migration, the precision-span calculation and normalized double conversion lived in `exact_candidate_solver`, while the fast
filter called back into the exact solver to prepare its own matrix. That ownership was backwards.

Precision-span classification and double conversion now live in `fast_candidate_filter`, which reads the shared integer game
directly.

The unused overload that included the common game denominator was deleted. The reduced-system fast path deliberately ignores a
common positive payoff scale.

## Do Not Share More Yet

The shared boundary should stop here. The following remain FracESSA-specific:

- `matrix_dbl`;
- candidate and ESS result types;
- circular support and bracelet generation;
- affine cyclic candidate reduction;
- fast-mode cutoffs and precision-span decisions;
- candidate-specific KKT LDLT;
- candidate CSV and JSON formatting;
- logging;
- the FracESSA test database.

The packed multiword support representation may eventually be shareable, but changing it now would enlarge an already substantial
migration and touch a proven hot path without an immediate need.

## Build And Release Consequences — Implemented

The release infrastructure now complements the Git submodule by:

- checking out submodules in every CLI, wheel, source-distribution, and documentation build that needs coposit;
- including the required coposit source in the PyPI source distribution;
- ensuring standalone source packages contain the pinned coposit source;
- no longer describing GitHub's automatically generated FracESSA tag archive as complete corresponding source, because an ordinary
  GitHub archive does not contain submodule contents; and
- publishing a complete source archive containing the exact coposit revision.

This is required for reproducible builds and clean GPL source distribution.

## Checkpoint Order And Status

### 1. coposit foundation

- [x] Return `parsed_matrix` from both compact and Matrix Market parsers.
- [x] Preserve the common denominator and compact-circular flag.
- [x] Add generic integer string conversion.
- [x] Add `coposit::core` and `coposit::safe`.
- [x] Make coposit embeddable without building all applications, Python modules, and tests.
- [x] Add focused parser and embedded-library checks.
- [x] Commit coposit first.

### 2. FracESSA representation migration

- [x] Update the coposit submodule pointer.
- [x] Use coposit's parser and exact-type aliases.
- [x] Store one integer game plus its denominator.
- [x] Pass references to the exact solver and convert the fast matrix directly from the shared integer game.
- [x] Change candidate vectors to `std::vector<fraction>`.
- [x] Materialize all rational output only when requested.
- [x] Delete `matrix_fraction.hpp`, the duplicate parser, and rational-to-integer conversion.
- [x] Slim `fraction.hpp` to the remaining output operations.
- [x] Preserve the preceding copositivity algorithm through the representation benchmark.
- [x] Verify outputs and benchmark against the preceding FracESSA revision.

### 3. Copositivity integration

- [x] Replace FracESSA's complete copositivity preprocessing and final decision with one `coposit::safe` call, rather than putting
  coposit behind FracESSA's former preliminary paths.
- [x] Verify the complete correctness corpus and documented difficult matrices.
- [x] Benchmark this algorithm change separately.
- [x] Delete FracESSA's obsolete copositivity preprocessing, Hadeler implementation, and generic LDLT code only after verification.

Separating the representation migration from the algorithm replacement preserves a useful performance baseline and makes regressions
much easier to diagnose.

## Final Recommendation

Eliminate `matrix_frc`, retain `matrix_dbl`, keep a minimal output-only `fraction`, and move generic exact-integer infrastructure to
coposit. Represent parsed input explicitly as `{integer matrix, positive denominator, compact-circular flag}`. Do not hide the
denominator inside the integer or matrix class, and do not move FracESSA's specialized candidate LDLT into coposit.
