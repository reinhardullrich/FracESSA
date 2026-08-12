# Coposit Integration And Exact-Type Ownership

Status: proposed architecture; no implementation has started.

## Decision

The general direction is correct, with two important qualifications:

1. The common input denominator must belong to the parsed input as a whole. It must not be hidden inside `integer` or
   `matrix_integer`.
2. The representation cleanup and the replacement of FracESSA's strict-copositivity algorithm must be separate measured changes.
   Otherwise, a timing or behavior change could not be attributed reliably.

The intended dependency direction is:

```text
FracESSA -> Coposit
```

Coposit must never depend on FracESSA. A third common repository is unnecessary.

## Final Ownership

| Component | Decision | Owner |
|---|---|---|
| `matrix_dbl` | Keep unchanged | FracESSA |
| `matrix_frc` | Remove | None |
| `fraction` | Keep as a small output-only rational value | FracESSA |
| `integer` | Use Coposit's implementation | Coposit |
| `matrix_int` | Alias Coposit's `matrix_integer` | Coposit |
| Generic fraction-free LDLT | Use Coposit's implementation | Coposit |
| Candidate-specific fraction-free KKT LDLT | Keep | FracESSA |
| Exact matrix parsers | Use Coposit's parsers | Coposit |
| Strict-copositivity implementation | Call one public Coposit function | Coposit |

## Keep `matrix_dbl`

`matrix_dbl` is a small row-major dense matrix used by the fast candidate filter for:

- the complete converted double game;
- the reusable reduced candidate system;
- in-place Bunch-Kaufman factorization.

Replacing it with `std::vector<double>` would not eliminate meaningful code. It would merely distribute expressions such as
`row * dimension + column` throughout the numerical solver.

Coposit does not need a corresponding dense double factorization matrix. A shared or templated matrix class would therefore add an
abstraction without removing duplication.

The fast-filter constructor currently receives a complete `matrix_frc` but uses only its dimension. After the input migration, it
should receive the dimension and prepare its double matrix directly from the shared integer matrix.

## Remove `matrix_frc`

`matrix_frc` currently combines four unrelated responsibilities:

| Current responsibility | Replacement |
|---|---|
| Parsed game storage | Coposit `parsed_matrix` |
| Compact circular-matrix expansion | Coposit parser |
| Matrix retained for logging and symmetry detection | Shared integer matrix plus denominator |
| Candidate probability vector | `std::vector<fraction>` |

The current input path is:

```text
text
  -> complete matrix_frc
  -> copy retained by FracESSA
  -> temporary FLINT rational matrix
  -> denominator clearing
  -> complete matrix_int
```

The replacement is:

```text
text
  -> Coposit parser
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
Pretty diagnostic output is a FracESSA responsibility and must not become a method of Coposit's generic integer matrix.

### Circular symmetry

Exact symmetry detection can compare entries of `Z` directly because the entire matrix has one common positive denominator:

$$
A_{ij}=A_{kl}\quad\Longleftrightarrow\quad Z_{ij}=Z_{kl}.
$$

The affine-symmetry detector should therefore receive the integer matrix. No rational matrix is required.

## Represent The Parsed Input Explicitly

Coposit should return one small aggregate:

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

Coposit uses `parsed.matrix` and ignores the denominator because copositivity is invariant under multiplication by a positive
scalar. FracESSA uses:

- `matrix` for exact candidate and stability calculations;
- `matrix` for fast double conversion;
- `denominator` when constructing exact payoffs;
- `compact_circular` when selecting circular support generation.

Candidate probabilities do not need the game denominator. A common positive scaling of the payoff matrix does not change the
equilibrium vector.

## Store The Complete Integer Game Once

Currently the analyzer owns the rational game and the exact solver owns another complete integer game. The new ownership should be:

```text
basic_fracessa
  owns parsed_matrix { Z, d, compact_circular }

  exact_candidate_solver
    holds const references to Z and d

  fast_candidate_filter
    reads Z once and owns its converted matrix_dbl
```

The parser result should be moved into the analyzer. `matrix_integer` already has an inexpensive move operation implemented by
swapping FLINT matrix storage.

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
std::vector<linalg::fraction> vector;
```

That reflects their actual meaning as an `n`-element vector rather than an `n x 1` matrix.

Coposit should not acquire this rational output class. Its maintained numerical core remains integer-only.

## Materialize Rational Output Only When Requested

The exact solver currently constructs the exact rational payoff and `payoff_dbl` for every successful candidate, even when candidate
output was not requested. Only the dense probability vector is conditional.

The same condition should cover:

- the exact probability vector;
- the exact payoff;
- `payoff_dbl`.

Counting, pruning, and stability already use integer numerators and denominators. They do not require the public rational objects.
The current `materialize_vector` argument should therefore become `materialize_output`.

This may improve runs with many candidates and no candidate output, but the improvement must be measured rather than assumed.

## Use Coposit's Integer Types Through Aliases

The FracESSA and Coposit integer wrappers are already nearly identical. FracESSA additionally has `set_string()` and `to_string()`;
these are generic integer operations and should move into Coposit's `integer`.

FracESSA's existing headers can then become zero-cost aliases:

```cpp
namespace linalg {
using integer = coposit::integer;
using matrix_int = coposit::matrix_integer;
}
```

These aliases introduce no object, copy, function call, or runtime layer. Retaining them permanently is reasonable. Replacing every
`linalg::matrix_int` spelling with `coposit::matrix_integer` would create a large mechanical diff without a speed or correctness
benefit.

## Share The Generic LDLT, Not The Candidate-Specific LDLT

FracESSA's generic `fraction_free_ldlt.hpp` should eventually be removed in favor of Coposit's more capable generic implementation.
Generic rank, inertia, solve, and nullspace operations belong in Coposit.

FracESSA's `fraction_free_ldlt_kkt.hpp` must remain. It is specialized for the hot candidate path:

- one right-hand side;
- border-reduced equilibrium equations;
- reusable workspaces;
- exact reduced-Hessian inertia;
- stability reuse.

Moving that solver into Coposit would place ESS-specific candidate-solving behavior in a generic copositivity library and reverse the
intended conceptual dependency.

## Provide An Embeddable Coposit Target

Coposit currently builds its command-line programs, analysis companions, Python modules, and potentially its tests. FracESSA must
not trigger all of those merely by adding the Coposit source directory.

Coposit should provide:

```text
Coposit::core
    integer and matrix types
    parsers
    generic LDLT
    model-independent preprocessing headers

Coposit::safe
    selected exact strict-copositivity implementation
```

It should also expose one small public function conceptually equivalent to:

```cpp
bool coposit::safe::is_strictly_copositive(const matrix_integer& matrix);
```

This function owns the selected preprocessing and Dickinson Final implementation. FracESSA must not know which Coposit model source
file supplies the decision.

Embedding options should disable unrelated products when Coposit is a subproject:

```text
COPOSIT_BUILD_APPS=OFF
COPOSIT_BUILD_PYTHON=OFF
COPOSIT_BUILD_TESTS=OFF
```

Those options should default to enabled when Coposit is built as the top-level project and disabled when another project embeds it.

## Parser Behavior

Coposit's parser already supports:

- integers;
- exact fractions;
- exact decimals;
- scientific notation;
- compact circular input;
- full upper-triangular input;
- Matrix Market input.

It already computes a common positive denominator but currently discards it. It should return that denominator in `parsed_matrix`.

FracESSA should adopt the broader exact syntax. Decimal and scientific input still become exact rational values; no floating-point
parsing is involved.

For the CLI, the simplest shared rule is the existing Coposit behavior:

- an argument beginning like `dimension#...` is inline matrix input;
- another argument is treated as a file path;
- Matrix Market files are selected by their banner.

Python does not initially need another native file-reading API. Python can read a file with `Path.read_text()` and pass the contents
to the same parser. A convenience `Matrix.from_file()` should be added only if real use requires it.

## Move Fast Preparation Out Of The Exact Solver

The precision-span calculation and normalized double conversion currently live in `exact_candidate_solver`, while the fast filter
calls back into the exact solver to prepare its own matrix. This ownership is backwards.

Precision-span classification and double conversion are fast-mode preparation. They should move into `fast_candidate_filter`, which
will read the shared integer game directly.

The current overload that includes the common game denominator appears unused. The reduced-system fast path deliberately ignores a
common positive payoff scale. That unused overload can be deleted during the migration.

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

## Build And Release Consequences

The Git submodule alone is insufficient. Before FracESSA depends on Coposit, the release infrastructure must also:

- check out submodules in every CLI, wheel, source-distribution, and documentation build that needs Coposit;
- include the required Coposit source in the PyPI source distribution;
- ensure standalone source packages contain the pinned Coposit source;
- stop describing GitHub's automatically generated FracESSA tag archive as complete corresponding source, because an ordinary
  GitHub archive does not contain submodule contents;
- either publish a complete source archive or link explicitly to the exact Coposit revision alongside the FracESSA source.

This is required for reproducible builds and clean GPL source distribution.

## Recommended Implementation Order

### 1. Coposit foundation

- Return `parsed_matrix` from both compact and Matrix Market parsers.
- Preserve the common denominator and compact-circular flag.
- Add generic integer string conversion.
- Add `Coposit::core` and `Coposit::safe`.
- Make Coposit embeddable without building all applications, Python modules, and tests.
- Add focused parser and embedded-library checks.
- Commit Coposit first.

### 2. FracESSA representation migration

- Update the Coposit submodule pointer.
- Use Coposit's parser and exact-type aliases.
- Store one integer game plus its denominator.
- Pass references to the exact solver and convert the fast matrix directly from the shared integer game.
- Change candidate vectors to `std::vector<fraction>`.
- Materialize all rational output only when requested.
- Delete `matrix_fraction.hpp`, the duplicate parser, and rational-to-integer conversion.
- Slim `fraction.hpp` to the remaining output operations.
- Preserve the current copositivity algorithm temporarily.
- Verify outputs and benchmark against the preceding FracESSA revision.

### 3. Copositivity integration

- Replace FracESSA's final Hadeler call with `Coposit::safe`.
- Verify the complete correctness corpus and documented difficult matrices.
- Benchmark this algorithm change separately.
- Delete FracESSA's obsolete generic copositivity and generic LDLT code only after verification.

Separating the representation migration from the algorithm replacement preserves a useful performance baseline and makes regressions
much easier to diagnose.

## Final Recommendation

Eliminate `matrix_frc`, retain `matrix_dbl`, keep a minimal output-only `fraction`, and move generic exact-integer infrastructure to
Coposit. Represent parsed input explicitly as `{integer matrix, positive denominator, compact-circular flag}`. Do not hide the
denominator inside the integer or matrix class, and do not move FracESSA's specialized candidate LDLT into Coposit.
