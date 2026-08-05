# Exact Stability Schur-Complement Plan

Status: implementation-ready plan; no source changes have been made.

## 1. Goal

When an exact candidate has outside best replies, reuse the exact fraction-free factorization of its reduced Hessian instead of:

1. constructing the complete rational Bomze matrix `Bee`;
2. factoring that larger matrix for the positive-definite shortcut; and
3. recursively eliminating the unrestricted support coordinates one at a time.

The replacement solves the already-factored reduced Hessian for one new matrix right-hand side and constructs the smaller Schur
complement directly. The existing exact positive-definiteness and strict-copositivity routines then decide stability on that smaller
matrix.

Correctness remains the first requirement. The optimization is acceptable only if every candidate and ESS decision is unchanged.

## 2. What the current solver preserves

For a support of size `s`, `find_candidate_safe::find()` constructs the integer-scaled reduced Hessian

$$
\widehat H=dH\in\mathbb Z^{r\times r},
\qquad r=s-1,
$$

where `d > 0` is the common game denominator. It then calls
`fraction_free_ldlt_workspace::solve_inplace()`.

That call overwrites three different objects in three different ways:

| Object | State after the call |
|---|---|
| `reduced_system_` | Contains the complete fraction-free factorization of $\widehat H$ in its lower triangle and diagonal. |
| `right_hand_side_` | Contains the forward-transformed candidate right-hand side. Its original values are no longer available here. |
| `solution_numerators_` | Contains the candidate solution numerators after back substitution. |

The back-substitution loop reads `reduced_system_`; it does **not** overwrite it. Therefore, the factorization is still valid when
`check_stability()` is called immediately after the successful exact candidate solve.

The original right-hand side is the disposable object. This is not a problem for the proposed stability solve: the new matrix
right-hand side $\widehat G$ can be rebuilt from the existing reduced-entry cache, and those cached entries remain available after
the solve.

The recommended boundary is therefore **not** a getter that returns or copies `reduced_system_`. The factor representation should
remain private. The workspace receives the new matrix right-hand side $\widehat G$, reads the retained factor in place, and overwrites
only that disposable right-hand-side matrix with its solution numerators.

## 3. Block matrix and dimensions

Let:

- $I$ be the candidate support;
- $J$ be the extended support, including unused pure strategies that tie the candidate payoff;
- $m\in I$ be the same reference strategy used by the candidate solve;
- $U=I\setminus\{m\}$;
- $K=J\setminus I$.

Define

$$
r=|U|=s-1,
\qquad
k=|K|=|J|-|I|.
$$

In the coordinate order $(U,K)$, the reduced Hessian of the extended support is

$$
R=
\begin{pmatrix}
H&G\\
G^T&Q
\end{pmatrix},
$$

with

$$
H\in\mathbb Q^{r\times r},
\qquad
G\in\mathbb Q^{r\times k},
\qquad
Q\in\mathbb Q^{k\times k}.
$$

The Bomze matrix is `Bee = -2R`. Since candidate stability reaches this path only after proving $H\prec0$, its support block
$P=-2H$ is positive definite.

The final matrix that must be strictly copositive is the Schur complement

$$
S=-2\left(Q-G^TH^{-1}G\right)\in\mathbb Q^{k\times k}.
$$

No inverse will be constructed. We solve

$$
HX=G
$$

with the retained factorization and then form

$$
S=-2(Q-G^TX).
$$

## 4. Keep the construction integral

The integer game already provides

$$
\widehat H=dH,
\qquad
\widehat G=dG,
\qquad
\widehat Q=dQ.
$$

Let the reused fraction-free solve return an integer matrix $N$ and positive denominator $\delta$ satisfying

$$
\widehat HN=\widehat G\,\delta.
$$

Because the same positive factor `d` occurs on both sides,

$$
\frac{N}{\delta}=H^{-1}G.
$$

Construct the integer matrix

$$
\widehat S=-2\left(\delta\widehat Q-\widehat G^TN\right).
$$

Then

$$
\widehat S=d\delta S.
$$

Both `d` and $\delta$ are positive. Multiplication by `d * delta` changes neither positive definiteness nor strict copositivity.
The existing rational routines can therefore receive $\widehat S$ as a rational matrix with denominator one. No rational matrix
division is needed during Schur-complement construction.

## 5. Why the reused solve can be narrow

The existing exact solver can handle a general nonsingular symmetric matrix. If it encounters a zero active diagonal, it records a
symmetric swap or coordinate addition and later restores the original solution coordinates.

The stability reuse method does not need that general mechanism:

1. it is called only after the existing inertia result proved $H\prec0$;
2. every leading principal submatrix of a negative-definite matrix is negative definite and therefore nonsingular;
3. every active Bareiss diagonal pivot is consequently nonzero; and
4. the candidate factorization recorded no coordinate swap or coordinate addition.

The new method should therefore be named for this restricted contract and assert that `operation_count_ == 0`. It should not grow
into a second general solver or replay pivot operations that cannot occur on this path.

## 6. Proposed source changes

Future implementation should change only these existing source files:

| File | Minimal change |
|---|---|
| `cpp/include/linalg/flint_style_fraction_free_ldlt.hpp` | Add one multi-column solve that reuses an unmodified negative-definite factorization. |
| `cpp/include/fracessa/find_candidate_safe.hpp` | Declare the Schur builder and one reusable integer matrix for $N$. |
| `cpp/src/find_candidate_safe.cpp` | Build $\widehat G$ and $\widehat Q$ from the existing integer game/cache, reuse the factorization, and write $\widehat S$. |
| `cpp/include/fracessa/fracessa.hpp` | Rename `bee_matrix_` to the accurate `stability_matrix_`. |
| `cpp/src/checkstab.cpp` | Replace complete `Bee` construction and recursive elimination with the Schur builder and checks on the smaller matrix. |
| `cpp/tests/test_linear_solver.cpp` | Add one focused multi-right-hand-side factor-reuse test. |

No parser, CLI, Python, candidate-output, database, support-generator, fast-search, or copositivity source file needs to change.

## 7. Proposed factor-reuse method

The method should overwrite the supplied right-hand-side matrix with its solution numerators. This avoids adding a second temporary
matrix. The caller already has the positive denominator from the original factorization.

Proposed interface:

```cpp
// Solve factored_system * X = right_hand_sides * denominator for several columns.
// Preconditions: factored_system is the retained negative-definite factorization from solve_inplace(),
// and that factorization required no coordinate operations. right_hand_sides becomes X.
void solve_factored_negative_definite_inplace(
    matrix_int& right_hand_sides,
    const integer& denominator,
    const matrix_int& factored_system) const;
```

Core implementation shape:

```cpp
assert(operation_count_ == 0);
assert(factored_system.rows() == factored_system.cols());
assert(right_hand_sides.rows() == factored_system.rows());

const size_t dimension = factored_system.rows();
const size_t right_hand_side_count = right_hand_sides.cols();

// Replay only the right-hand-side part of the fraction-free forward elimination.
for (size_t pivot_position = 0; pivot_position + 1 < dimension; ++pivot_position) {
    const auto pivot = factored_system(pivot_position, pivot_position);

    for (size_t row = pivot_position + 1; row < dimension; ++row) {
        for (size_t column = 0; column < right_hand_side_count; ++column) {
            auto value = right_hand_sides(row, column);
            value.set_product(value, pivot);
            value.submul(factored_system(row, pivot_position), right_hand_sides(pivot_position, column));
            if (pivot_position > 0) {
                const auto previous_pivot = factored_system(pivot_position - 1, pivot_position - 1);
                if (!previous_pivot.is_one()) value.divide_exact(previous_pivot);
            }
        }
    }
}

// Back substitution is safe in place: rows below `row` already contain solved numerators,
// and no later step needs the old transformed value at those rows.
for (size_t column = 0; column < right_hand_side_count; ++column) {
    for (size_t row = dimension; row-- > 0;) {
        auto numerator = right_hand_sides(row, column);
        numerator.set_product(numerator, denominator);
        for (size_t solved_row = row + 1; solved_row < dimension; ++solved_row) {
            numerator.submul(factored_system(solved_row, row), right_hand_sides(solved_row, column));
        }
        numerator.divide_exact(factored_system(row, row));
    }
}
```

The first version should use the normal FLINT `fmpz_mul`, `fmpz_submul`, and `fmpz_divexact` operations already wrapped by
`integer::reference`. It should not duplicate the specialized immediate-integer machinery before a benchmark demonstrates that this
rare stability path needs it.

## 8. Proposed Schur builder

Add one method to the exact candidate object because it owns all required private state: the integer game, retained factorization,
common denominator, and reduced-entry cache.

```cpp
// Build a positive integer multiple of the exact stability Schur complement.
// Called immediately after find() succeeds for the same support and proves H negative definite.
void build_scaled_stability_schur(
    bitset64 support,
    bitset64 extended_support,
    linalg::matrix_frc& result);
```

Implementation outline:

```cpp
const bitset64 outside_best_replies = bs64::subtract(extended_support, support);

uint8_t support_indices[bs64::kMaxBitsetDimension];
uint8_t outside_indices[bs64::kMaxBitsetDimension];
const size_t support_size = bs64::extract_set_indices(support, dimension_, support_indices);
const size_t outside_size = bs64::extract_set_indices(outside_best_replies, dimension_, outside_indices);
const size_t reduced_dimension = support_size - 1;
const size_t reference = support_indices[0];

stability_solution_numerators_.resize(reduced_dimension, outside_size);

// Before the solve this matrix is G-hat. Afterwards it is N.
for (size_t row = 0; row < reduced_dimension; ++row) {
    const size_t i = support_indices[row + 1];
    for (size_t column = 0; column < outside_size; ++column) {
        stability_solution_numerators_(row, column) = reduced_entry(reference, i, outside_indices[column]);
    }
}

if (reduced_dimension > 0) {
    ffldlt_workspace_.solve_factored_negative_definite_inplace(
        stability_solution_numerators_, solution_denominator_, reduced_system_);
}

if (result.rows() != outside_size || result.cols() != outside_size) {
    result = linalg::matrix_frc(outside_size, outside_size);
}

linalg::integer scaled_entry;
linalg::integer one(1);
for (size_t row = 0; row < outside_size; ++row) {
    for (size_t column = 0; column <= row; ++column) {
        scaled_entry.set_product(
            solution_denominator_, reduced_entry(reference, outside_indices[row], outside_indices[column]));

        for (size_t inner = 0; inner < reduced_dimension; ++inner) {
            const size_t i = support_indices[inner + 1];
            scaled_entry.submul(
                reduced_entry(reference, i, outside_indices[row]), stability_solution_numerators_(inner, column));
        }

        scaled_entry.multiply(2);
        scaled_entry.negate();
        result(row, column).set_ratio(scaled_entry, one);
        if (row != column) result(column, row) = result(row, column);
    }
}
```

`reduced_entry(reference, row, column)` should be a small private cache accessor implementing the same formula already used by
`build_reduced_system()`:

$$
d\left(A_{ij}-A_{mj}-A_{im}+A_{mm}\right).
$$

The accessor should populate the existing dense cache on the first request and return the stored integer thereafter. It need not
replace the current row-specialized loop in `build_reduced_system()`: changing that proven hot loop is unnecessary.

The solver overwrites $\widehat G$ with $N$. This does not require a second $r\times k$ matrix because every original
$\widehat G_{ia}$ was inserted through `reduced_entry()` and remains in the existing cache for the later product
$\widehat G^TN$.

For a pure support, `reduced_dimension == 0`. The solve is skipped, `solution_denominator_` is already one, and the formula reduces
correctly to

$$
\widehat S=-2\widehat Q.
$$

This explicit case also prevents accidental reuse of a stale nonempty factorization left by an earlier support.

## 9. Proposed `check_stability()` tail

The early decisions remain unchanged:

1. reject when the retained reduced Hessian is not negative definite;
2. accept when the extended support equals the support.

Only the remaining outside-best-reply path changes:

```cpp
find_candidate_safe_.build_scaled_stability_schur(
    candidate_.support, candidate_.extended_support, stability_matrix_);
auto& schur = stability_matrix_;

if (schur.is_positive_definite()) {
    candidate_.stability = "T_pd_frc";
    candidate_.is_ess = true;
    return;
}

if (kay_size <= 1) {
    candidate_.stability = "F_not_pd_kay_0_1";
    candidate_.is_ess = false;
    return;
}

candidate_.is_ess = linalg::is_strictly_copositive(schur);
candidate_.stability = candidate_.is_ess ? "T_copos" : "F_not_copos";
```

The current complete `Bee` construction, `r`-step recursive reduction, rolling rational matrices, and partial-pivot checks are then
deleted. They are mathematically represented by the one block Schur complement. Existing externally visible stability labels remain
unchanged. Diagnostic logging should call the matrix `scaled stability Schur complement`, not `Bee`.

## 10. Correctness checks

### 10.1 Focused factor-reuse test

Use

$$
H=
\begin{pmatrix}
-2&1\\
1&-3
\end{pmatrix},
\qquad
G=
\begin{pmatrix}
1&2\\
3&4
\end{pmatrix}.
$$

The matrix $H$ is negative definite, `det(H) = 5`, and

$$
H^{-1}G=
\frac15
\begin{pmatrix}
-6&-10\\
-7&-10
\end{pmatrix}.
$$

Factor `H` once through the existing solver, call the new two-column reuse method, and verify the positive denominator and all four
numerators exactly. This test catches incorrect previous-pivot division, column handling, and in-place back substitution.

### 10.2 Stability regression coverage

Run the existing complete C++ and matrix-correctness suites. In particular, preserve examples that reach all relevant final states:

- positive-definite acceptance (`T_pd_frc`);
- one outside best reply rejected after positive definiteness fails (`F_not_pd_kay_0_1`);
- strict-copositivity acceptance (`T_copos`; canonical database matrix 29 contains such a candidate);
- strict-copositivity rejection (`F_not_copos`; canonical database matrices 4 and 7 contain examples); and
- rejection before Schur construction because $H$ is not negative definite (`F_not_part_copos`).

Candidate count, ESS count, support structure, extended supports, exact vectors, exact payoffs, and stability labels must match the
current baseline. Candidate IDs are not a mathematical requirement, but this change does not alter support enumeration, so they
should also remain unchanged.

### 10.3 Direct equivalence check during development

In an experiment-only build, retain the old stability function under a different name and compare old and new decisions for every
stored matrix that reaches the outside-best-reply path. Remove that duplicate code before the production patch is accepted.

## 11. Performance measurement

Benchmark only after correctness passes. Compare current code with the Schur version on:

1. matrices whose candidates almost always have `extended_support == support`, to confirm no measurable regression in the common path;
2. matrices 4, 7, 19, 20, 29, and 32, which exercise outside best replies and copositivity paths; and
3. representative small, medium, and large dimensions, with circular and non-circular games where available.

Expected work on the changed path:

| Current | Proposed |
|---|---|
| Construct a rational $(r+k)\times(r+k)$ `Bee`. | Build one integer $r\times k$ right-hand side. |
| Rational positive-definite factorization of size $r+k$. | Reuse the existing factorization of $H$. |
| Up to `r` rational rank-one reductions with new matrices. | One multi-RHS solve, $O(r^2k)$ integer operations. |
| Final strict copositivity on a $k\times k$ matrix. | Same final strict copositivity on a $k\times k$ matrix. |

The retained implementation should show a repeatable end-to-end improvement on the affected matrices and no meaningful regression
on the common no-outside-reply path. If it does not, keep this document as a rejected experiment and revert the source patch.

## 12. Deliberately excluded work

Do not include any of the following in the first implementation:

- explicit construction of $H^{-1}$;
- a generic reusable factorization class;
- support for replaying symmetric swaps or coordinate additions;
- a second copy of $\widehat G$;
- a new integer copositivity implementation;
- a specialized immediate-integer fast path for the reused right-hand-side solve;
- changes to the existing candidate factorization loop; or
- unrelated solver, parser, Python, CLI, or support-generator changes.

These additions are unnecessary for correctness and would make the patch harder to audit. Add one only after a measured result shows
that this smaller implementation is insufficient.
