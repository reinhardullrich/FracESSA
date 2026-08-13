# Exact Stability Schur Complement

Status: current mathematical design and implementation record. Implemented and verified on 2026-08-05; later exact-copositivity
shortcuts changed only the decisions made after this reduction.

## 1. Implemented result

When an exact candidate has outside best replies, reuse the exact fraction-free factorization of its reduced Hessian instead of:

1. constructing the complete rational Bomze matrix `Bee`;
2. factoring that larger matrix for the positive-definite shortcut; and
3. recursively eliminating the unrestricted support coordinates one at a time.

The replacement solves the already-factored reduced Hessian for one new matrix right-hand side and constructs the smaller Schur
complement directly. The current implementation passes that exact integer matrix to `coposit::safe` for the final strict-copositivity
decision.

Correctness remains the first requirement. The optimization is acceptable only if every candidate and ESS decision is unchanged.

## 2. What the current solver preserves

For a support of size `s`, `exact_candidate_solver::find()` constructs the integer-scaled reduced Hessian

$$
\widehat H=dH\in\mathbb Z^{r\times r},
\qquad r=s-1,
$$

where `d > 0` is the common game denominator. It then calls
`kkt_fraction_free_ldlt_workspace::solve_inplace()`.

That call overwrites three different objects in three different ways:

| Object | State after the call |
|---|---|
| `reduced_system_` | Contains the complete fraction-free factorization of $\widehat H$ in its lower triangle and diagonal. |
| `right_hand_side_` | Contains the forward-transformed candidate right-hand side. Its original values are no longer available here. |
| `solution_numerators_` | Contains the candidate solution numerators after back substitution. |

The back-substitution loop reads `reduced_system_`; it does **not** overwrite it. Therefore, the factorization is still valid when
`check_stability()` is called immediately after the successful exact candidate solve.

The original right-hand side is the disposable object. This is not a problem for the retained stability solve: the new matrix
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
- $K=J\setminus I$ be the set of outside best replies.

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

Bomze's final matrix is the Schur complement

$$
B^{(r)}=-2\left(Q-G^TH^{-1}G\right)\in\mathbb Q^{k\times k}.
$$

Production stores the half-scaled matrix

$$
C=-\left(Q-G^TH^{-1}G\right).
$$

Since $B^{(r)}=2C$ and strict copositivity is unchanged by multiplication by a positive scalar, testing $C$ is equivalent and avoids
one unnecessary multiplication of every entry.

No inverse is constructed. Production solves

$$
HX=G
$$

with the retained factorization and then forms

$$
C=-(Q-G^TX).
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
\widehat C=-\left(\delta\widehat Q-\widehat G^TN\right).
$$

Then

$$
\widehat C=d\delta C.
$$

Both `d` and $\delta$ are positive. Multiplication by `d * delta` changes neither positive definiteness nor strict copositivity.
The exact integer routines therefore test $\widehat C$ directly. No rational matrix construction or division is needed.

## 5. Why the reused solve can be narrow

The existing exact solver can handle a general nonsingular symmetric matrix. If it encounters a zero active diagonal, it records a
symmetric swap or coordinate addition and later restores the original solution coordinates.

The stability reuse method does not need that general mechanism:

1. it is called only after the existing inertia result proved $H\prec0$;
2. every leading principal submatrix of a negative-definite matrix is negative definite and therefore nonsingular;
3. every active Bareiss diagonal pivot is consequently nonzero; and
4. the candidate factorization recorded no coordinate swap or coordinate addition.

The reuse method is named for this restricted contract and asserts that `operation_count_ == 0`. It must not grow
into a second general solver or replay pivot operations that cannot occur on this path.

## 6. Implemented source changes

The original patch changed only these source files:

| File | Minimal change |
|---|---|
| `cpp/include/fracessa/fraction_free_ldlt_kkt.hpp` | Added one multi-column solve that reuses an unmodified negative-definite factorization. |
| `cpp/include/fracessa/exact_candidate_solver.hpp` | Declared the reduced-$B$ builder and one reusable integer matrix for $N$. |
| `cpp/src/exact_candidate_solver.cpp` | Builds $\widehat G$ and $\widehat Q$ from the integer game/cache, reuses the factorization, and writes $\widehat C$. |
| `cpp/include/fracessa/fracessa.hpp` | Renamed `bee_matrix_` to the accurate `scaled_reduced_b_`. |
| `cpp/src/check_stability.cpp` | Replaced complete `Bee` construction and recursive elimination with the reduced-$B$ builder and checks on the smaller matrix. |
| `cpp/tests/test_linear_solver.cpp` | Added one focused multi-right-hand-side factor-reuse test. |

No parser, CLI, Python, candidate-output, database, support-generator, fast-search, or copositivity source file changed in the original
patch.

## 7. Factor-reuse method

The method overwrites the supplied right-hand-side matrix with its solution numerators. This avoids a second temporary
matrix. The caller already has the positive denominator from the original factorization.

Implemented interface:

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

The implementation uses the normal FLINT `fmpz_mul`, `fmpz_submul`, and `fmpz_divexact` operations already wrapped by
`integer::reference`. It does not duplicate the specialized immediate-integer machinery used by hotter paths.

## 8. Schur builder

The exact candidate object owns the builder because it already owns the integer game, retained factorization, common denominator, and
reduced-entry cache. One-word and multiword overloads extract indices and call the same internal implementation:

```cpp
size_t build_scaled_reduced_b(
    fracessa::support::bitset support,
    fracessa::support::bitset outside_best_replies,
    fracessa::numeric::matrix_int& result);
size_t build_scaled_reduced_b(
    const fracessa::support::bitset_multiword& support,
    const fracessa::support::bitset_multiword& outside_best_replies,
    fracessa::numeric::matrix_int& result);
```

The shared builder fills $\widehat G$, replays the retained factorization to obtain $N$, and writes only one triangle of
$\widehat C=-(\delta\widehat Q-\widehat G^TN)$ before mirroring it. The result is an integer matrix; the discarded factor two is a
positive global scale and has no effect on either exact decision.

`reduced_entry(reference, row, column)` is a small private cache accessor implementing the same formula used by
`build_reduced_system()`:

$$
d\left(A_{ij}-A_{mj}-A_{im}+A_{mm}\right).
$$

The accessor populates the existing dense cache on first request and returns the stored integer thereafter. It does not
replace the current row-specialized loop in `build_reduced_system()`: changing that proven hot loop is unnecessary.

The solver overwrites $\widehat G$ with $N$. This does not require a second $r\times k$ matrix because every original
$\widehat G_{ia}$ was inserted through `reduced_entry()` and remains in the existing cache for the later product
$\widehat G^TN$.

For a pure support, `reduced_dimension == 0`. The solve is skipped, `solution_denominator_` is already one, and the formula reduces
correctly to

$$
\widehat C=-\widehat Q.
$$

This explicit case also prevents accidental reuse of a stale nonempty factorization left by an earlier support.

## 9. Current `check_stability()` flow

The early decisions remain unchanged:

1. reject when the retained reduced Hessian is not negative definite;
2. accept when the extended support equals the support.

The remaining path builds `scaled_reduced_b_` once and passes it to `coposit::safe`. Coposit owns the exact prechecks,
negative-entry connected-component reduction, and finite Dickinson certificate traversal. FracESSA neither duplicates those
decisions nor exposes their internal route in its stability label.

The former complete `Bee` construction, $r$ recursive reductions, rolling rational matrices, and partial-pivot checks are gone. The
single block Schur complement represents all of those operations. The stored integer result is called the `scaled reduced B matrix`:
it is a positive multiple of Bomze's final $B^{(r)}$ and therefore has exactly the same strict-copositivity decision.

## 10. Correctness checks

### 10.1 Focused factor-reuse test

Use the non-diagonal negative-definite system

$$
H=
\begin{pmatrix}
-4&-1&-1\\
-1&-3&-1\\
-1&-1&-2
\end{pmatrix},
\qquad
G=
\begin{pmatrix}
-12&-18\\
-15&-20\\
-14&-18
\end{pmatrix}.
$$

The matrix $H$ is negative definite, $\det(H)=-17$, and

$$
H^{-1}G=
\begin{pmatrix}
1&2\\
3&4\\
5&6
\end{pmatrix}.
$$

The implemented test factors `H` once through the existing solver, calls the new two-column reuse method, and verifies denominator
17 and numerator matrix $17H^{-1}G$ exactly. Three rows are necessary here: unlike a 2-by-2 example, this reaches the Bareiss step
that divides by the previous pivot. The same test also covers column handling and in-place back substitution.

### 10.2 Current stability regression coverage

The maintained tests preserve examples that reach the relevant final states:

- pure-strategy acceptance (`T_pure_ess`);
- reduced-Hessian negative-definite acceptance (`T_reduced_hessian_nd`);
- rejection before Schur construction because $H$ is not negative definite (`F_reduced_hessian_not_nd`);
- final Coposit acceptance and rejection (`T_copos` and `F_not_copos`).

Candidate count, ESS count, support structure, extended supports, exact vectors, exact payoffs, and stability labels must match the
current baseline. Candidate IDs are not a mathematical requirement, but this change does not alter support enumeration, so they
should also remain unchanged.

### 10.3 Promotion-time direct equivalence check

The implementation was compared row for row with the canonical database on all 309 analyzed matrices with an outside best reply
candidate and a sub-second safe baseline. All 7,808 stored representative candidates matched in candidate ID, exact vector, support,
extended support, multiplier, ESS decision, stability label, and exact payoff. Their multipliers represent 14,659 candidates. This
includes matrices 4, 5, 6, 7, 19, 20, 29, and 32, which together exercised acceptance, rejection, the zero-dimensional support
block, outside best replies, and rejection before Schur construction under the labels in force at promotion time.

## 11. Promotion-time performance evidence

The retained comparison covered:

1. matrices whose candidates almost always have `extended_support == support`, to confirm no measurable regression in the common path;
2. matrices 4, 7, 19, 20, 29, and 32, which exercise outside best replies and copositivity paths; and
3. representative small, medium, and large dimensions, with circular and non-circular games where available.

Expected work on the changed path:

| Previous | Implemented |
|---|---|
| Construct a rational $(r+k)\times(r+k)$ `Bee`. | Build one integer $r\times k$ right-hand side. |
| Rational positive-definite factorization of size $r+k$. | Reuse the existing factorization of $H$. |
| Up to `r` rational rank-one reductions with new matrices. | One multi-RHS solve, $O(r^2k)$ integer operations. |
| Final strict copositivity on a $k\times k$ matrix. | Same final strict copositivity on a $k\times k$ matrix. |

The retained Release benchmark used CPU 2 and one persistent Pybind process per build. Repeated measurements showed approximately
45.7% less time for matrix 19, 22.8% less for matrix 20, and 9.5% less for matrix 29. Matrix 32 was unchanged within noise because its
candidates with outside best replies fail the reduced-Hessian inertia check before Schur construction. Small matrices 4-7
improved by roughly 8-20% in the first run. The complete C++ suite, all 66 Python tests, and the exact database comparison passed
before retention.

## 12. Work excluded from the original patch

The original patch deliberately excluded:

- explicit construction of $H^{-1}$;
- a generic reusable factorization class;
- support for replaying symmetric swaps or coordinate additions;
- a second copy of $\widehat G$;
- a new integer copositivity implementation;
- a specialized immediate-integer fast path for the reused right-hand-side solve;
- changes to the existing candidate factorization loop; or
- unrelated solver, parser, Python, CLI, or support-generator changes.

These additions were unnecessary for correctness and would have made the patch harder to audit. Add one only after a measured result
shows that the smaller implementation is insufficient.
