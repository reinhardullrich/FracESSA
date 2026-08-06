# Integer Stability Matrix Migration Only

Status: proposed for review; not approved and not implemented.

Date: 2026-08-06.

## 1. Exact scope

Change only the final exact stability matrix from denominator-one rational storage (`fmpq`) to integer storage (`fmpz`). Keep the
current mathematical algorithm and the current decision order.

The affected value is the scaled reduced B matrix built after an exact candidate has already been found. Its entries are already
calculated as integers. The current code needlessly converts every entry to a fraction with denominator one before checking
positive definiteness and strict copositivity.

The intended data flow is:

```text
current:
integer scaled entry -> denominator-one fraction matrix -> rational stability checks

proposed:
integer scaled entry -> integer matrix -> integer stability checks
```

This is not a project-wide removal of rational arithmetic. Rational arithmetic remains where exact fractions are genuine values:

- parsed game entries;
- candidate probabilities;
- candidate payoffs; and
- public candidate output.

## 2. Explicit non-goals

This migration must not introduce any of the algorithm changes proposed in
`INTEGER_STABILITY_COPOSITIVITY.md`.

In particular, do not:

- replace the complete inverse by one solve with the right-hand side `-1`;
- replace the singular adjugate by a nullspace calculation;
- introduce Danninger's recursive algorithm;
- introduce connected-component or row-removal reductions;
- change Hadeler's principal-subset enumeration;
- change cardinality order or Gosper generation;
- add reusable solver frameworks or generic numeric templates;
- change any stability reason string; or
- combine this migration with unrelated cleanup or optimization.

## 3. Behaviour that must remain identical

The following path in `fracessa::check_stability()` remains in exactly the same order:

1. decide the one-, two-, and three-dimensional cases with the existing formulas;
2. perform the existing shared sign scan;
3. reject a nonpositive diagonal;
4. accept nonnegative off-diagonal entries;
5. accept negative-part diagonal dominance;
6. reject a nonpositive all-ones quadratic value;
7. accept a positive-definite matrix;
8. reject a non-positive-definite Z-matrix; and
9. run the current Hadeler determinant-and-adjugate test.

Candidate enumeration, candidate IDs, support order, ESS decisions, candidate rows, and stability reason strings must be unchanged.
The historically named reasons `T_pd_frc` and `F_not_copos` remain unchanged in this task even though `frc` will no longer describe
the internal stability-matrix representation. Renaming them would create an unrelated output and database migration.

## 4. Mechanical arithmetic translation

### 4.1 Scaled reduced B construction

`find_candidate_safe::build_scaled_reduced_b()` already finishes each entry in `linalg::integer scaled_entry`.

Replace only this final conversion:

```cpp
result(row, column).set_ratio(scaled_entry, one);
```

with assignment to an `linalg::matrix_int` entry. Mirror the other triangle exactly as today. Do not change the Schur-complement
formula, scaling, loop order, retained factorization, or returned `kay_size`.

### 4.2 Small dimensions and the sign scan

Translate the existing one-, two-, and three-dimensional formulas operation for operation from `fraction` to `integer`. Every
product and comparison remains exact. The shared sign scan reads integer entries directly instead of asserting that an `fmpq`
denominator equals one and then extracting its numerator.

No new formulas or decisions are introduced.

### 4.3 Positive-definiteness check

Use FLINT's exact `fmpz_mat_is_spd()` predicate on the integer matrix. It decides the same mathematical question as the current
rational `LDL^T` method: whether the complete symmetric matrix is positive definite.

This changes the internal factorization implementation, because an integer matrix should use FLINT's integer predicate instead of
being converted back to rational form. It does not add, remove, or reorder a stability decision. FLINT certifies the result exactly;
an inconclusive numerical attempt falls back to an exact method internally.

Writing a second custom fraction-free positive-definiteness factorization merely to imitate the current hand-written rational
`LDL^T` steps is deliberately excluded. It would add substantial code without preserving any externally visible behaviour.

### 4.4 Hadeler determinant and adjugate

Retain the current Hadeler implementation literally at the mathematical level.

For every principal subset of cardinality at least four:

1. copy the same principal submatrix into an integer matrix;
2. calculate its exact determinant with `fmpz_mat_det()`;
3. pass immediately when the determinant is positive;
4. when the determinant is negative, calculate the complete inverse with `fmpz_mat_inv()`;
5. construct the complete exact adjugate from

   $$
   \operatorname{adj}(A)=\det(A)A^{-1};
   $$

6. when the determinant is zero, construct the adjugate from all the same signed cofactor determinants as today; and
7. reject exactly when every adjugate entry is strictly positive.

`fmpz_mat_inv()` returns an integer numerator matrix `N` and an integer denominator `q` with `A^{-1}=N/q`. FLINT guarantees that
`q` divides `det(A)`, so the integer adjugate is

$$
\operatorname{adj}(A)=\frac{\det(A)}{q}N.
$$

Use exact integer division. Do not replace this complete inverse with the cheaper one-solve method in this task.

For a singular matrix, keep the current expensive cofactor construction. Each minor determinant changes from rational LU to
`fmpz_mat_det()`, but the set of minors and the adjugate test remain identical. The proposed nullspace replacement is explicitly
deferred.

## 5. Minimal source-file scope

The implementation should change only these source files:

1. `cpp/include/linalg/matrix_integer.hpp`
   - add exact positive-definiteness and log-string support needed by the existing caller;
2. `cpp/include/linalg/copositive_fraction.hpp` -> `cpp/include/linalg/copositive_integer.hpp`
   - replace the rational checker with the operation-for-operation integer version;
3. `cpp/include/fracessa/find_candidate_safe.hpp`
   - change the scaled reduced B output type to `matrix_int`;
4. `cpp/src/find_candidate_safe.cpp`
   - store the already-calculated integer entries directly;
5. `cpp/include/fracessa/fracessa.hpp`
   - store `scaled_reduced_b_` as `matrix_int` and correct the directly affected comments;
6. `cpp/src/checkstab.cpp`
   - include and call the integer checker without changing its control flow; and
7. `cpp/tests/test_copositivity.cpp`
   - translate checker tests to integer matrices and add coverage for the genuine dimension-four Hadeler branches.

No Python file, parser, CLI interface, Pybind interface, database schema, support generator, or candidate solver should change.
No CMake change is required for the arithmetic APIs in this plan: `fmpz_mat_det`, `fmpz_mat_inv`, and `fmpz_mat_is_spd` already
exist in the currently installed FLINT 3.4 as well as the project-local FLINT 3.6.

The old rational `copositive_fraction.hpp` should not remain beside the integer checker after the switch. Two production checkers
would create drift. The separate rational LU file and its direct unit test are outside this migration and should remain untouched
until a later, separately approved unused-code cleanup.

## 6. Verification plan

### 6.1 Preserve a before-change reference

Before editing source, build the current rational version in an ignored local build directory and retain its executable and Python
extension. Do not create another worktree and do not copy production source.

### 6.2 Focused unit tests

Port the current copositivity tests to integer inputs. Add true dimension-four cases because the existing tests named after Hadeler
mostly stop in the special two-dimensional formulas:

- positive determinant: a positive-definite integer matrix;
- negative determinant and nonsingular full-inverse branch: diagonal `5`, every off-diagonal `-2`;
- zero determinant and singular cofactor branch: diagonal `3`, every off-diagonal `-1`;
- arbitrary-precision integer entries beyond machine-word range; and
- the existing sign-scan and low-dimensional boundary cases.

The two equicorrelation examples have all proper principal submatrices strictly copositive, so they reach the intended full-size
Hadeler branches rather than failing earlier.

### 6.3 Complete project checks

Run the clean Release build, all CTest tests, and the maintained Python test suite with the project-local FLINT 3.6 installation.
Also configure and run the focused C++ tests against system FLINT 3.4 to verify that this limited migration did not accidentally
raise the dependency floor.

### 6.4 Rational-versus-integer equivalence

Run the retained before-change build and the new build over the same canonical SQLite matrices in safe mode. Compare:

- candidate count and ESS count;
- candidate and ESS support-size structures;
- every candidate row;
- support, extended support, multiplier, ESS flag, and stability reason; and
- exact probability and payoff strings when candidates are materialized.

The comparison must be exact, not tolerance-based. Timing fields are excluded. Any mismatch blocks the migration.

### 6.5 Performance gate

Benchmark the old rational and new integer builds only after correctness passes. Use matrices that actually reach the scaled reduced
B and full copositivity paths, covering small and large, circular and non-circular cases. Retain the migration only if it does not
cause a material regression; make no speed claim before measurement.

## 7. Implementation order and review gates

1. Build and retain the current rational reference.
2. Make only the seven-file source/test change listed above.
3. Show the complete diff for direct review before any cleanup.
4. Run focused unit tests and the complete local suites.
5. Run exact rational-versus-integer matrix equivalence.
6. Run the representative performance comparison.
7. Update current project knowledge and append the result to `CHANGES.md` only after the implementation is accepted.
8. Consider removing the now-production-unused rational LU and its tests only as a separate later task.

No commit, push, release, algorithm replacement, or unrelated cleanup belongs to this migration.

## 8. Approval requested later

Implementation approval would cover exactly the seven source/test paths in Section 5 and the behaviour described above. If the
implementation reveals that another source file is necessary, stop and request approval before changing it.
