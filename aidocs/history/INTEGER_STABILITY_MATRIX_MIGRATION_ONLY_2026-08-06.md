# Integer Stability Matrix Migration Only

Status: implemented and verified on 2026-08-06.

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
`../reference/INTEGER_HADELER_COPOSITIVITY.md`.

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

Use the retained `fraction_free_ldlt_factorization` on an integer copy. Each Bareiss pivot divided by the previous pivot is the
corresponding ordinary $LDL^T$ diagonal entry, so the factorization records exact inertia while it runs. The matrix is positive
definite exactly when it is nonsingular and every recorded diagonal sign is positive.

This reuses the general integer factorization requested after the original plan was written. It does not add, remove, or reorder a
stability decision and does not depend on floating-point intervals.

### 4.4 Hadeler determinant and adjugate

Retain the current Hadeler implementation literally at the mathematical level.

For every principal subset of cardinality at least four:

1. copy the same principal submatrix into an integer matrix;
2. factor it with `fraction_free_ldlt_factorization` and read its exact determinant;
3. pass immediately when the determinant is positive;
4. when the determinant is negative, solve all identity columns from the retained factorization to obtain
   $|\det(A)|A^{-1}$, then negate it to obtain the complete adjugate;
5. retain the identity

   $$
   \operatorname{adj}(A)=\det(A)A^{-1};
   $$

6. when the determinant is zero, construct the adjugate from all the same signed cofactor determinants as before; and
7. reject exactly when every adjugate entry is strictly positive.

The retained solve uses the positive denominator $|\det(A)|$. The negative-determinant branch therefore returns
$- |\det(A)|A^{-1}=\det(A)A^{-1}=\operatorname{adj}(A)$. This is still the complete inverse-equivalent calculation; it is not the
deferred one-right-hand-side optimization.

For a singular matrix, keep the current expensive cofactor construction. Off-diagonal cofactor minors are not symmetric and cannot
use $LDL^T$, so their determinants use `fmpz_mat_det()`. The set of minors and the adjugate test remain identical. The proposed
nullspace replacement is explicitly deferred.

## 5. Minimal source-file scope

The implementation should change only these source files:

1. `cpp/include/linalg/matrix_integer.hpp`
   - add identity and log-string support needed by the integer checker;
2. `cpp/include/linalg/copositive_fraction.hpp` -> `cpp/include/linalg/copositive_integer.hpp`
   - replace the rational checker with the operation-for-operation integer version;
3. `cpp/include/linalg/fraction_free_ldlt.hpp`
   - retain exact inertia and expose its positive-definiteness decision;
4. `cpp/include/fracessa/find_candidate_safe.hpp`
   - change the scaled reduced B output type to `matrix_int`;
5. `cpp/src/find_candidate_safe.cpp`
   - store the already-calculated integer entries directly;
6. `cpp/include/fracessa/fracessa.hpp`
   - store `scaled_reduced_b_` as `matrix_int` and correct the directly affected comments;
7. `cpp/src/checkstab.cpp`
   - include and call the integer checker without changing its control flow; and
8. `cpp/tests/test_copositivity.cpp`
   - translate checker tests to integer matrices and add coverage for the genuine dimension-four Hadeler branches;
9. `cpp/tests/test_linear_solver.cpp`
   - cross-check retained inertia against FLINT positive definiteness.

No Python file, parser, CLI interface, Pybind interface, database schema, support generator, or candidate solver should change.
No CMake change is required. The only direct FLINT matrix operation newly used by the checker is `fmpz_mat_det`, which exists in
both the installed FLINT 3.4 and project-local FLINT 3.6.

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

## 7. Verified result

- The focused low-dimensional, sign-scan, positive-determinant, negative-determinant, singular-cofactor, and arbitrary-precision
  tests pass.
- All ten C++/CLI tests pass with project-local FLINT 3.6 and system FLINT 3.4.
- All 66 Python tests and 49 subtests pass with the new native module.
- Old and new safe-mode candidate output matched byte-for-byte for 883 stored games after removing only `elapsed_ns`. Matrix 2998
  was skipped after the old executable exceeded a 30-second guard despite its stale calibration value.
- Nine pinned repetitions on CPU 2 gave new/old median runtime ratios of 0.292, 0.986, 1.003, and 1.006 for matrix IDs 137, 29, 53,
  and 54 respectively. The migration shows no material regression and a large gain on the copositivity-heavy case.

The rational LU implementation and its self-focused tests remain for a later separately approved cleanup. No algorithm replacement,
database migration, reason rename, commit, push, or release is part of this migration.
