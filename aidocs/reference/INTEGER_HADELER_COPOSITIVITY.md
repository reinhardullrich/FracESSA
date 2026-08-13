# Integer Hadeler Stability And Copositivity

Status: historical implementation and mathematical reference. Initially implemented and verified on 2026-08-06, restored as the
production exact fallback on 2026-08-08, and removed from FracESSA when coposit became the complete strict-copositivity backend on
2026-08-12. The mathematics and measurements below document that retired implementation; they do not describe current source paths.

Initial implementation and verification: 2026-08-06. Production restoration: 2026-08-08.

## Historical Implemented Boundary

At this checkpoint, the scaled reduced B matrix and every stability decision already used exact integer arithmetic. The Hadeler checker also
enumerates principal subsets in the required cardinality order and obtains each exact determinant from the retained general
fraction-free $LDL^T$ factorization.

The formerly inefficient work after a nonpositive determinant was:

- a negative determinant constructs the complete adjugate by solving for every identity column;
- a zero determinant constructs the complete adjugate from all signed cofactors.

The implementation below replaced both branches without changing the surrounding checker at that checkpoint.

## 5. Integer Hadeler decision

### 5.1 Exact scope

For every principal matrix

$$
C=M[L,L]\in\mathbb Z^{\ell\times\ell},
$$

replace the explicit-adjugate decision with:

1. $det(C)>0$: pass immediately, unchanged;
2. $det(C)<0$: solve one system $Cy=-\mathbf 1$ from the factorization already calculated for $C$;
3. $det(C)=0$: compute one exact right-nullspace basis of the restored, unfactored $C$.

The outer cardinality order, Gosper subset generator, direct formulas for dimensions one through three, public checker interface,
candidate results, and stability reasons must remain unchanged.

### 5.2 Required mathematical invariant

Principal subsets are visited by increasing cardinality. Therefore, when the checker reaches an $\ell\times\ell$ principal matrix
$C$, every proper principal submatrix of $C$ has already passed strict copositivity.

This is the hypothesis of Hadeler's Theorem 3. Under this hypothesis, the following statements are equivalent:

1. $C$ is not strictly copositive.
2. For every $b>0$, there are $x>0$ and $\lambda\leq0$ such that $Cx=\lambda b$.
3. $\det(C)\leq0$ and $\operatorname{adj}(C)>0$ entrywise.

The replacement below is valid only inside this cardinality-ordered Hadeler enumeration. It is not a general standalone rule for
an arbitrary symmetric matrix.

### 5.3 The three determinant branches

#### Positive determinant

If

$$
\det(C)>0,
$$

Hadeler's third condition is false, so $C$ passes. The code at this checkpoint already did this and remained unchanged.

#### Negative determinant: one retained solve

If $\det(C)<0$, then $C$ is nonsingular. Choose $b=\mathbf 1$ in Hadeler's second statement. A failing matrix has
$x>0$ and $\lambda<0$ with

$$
Cx=\lambda\mathbf 1.
$$

After setting $y=x/(-\lambda)$,

$$
Cy=-\mathbf 1,
\qquad y>0.
$$

Conversely, any positive solution of this system is already a direct rejection witness because

$$
y^TCy=-y^T\mathbf 1<0.
$$

Consequently, under the Hadeler invariant,

$$
C\text{ fails strict copositivity}
\quad\Longleftrightarrow\quad
C^{-1}(-\mathbf 1)>0.
$$

The existing `fraction_free_ldlt_factorization::solve_inplace()` must perform this solve. It reuses the factorization that already
provided $\det(C)$ and returns integer numerators with the positive denominator $|\det(C)|$. Therefore:

- fill one $\ell\times1$ integer matrix with `-1`;
- solve it through the retained factorization;
- reject exactly when every resulting numerator is strictly positive;
- do not construct, normalize, or inspect a complete inverse or adjugate.

Do not call `fmpz_mat_solve()`: that would repeat the elimination already performed by the retained factorization. The earlier plan's
warning about an unknown denominator sign also no longer applies because the project's solver explicitly returns a positive
denominator.

#### Zero determinant: one exact nullspace

If $\det(C)=0$, Hadeler's failing singular case has rank $\ell-1$ and a strictly positive null vector. Thus:

- nullity one with a strictly positive basis column: reject;
- nullity one with a strictly negative basis column: reject after reversing its sign conceptually;
- nullity one with mixed signs or any zero component: pass;
- nullity at least two: pass.

The last case follows because rank at most $\ell-2$ makes every $(\ell-1)$-minor zero, hence
$\operatorname{adj}(C)=0$, while Hadeler requires an entrywise strictly positive adjugate for rejection. These conclusions depend on
the proper-principal-submatrix invariant above.

Call

```cpp
fmpz_mat_nullspace(nullspace.native_handle(), principal_submatrix.native_handle());
```

after restoring `principal_submatrix` to the original principal matrix. FLINT returns the nullity and writes independent right-nullspace basis
vectors into the first `nullity` columns of an output allocated with at least $\ell\times\ell$ entries. When the nullity is one,
inspect only column zero. The determinant has already proved singularity, so retain only a debug assertion that the returned nullity
is positive.

### 5.4 Why the six-step algorithm is correct

This subsection separates Hadeler's published theorem from the consequences used by FracESSA.

Let the complete matrix being tested be

$$
M\in\mathbb{R}^{k\times k}.
$$

Strict copositivity means

$$
z^T M z>0
$$

for every nonzero vector \(z\geq0\).

#### Step 1: enumerate principal submatrices only

Let \(L\) be the support of \(z\), meaning the indices where \(z_i>0\). Remove the zero components and call the remaining strictly
positive vector \(x=z_L\). Then

$$
z^T M z=x^T M[L,L]x.
$$

The same index set \(L\) occurs on both sides of \(M\), so \(M[L,L]\) is a principal submatrix. Thus any copositivity failure of
\(M\) is already a failure of one principal submatrix. A submatrix \(M[I,J]\) with different row and column sets cannot arise from
the quadratic form and need not be enumerated.

#### Step 2: compute the determinant after smaller cardinalities have passed

Hadeler defines a matrix of order \(\ell\) as strictly copositive of order \(\ell-1\) when every principal submatrix of order
\(\ell-1\) is strictly copositive. FracESSA enumerates cardinalities \(1,2,\ldots,k\). Therefore, when it reaches

$$
C=M[L,L]\in\mathbb{R}^{\ell\times\ell},
$$

every proper principal submatrix of \(C\) has already passed. This establishes the hypothesis of Hadeler's Theorem 3.

Under that hypothesis, the theorem states that the following three claims are equivalent:

1. \(C\) is not strictly copositive.
2. For every \(b>0\), there are \(x>0\) and \(\lambda\leq0\) such that \(Cx=\lambda b\).
3. \(\det(C)\leq0\) and \(\operatorname{adj}(C)>0\) entrywise.

This equivalence is the mathematical decision rule. Computing \(\det(C)\) first separates its three possible cases.

#### Step 3: a positive determinant passes

If

$$
\det(C)>0,
$$

Hadeler's third condition is false. Therefore the first condition is false, and \(C\) is strictly copositive. This is a direct use
of the published theorem, not a new FracESSA result.

#### Step 4: a negative determinant needs one solve

If \(\det(C)<0\), then \(C\) is invertible. Choose the one fixed positive vector \(b=\mathbf{1}\).

If \(C\) is not strictly copositive, Hadeler's second statement gives \(x>0\) and \(\lambda\leq0\) with

$$
Cx=\lambda\mathbf{1}.
$$

Invertibility excludes \(\lambda=0\), because that would imply \(Cx=0\) with \(x\neq0\). Hence \(\lambda<0\). Define

$$
y=\frac{x}{-\lambda}>0.
$$

Then

$$
Cy=-\mathbf{1}.
$$

Conversely, if the unique solution of \(Cy=-\mathbf{1}\) satisfies \(y>0\), then

$$
y^TCy=-y^T\mathbf{1}=-\sum_{i=1}^{\ell}y_i<0.
$$

The vector \(y\) is a direct failure witness. Consequently, under Hadeler's hypothesis,

$$
C\text{ is not strictly copositive}
\quad\Longleftrightarrow\quad
C^{-1}(-\mathbf{1})>0.
$$

Hadeler states the quantified equation in Step 2. Choosing \(b=\mathbf{1}\) and reducing it to one solve is the FracESSA
derivation.

#### Step 5: a zero determinant needs only the nullspace

If \(\det(C)=0\), then \(C\) has a nonzero nullspace.

If that nullspace is one-dimensional and its basis vector \(v\) is strictly one-sign, choose its sign so that \(v>0\). Then

$$
Cv=0
\qquad\text{and}\qquad
v^TCv=0.
$$

Therefore \(C\) is not strictly copositive.

For the converse, Hadeler's proof shows that a singular matrix which fails strict copositivity under the proper-principal-submatrix
hypothesis must be positive semidefinite of rank \(\ell-1\), with a strictly positive null vector. Thus its nullspace must be
one-dimensional and one-sign. A one-dimensional mixed-sign basis, or a basis containing a zero component, cannot be the singular
failure described by the theorem.

If the nullity is at least two, then

$$
\operatorname{rank}(C)\leq\ell-2.
$$

Every \((\ell-1)\times(\ell-1)\) minor is then zero, so

$$
\operatorname{adj}(C)=0.
$$

Hadeler requires the adjugate to be strictly positive entrywise for failure. Therefore this case passes. Replacing all cofactors by
one nullspace calculation is a derived implementation of Hadeler's criterion, not a separate copositivity theorem.

#### Step 6: no general submatrices or explicit adjugate remain

The old implementation evaluates Hadeler's third statement literally by constructing \(\operatorname{adj}(C)\). Off-diagonal
adjugate entries are cofactors obtained by deleting different row and column indices, so that implementation creates non-principal
minors internally.

Steps 4 and 5 decide exactly the same condition through one solve or one nullspace. The outer algorithm still enumerates all
principal submatrices, but the replacement no longer constructs general submatrices, individual cofactors, a complete inverse, or
an explicit adjugate.

### 5.5 Independent literature cross-check

An internet literature search on 2026-08-06 found independent statements supporting every mathematical part of the decision. It did
not find the complete six-step implementation written in exactly this form. The distinction is:

| FracESSA step | Independent source | Status |
|---|---|---|
| Only principal submatrices matter | Hiriart-Urruty and Seeger, Proposition 3.1; Kaplan's principal-submatrix eigenvalue criterion | Published directly |
| Proper principal submatrices establish the recursive hypothesis | Cottle-Habetler-Lemke's order-\(n-1\) criterion; *Matrix Positivity*, Definition 6.3.1 | Published directly |
| \(\det(C)>0\) passes | Cottle-Habetler-Lemke strict criterion; *Matrix Positivity*, Theorem 6.3.5 | Immediate contrapositive |
| \(\det(C)<0\) is decided by \(Cy=-\mathbf 1\) | The same criterion gives \(C^{-1}<0\) for a failing nonsingular matrix | FracESSA's one-solve specialization |
| \(\det(C)=0\) fails exactly through a positive one-dimensional nullspace | Väliaho, Theorems 4.2 and 4.4; *Matrix Positivity*, Theorem 6.3.5 | Published essentially exactly |
| General minors and the explicit adjugate can be omitted | The solve and nullspace equivalents replace the published adjugate condition | FracESSA implementation consequence |

The main independent sources are:

1. R. W. Cottle, G. J. Habetler, and C. E. Lemke,
   [*On classes of copositive matrices*](https://doi.org/10.1016/0024-3795(70)90002-9), *Linear Algebra and its Applications* 3
   (1970), 295-310. This is the earlier source of the determinant/adjugate induction used by Hadeler. For a matrix strictly
   copositive of order \(n-1\), failure is characterized by \(\det(C)\leq0\) and \(\operatorname{adj}(C)>0\).
2. H. Väliaho, [*Almost copositive matrices*](https://doi.org/10.1016/0024-3795(89)90402-3), *Linear Algebra and its Applications*
   116 (1989), 121-134. Theorems 4.2 and 4.4 state that an almost strictly copositive matrix is positive semidefinite, has rank
   \(n-1\), and has a strictly positive eigenvector for eigenvalue zero. This validates the singular failure shape in Step 5.
3. W. Kaplan, [*A test for copositive matrices*](https://doi.org/10.1016/S0024-3795(00)00138-5), *Linear Algebra and its
   Applications* 313 (2000), 203-206. Kaplan proves that a symmetric matrix is strictly copositive exactly when no principal
   submatrix has a strictly positive eigenvector with a nonpositive eigenvalue. This independently confirms principal-only
   enumeration and positive-nullvector rejection.
4. J.-B. Hiriart-Urruty and A. Seeger,
   [*A variational approach to copositive matrices*](https://www.math.univ-toulouse.fr/~mongeau/JBHU-copositive.pdf), *SIAM
   Review* 52 (2010), 593-629. Section 3.2 surveys the Cottle-Habetler-Lemke determinant/adjugate test. Section 4.1 derives the
   principal-submatrix eigenvalue view and explains why all \(2^n-1\) nonempty principal subsets appear.
5. C. R. Johnson, R. L. Smith, and M. J. Tsatsomeros,
   [*Matrix Positivity*, Chapter 6](https://doi.org/10.1017/9781108778619.007), Cambridge University Press, 2020. Theorem 6.3.5
   restates the strict Cottle-Habetler-Lemke criterion and separates the failing cases into negative determinant, with one negative
   eigenvalue and a positive eigenvector, and zero determinant, with rank \(n-1\) and a positive null eigenvector. Theorem 6.3.11
   restates Kaplan's strict principal-submatrix criterion.

The one-solve reduction is shorter than the published full-inverse or adjugate tests but is not an additional assumption. In the
negative-determinant case, the published criterion gives

$$
C^{-1}=\frac{\operatorname{adj}(C)}{\det(C)}<0
$$

for every failing matrix. Hence \(-C^{-1}\mathbf 1>0\), which is exactly the solution of \(Cy=-\mathbf 1\). Conversely, any
positive solution is already a direct quadratic-form witness. No source found in this search states this fixed-right-hand-side
optimization as an algorithmic recommendation; it is FracESSA's algebraic specialization of the classical theorem.

Likewise, the nullity-at-least-two pass is not a new copositivity theorem. It is the contrapositive of the published singular
classification: every singular failure has rank \(n-1\), so a matrix of nullity at least two cannot be such a failure once all
proper principal submatrices have passed.

### 5.6 Historical core implementation

The then-current `is_strictly_copositive_hadeler()` generic branch:

1. extracts `subset_indices` into the fixed stack buffer;
2. copies the current principal matrix into `principal_submatrix`;
3. factorizes `principal_submatrix` in place;
4. reads the exact determinant sign.

The determinant decision after step 4 is:

```cpp
if (determinant_sign > 0) return true;

if (nonsingular) { // Here determinant_sign < 0.
    matrix_int solution(current_dim, 1);
    const integer minus_one(-1);
    for (size_t row = 0; row < current_dim; ++row) solution(row, 0) = minus_one;

    integer denominator;
    factorization_.solve_inplace(solution, denominator, principal_submatrix);
    assert(denominator.sign() > 0);

    for (size_t row = 0; row < current_dim; ++row)
        if (solution(row, 0).sign() <= 0) return true;
    return false;
}

// factorize_inplace() overwrote principal_submatrix. Restore the original principal matrix only for this singular branch.
for (size_t row = 0; row < current_dim; ++row)
    for (size_t column = 0; column < current_dim; ++column)
        principal_submatrix(row, column) = A(subset_indices[row], subset_indices[column]);

matrix_int nullspace(current_dim, current_dim);
const slong nullity = fmpz_mat_nullspace(nullspace.native_handle(), principal_submatrix.native_handle());
assert(nullity > 0);
if (nullity != 1) return true;

const int basis_sign = nullspace(0, 0).sign();
if (basis_sign == 0) return true;
for (size_t row = 1; row < current_dim; ++row)
    if (nullspace(row, 0).sign() != basis_sign) return true;
return false;
```

The implementation deliberately has no helper class, generic witness abstraction, solver wrapper, or persistent scratch framework.
The local objects are already smaller than the adjugate matrices they replaced. Reusable allocation remains unjustified unless a
representative benchmark shows that these allocations matter.

### 5.7 Historical one-solve patch scope

The original 2026-08-06 replacement of explicit adjugates changed only two source files:

| File | Minimal change |
|---|---|
| `cpp/include/linalg/copositive_integer.hpp` | Include `<cassert>` directly, replace the negative and singular adjugate branches, and delete the two now-unused private helpers. |
| `cpp/tests/test_copositivity.cpp` | Rename the affected branch tests and add the missing passing witnesses. |

That isolated replacement required no change in:

- `fraction_free_ldlt.hpp`: it already solves one or many right-hand sides and guarantees a positive denominator;
- `matrix_integer.hpp`: it already exposes the FLINT matrix handle and required entry operations;
- CMake or the FLINT minimum: `fmpz_mat_nullspace` exists in both supported local FLINT 3.4 and project-local FLINT 3.6;
- `check_stability.cpp`, candidate search, support generation, CLI, Pybind, Python, SQLite, or release workflow.

The 2026-08-08 production restoration necessarily moved the already-verified general factorization from `archive/` back to active
`linalg`, reconnected `check_stability.cpp`, and retained the independently implemented negative-component reduction. Those restoration
changes do not alter the one-solve/nullspace proof.

### 5.8 Focused regression cases

The focused matrices use dimension four so they reach the generic Hadeler branch after their proper principal submatrices pass.

1. **Negative determinant, reject:** diagonal `5`, off-diagonal `-2`. The all-ones vector solves $Cy=-\mathbf 1$ and is positive.
2. **Negative determinant, pass:** diagonal `1`, off-diagonal `2`. The matrix is strictly copositive, while the unique solution of
   $Cy=-\mathbf 1$ is strictly negative.
3. **Nullity one, reject:** diagonal `3`, off-diagonal `-1`. The nullspace is spanned by the positive all-ones vector.
4. **Nullity one, pass:** $4I-vv^T$ for $v=(1,-1,1,-1)^T$. Its only null direction has mixed signs.
5. **Nullity three, pass:** the all-ones matrix. Its nullspace has dimension three, but $x^TCx=(\mathbf 1^Tx)^2>0$ for every
   nonzero $x\geq0$.
6. **Arbitrary precision:** multiply at least one negative-determinant and one singular case by a positive integer larger than 64
   bits; the decisions must be unchanged.

The suite also retains the one-, two-, and three-dimensional, sign-scan, and early-termination tests.

### 5.9 Historical verification procedure

The original one-solve/nullspace change used this acceptance procedure:

1. Build the current commit in a separate ignored Release directory and preserve it as the before-change checker.
2. Run the before and after checkers over the then-retained 1,069-row copositivity corpus. All permutation-inequivalent matrices had
   to retain their stored strict-copositivity result: 41 positive and 1,028 negative.
3. Compare complete safe candidate output before and after for source matrix IDs 4, 7, 19, 20, 29, 53, and 54. Exclude only timing;
   all candidate fields and order must match exactly.
4. Run the complete CTest and Python test suites against project-local FLINT 3.6.
5. Configure and run the focused C++ tests against system FLINT 3.4, proving that this task does not accidentally change the
   dependency floor.

Any mismatch blocks the change. No tolerance comparison is relevant because every operation and expected result is exact.

The corpus was deliberately not wired into production or CTest. One temporary ignored C++ driver performed Step 2 and the corpus
benchmark and was then removed; no permanent framework was added for this one implementation comparison.

### 5.10 Historical performance procedure

After correctness passed, the original change compared clean before and after Release/native/LTO builds with the same compiler,
FLINT 3.6, persistent process, CPU 2, and native nanosecond medians.

Measure:

1. the complete 1,069-row reduced-B corpus;
2. the negative-determinant and singular dimension-four cases separately, repeated long enough to remove startup noise;
3. the seven end-to-end source games listed above, to ensure the ordinary application path does not regress.

The expected gain is specific:

- negative determinant: one triangular solve with one column instead of solves for $\ell$ identity columns;
- zero determinant: one nullspace elimination instead of about $\ell^2/2$ cofactor determinants of order $\ell-1$.

Retain the change only if all exact results match, the affected branches improve repeatably, and the ordinary end-to-end path has no
material regression. Do not add factorization caching or reusable allocation merely because it looks faster; benchmark this minimal
replacement first.

### 5.11 Historical scope exclusions

The original one-solve/nullspace patch intentionally excluded:

- negative-edge connected-component decomposition;
- recursive row, duplicate-row, or Schur reductions;
- changed principal-subset generation or ordering;
- new early acceptance or rejection rules;
- reason-string or database migrations;
- rational-helper cleanup outside the two deleted adjugate methods;
- a FLINT version gate;
- a generic linear-algebra or copositivity interface;
- changes to the KKT candidate factorization.

Those were independent decisions. The negative-entry connected-component decomposition was implemented and verified later, before
the complete checker moved to coposit. The mathematical job described here remains unchanged: replace explicit adjugate construction
with one retained solve or one exact nullspace while preserving Hadeler's exact decision.

### 5.12 References

- K. P. Hadeler, *On copositive matrices*, Linear Algebra and its Applications 49 (1983), Theorem 3; audited local transcription:
  `research/papers/Hadeler_1983.md`.
- Detailed local derivation of the fixed-right-hand-side test: `research/HADELER_ONE_SOLVE_REPLACEMENT.md`.
- FLINT 3.6 integer-matrix documentation: <https://flintlib.org/doc/fmpz_mat.html>. `fmpz_mat_nullspace` returns the exact right
  nullity and writes basis vectors into columns; the official 3.6 source writes them into the first `nullity` columns.

### 5.13 Verified 2026-08-06 result

- Only `copositive_integer.hpp` and its focused test file changed in executable source.
- All 1,069 reduced-B corpus rows retained their stored exact result: 41 strictly copositive and 1,028 not strictly copositive.
- Complete safe candidate output for source IDs 4, 7, 19, 20, 29, 53, and 54 remained byte-identical after removing only timing.
- All 10 FLINT 3.6 CTest targets, all 66 Python tests, and the focused checker tests against system FLINT 3.4 passed.
- Nine alternating CPU-2 corpus repetitions reduced the uncontended wall-time median from 4.660 s to 0.995 s, a 4.68x speedup
  and 78.65% less time. Process CPU-time medians independently gave the same 4.68x result; scheduler delay stayed below 0.93%.
- CPU-2 end-to-end medians for the larger IDs 29, 53, and 54 changed by -0.08%, -0.09%, and -0.19%, respectively; no material
  ordinary-path regression was measured.
