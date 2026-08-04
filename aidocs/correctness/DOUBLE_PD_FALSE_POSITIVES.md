# Why Double Positive-Definiteness Was Not a Valid ESS Certificate

Last verified: 2026-07-26

## Summary

The old stability path treated a successful binary64 Cholesky calculation as a
final proof that an exact rational Bee matrix was positive definite. That was
mathematically unsound.

A floating-point Cholesky result can be useful as a heuristic, but it cannot
certify the sign of an exact matrix near the positive-definite boundary. Two
different rational matrices can round to exactly the same binary64 matrix even
when one is positive definite and the other is indefinite. No tolerance applied
after that information has been discarded can distinguish them.

Verification matrices 36 and 37 exposed the concrete failure:

- each game has zero ESS in exact arithmetic;
- the old double shortcut reported two ESS;
- all four incorrect candidates were classified as `T_pd_dbl`;
- removing the shortcut makes both `fast` and `safe` return zero ESS.

The defect was especially dangerous because only the positive floating result
was final. A negative floating result continued to exact arithmetic, but a
positive floating result skipped the exact check and changed the answer.

## How FracESSA Constructs Bee

Let `A` be the symmetric payoff matrix and let strategy 0 be the reference
strategy `m`. For reduced indices `i,j`, FracESSA constructs

```text
B[i,j] =
    A[m,j] + A[j,m] + A[i,m] + A[m,i]
  - A[i,j] - A[j,i] - 2*A[m,m].
```

For the examples below, the first row and column of `A` are zero. Therefore,
for the pure candidate `x = (1,0,0,0)`,

```text
B = -2 * A[1:4, 1:4].
```

All outside strategies tie the pure candidate at payoff zero, so the extended
support contains all four strategies. Stability therefore depends on this
`3 x 3` Bee matrix.

## Regression Matrix 36

The complete symmetric game is

```text
      0       0       0        0
      0     -53    33/2     73/2
      0    33/2   -85/2       26
      0    73/2      26   -125/2
```

CLI input:

```text
4#0,0,0,0,-53,33/2,73/2,-85/2,26,-125/2
```

For the pure candidate, Bee is

```text
B36,pure =
    106   -33   -73
    -33    85   -52
    -73   -52   125
```

Its exact leading principal minors are

```text
Delta1 = 106
Delta2 = 106*85 - 33^2 = 7921
Delta3 = det(B36,pure) = 0.
```

The vector `(1,1,1)^T` is an exact null vector. The matrix is singular and is
not positive definite. It also fails strict copositivity because
`(1,1,1) B (1,1,1)^T = 0`.

The second candidate has support `{1,2,3}`. Its Bee matrix is

```text
B36,mixed =
    106   139   179
    139   257   160
    179   160   377
```

It has the same first two exact pivots and determinant zero. Its null vector is
`(-3,1,1)^T`, so it is also not positive definite.

## Regression Matrix 37

Matrix 37 is the smaller-integer version:

```text
      0      0      0      0
      0   -7/2      2    3/2
      0      2   -9/2    5/2
      0    3/2    5/2     -4
```

CLI input:

```text
4#0,0,0,0,-7/2,2,3/2,-9/2,5/2,-4
```

For the pure candidate,

```text
B37,pure =
      7   -4   -3
     -4    9   -5
     -3   -5    8
```

The exact leading principal minors are

```text
Delta1 = 7
Delta2 = 7*9 - 4^2 = 47
Delta3 = det(B37,pure) = 0.
```

Again, `(1,1,1)^T` is an exact null vector. Exact LDL^T therefore has diagonal
pivots

```text
D = (7, 47/7, 0),
```

and must reject the matrix.

The second candidate produces

```text
B37,mixed =
      7   11   10
     11   24    9
     10    9   21
```

Its leading minors are also `7`, `47`, and `0`, and
`(-3,1,1)^T` is a null vector.

## What Binary64 Computed

The old code formed a Cholesky factor in `double` and rejected a diagonal
only when

```text
computed_diagonal <= epsilon_binary64 * ||B||_infinity,
```

where `epsilon_binary64 = std::numeric_limits<double>::epsilon() = 2^-52`.

On the optimized GCC 16.1.1 build used for the regression, the exact final
diagonal was zero in all four cases, but rounding produced:

| Game | Candidate | `||B||_inf` | Rejection threshold | Computed final diagonal | Old result |
|---|---|---:|---:|---:|---|
| 36 | pure | 250 | 5.551115123125783e-14 | 5.684341886080802e-14 | accepted |
| 36 | mixed | 716 | 1.589839371263224e-13 | 1.705302565824240e-13 | accepted |
| 37 | pure | 18 | 3.996802888650564e-15 | 5.329070518200751e-15 | accepted |
| 37 | mixed | 44 | 9.769962616701378e-15 | 1.065814103640150e-14 | accepted |

Every rounded residual was slightly larger than the tolerance. The old method
therefore returned `true` four times even though every exact determinant was
zero.

The rounded residual depends on compiler, optimization, instruction selection,
and evaluation order. Another build may round one of these examples down to
zero. That variability is not a defense of the method; it is another reason a
floating result cannot be a platform-independent mathematical certificate.

## A Normal-Looking All-Integer Counterexample

Matrix 37 already has ordinary-sized values, but its game contains halves. An
even cleaner example is obtained by multiplying the entire game by four:

```text
       0    0     0     0
       0  -14     8     6
Aint = 0    8   -18    10
       0    6    10   -16
```

CLI input:

```text
4#0,0,0,0,-14,8,6,-18,10,-16
```

There is no small input value. Every payoff is an integer between `-18` and
`10`. Nevertheless, the pure candidate produces

```text
Bint = 4*B37,pure =
     28   -16   -12
    -16    36   -20
    -12   -20    32
```

Its exact leading principal minors are

```text
Delta1 = 28
Delta2 = 28*36 - 16^2 = 752
Delta3 = det(Bint) = 0.
```

The rows sum to zero, so `(1,1,1)^T` is an exact null vector. The matrix is
singular even though none of its entries looks numerically suspicious. This is
the central lesson: ordinary-sized entries do not imply a well-conditioned
matrix. Exact algebraic cancellation can place a normal-looking matrix exactly
on the positive-definite boundary.

With the old optimized double calculation:

```text
||Bint||_infinity        = 72
exact final pivot        = 0
computed final pivot     = 2.1316282072803006e-14
old code threshold       = 1.5987211554602254e-14
10^-20 relative threshold = 7.2e-19
```

The computed pivot is larger than both thresholds. The old program therefore
returned `T_pd_dbl` instead of continuing to exact arithmetic.

Observed complete-game result:

- pre-fix program: 2 ESS, both incorrectly classified `T_pd_dbl`;
- corrected exact-PD program: 0 ESS.

The factor four is deliberate. Cholesky square roots then scale by exactly two
in binary arithmetic, preserving the rounding pattern from matrix 37. More
generally, `4^k*A37` gives the same type of all-integer counterexample while the
values remain in binary64's normal finite range.

Thus a tolerance of `10^-20` does not require an input containing `10^-20`.
This small-integer game already defeats it. The tiny number is created by
rounding and cancellation inside the factorization, not supplied in the matrix.

## Normal Inputs with a Rounded Pivot Around 10^-20

There is one binary64 limitation to state first. If the final Cholesky
subtraction operates near a diagonal value of order one, adjacent doubles are
about `10^-16` apart. A nonzero stored result near that scale cannot be
`10^-20`; it is either zero or at least roughly one local binary64 spacing.

To obtain a rounded pivot around `10^-20`, the Bee matrix scale must therefore
be around `10^-5` to `10^-4`. The input payoff matrix need not be small because
adding one constant to every payoff leaves Bee unchanged.

Use the exact power-of-four scale

```text
s = 4^-9 = 2^-18 = 1/262144
```

and multiply matrix 37 by `s`:

```text
          0           0           0           0
          0   -7/524288    1/131072    3/524288
A20 =     0    1/131072   -9/524288    5/524288
          0    3/524288    5/524288     -1/65536
```

CLI input:

```text
4#0,0,0,0,-7/524288,1/131072,3/524288,-9/524288,5/524288,-1/65536
```

The nonzero payoff magnitudes range only from

```text
5.7220458984375e-6 to 1.71661376953125e-5.
```

All denominators are powers of two, so every input is represented exactly in
binary64. There is no rational-to-double conversion error in this example.

To make every displayed payoff look completely ordinary, add `1` to every
entry of `A20`:

```text
             1                 1                 1                 1
             1      524281/524288      131073/131072      524291/524288
Anear1 =     1      131073/131072      524279/524288      524293/524288
             1      524291/524288      524293/524288        65535/65536
```

CLI input:

```text
4#1,1,1,1,524281/524288,131073/131072,524291/524288,524279/524288,524293/524288,65535/65536
```

Every payoff now lies between

```text
0.99998283386230469 and 1.0000095367431641.
```

The common `1` cancels from the Bee construction. These dyadic inputs are also
exactly representable in binary64, and the optimized old build produces the
same `B20` bit for bit.

For the pure candidate,

```text
B20 = s*B37,pure =
     2.6702880859375e-5   -1.52587890625e-5   -1.1444091796875e-5
    -1.52587890625e-5      3.4332275390625e-5   -1.9073486328125e-5
    -1.1444091796875e-5   -1.9073486328125e-5    3.0517578125e-5
```

Exact arithmetic still gives

```text
Delta1 = 7*s > 0
Delta2 = 47*s^2 > 0
Delta3 = 0,
```

and `(1,1,1)^T` remains an exact null vector. Positive scaling cannot change
the exact ESS answer.

The optimized old double path produced

```text
||B20||_infinity          = 6.866455078125e-5
exact final pivot         = 0
computed final pivot      = 2.0328790734103208e-20
old code threshold        = 1.5246593050577406e-20
absolute 10^-20 threshold = 1.0e-20
relative 10^-20 threshold = 6.866455078125e-25
```

The rounded pivot is larger than all three thresholds, so each version accepts
the false positive. Observed complete-game result:

- pre-fix program: 2 false ESS;
- corrected exact-PD program: 0 ESS.

`Anear1` is the requested normal-input example with a rounded pivot around
`10^-20`: all visible inputs are approximately one, while the exact final pivot
is zero and the computed one is `2.0328790734103208e-20`. The tiny value is
generated entirely inside Cholesky by cancellation and rounding.

## An Indefinite Counterexample for Arbitrarily Small Eta

Start with the smaller pure Bee matrix and decrease only its last diagonal:

```text
             7   -4      -3
B(eta) =    -4    9      -5
            -3   -5   8-eta
```

For every exact rational `eta > 0`,

```text
Delta1 = 7
Delta2 = 47
Delta3 = det(B(eta)) = -47*eta < 0.
```

Thus `B(eta)` is indefinite for every positive eta, however small. The
positive vector `q = (1,1,1)^T` gives the even more direct witness

```text
q^T B(eta) q = -eta < 0.
```

The corresponding FracESSA game is obtained by setting the first row and column
to zero and using `A[1:4,1:4] = -B(eta)/2`:

```text
          0      0      0          0
          0   -7/2      2        3/2
A(eta) =  0      2   -9/2        5/2
          0    3/2    5/2   -4 + eta/2
```

The pure strategy `(1,0,0,0)` remains a candidate with payoff zero and all
four strategies in its extended support. Exact arithmetic rejects it because
its Bee matrix is indefinite.

### Why the perturbation disappears in double

Near `8`, the binary64 value immediately below `8` is `8 - 2^-50`.
Consequently, if

```text
0 < eta < 2^-51 = 4.440892098500626e-16,
```

round-to-nearest converts `8-eta` back to exactly `8`. Equivalently, the
payoff `-4 + eta/2` converts back to exactly `-4`.

Therefore all such exact matrices `B(eta)` become the same binary64 matrix
`B37,pure`. The old code then reproduced its positive rounding residual and
declared the pure strategy ESS.

This interval contains arbitrarily small positive rationals. Making eta smaller
does not make the floating decision safer; it makes the exact distinction more
likely to disappear during conversion.

## Concrete Eta = 10^-20 Example

Choose

```text
eta = 10^-20 = 1/100000000000000000000.
```

Then

```text
-4 + eta/2 =
-799999999999999999999 / 200000000000000000000.
```

The complete CLI input is

```text
4#0,0,0,0,-7/2,2,3/2,-9/2,5/2,-799999999999999999999/200000000000000000000
```

Exact facts:

```text
det(B(10^-20)) = -47/100000000000000000000 < 0
(1,1,1) B(10^-20) (1,1,1)^T = -10^-20.
```

Observed behavior:

- pre-fix double shortcut: the pure candidate was incorrectly `T_pd_dbl`;
- exact stability path: the pure candidate is correctly `F_not_copos`;
- the game has one separate, genuine mixed ESS, so the old total was 2 and the
  corrected total is 1.

If the Cholesky rejection threshold itself were changed to `10^-20`, the old
method would still accept the false pure ESS:

- absolute threshold: `5.329070518200751e-15 > 10^-20`;
- relative threshold `10^-20 * ||B||_inf`:
  `5.329070518200751e-15 > 1.8e-19`.

Lowering the tolerance therefore makes this false positive easier to accept.

## Scaled Construction

The construction can also defeat any chosen absolute threshold while keeping
an arbitrarily small exact negative direction.

Choose a power-of-four scale `S = 4^k` and define

```text
B(S,eta) = S*B37,pure - eta*E33,
```

where `E33` has one in the final diagonal position and zero elsewhere. Then

```text
det(B(S,eta)) = -47*S^2*eta < 0.
```

Choose `k` large enough that

```text
eta < S*2^-51
```

and the perturbation is lost when converted to binary64. Double sees exactly
`S*B37,pure`. Its spurious Cholesky residual scales with `S`, so `S` can
also make that residual larger than any fixed absolute rejection threshold,
subject only to the normal finite range of binary64.

For a relative rule `tau*||B||_inf`, the unscaled example already fails for

```text
tau < 5.329070518200751e-15 / 18
    = 2.960594732333751e-16.
```

In particular, `tau = 10^-20` is far inside the failing range.

## Why Tolerance Tuning Cannot Repair This Boundary Decision

For any sufficiently small `eta > 0`, consider the pair

```text
B_plus  = B37,pure + eta*E33
B_minus = B37,pure - eta*E33.
```

Their exact determinants are

```text
det(B_plus)  =  47*eta > 0
det(B_minus) = -47*eta < 0.
```

By Sylvester's criterion, `B_plus` is positive definite and `B_minus` is
indefinite. Yet both matrices round to the identical binary64 matrix
`B37,pure` when `eta < 2^-51`.

A double-only function receives identical bits for these two mathematically
different inputs and must return the same answer:

- if it returns `true`, it is wrong for `B_minus`;
- if it returns `false`, it cannot certify `B_plus`.

Changing one epsilon cannot recover information discarded by rational-to-double
conversion. A rigorous one-sided floating certificate would need tracked input
rounding bounds and interval/error-bound arithmetic. That is substantially more
complex than the exact rational PD check already available in FracESSA.

## Correct Resolution

The minimal correctness-first resolution is the implemented one:

1. Do not construct `Bee_dbl` for stability.
2. Do not allow a floating PD result to declare ESS.
3. Construct the exact rational Bee matrix directly.
4. Use exact LDL^T signs and the existing exact copositivity path.

The previously measured cost was approximately 10% for repeated small matrices
and approximately 1% over the 35 historical benchmark matrices. That cost buys
a deterministic mathematical decision rather than a compiler-dependent guess.

Related files:

- `../experiments/speed_comparison_2026-07-26/MICROBENCHMARK_COMPARISON.md`:
  persistent-process speed comparison.
- `../../testdata/fracessa_testdata.sqlite3`: regression inputs 36 and 37.
