# Historical Fast Candidate False Rejections

Status: correctness reference for historical failures and the non-certifying boundary that still applies: fast remains heuristic;
safe is the exact method. Reproduction results below are a dated snapshot, not current benchmark data.

Last verified: 2026-08-01 against historical fast revision `8697ebaf` and the current working tree

## Result

Fast revision `8697ebaf` can reject a real candidate and ESS even when none of its six matrix-wide input checks requests safe
fallback. One counterexample exists for each of its three per-support rejection tests:

1. a nonsingular system can produce a pivot below the fixed cutoff;
2. a strictly positive exact probability can be computed as negative; and
3. an exact outside payoff below the equilibrium payoff can be computed above the permitted margin.

All three matrices are exact, pass those six matrix-wide input checks, and lose a real ESS in historical fast. Current fast sends
small pivots to exact checking and selects matrix-wide safe search for $P\geq10^9$. It therefore recovers all three stored
regressions. Its remaining double inequality rejections are still heuristic, so this is not a general correctness proof.

## Counterexample 1: Pivot Cutoff

Let

$$
\delta=10^{-12}
$$

and

$$
A=
\begin{pmatrix}
-3 & 1 & 2\\
1 & -\dfrac{1+\delta}{3} & \dfrac{-2+\delta}{3}\\
2 & \dfrac{-2+\delta}{3} & -\dfrac{4+\delta}{3}
\end{pmatrix}.
$$

FracESSA's exact upper-triangle input is

```text
3#-3,1,2,-1000000000001/3000000000000,-1999999999999/3000000000000,-4000000000001/3000000000000
```

The actual results are:

| Method | Candidates | ESS |
|---|---:|---:|
| historical `fast` | 0 | 0 |
| current `fast` | 1 | 1 |
| `safe` | 1 | 1 |

The `safe` result is

```text
1;1/3,1/3,1/3;7;3;7;3;;1;T_reduced_hessian_nd;0;0.000000
```

## Why It Is an Exact Candidate

Every row of $A$ sums to zero. Therefore

$$
x=
\begin{pmatrix}
1/3\\
1/3\\
1/3
\end{pmatrix}
$$

satisfies

$$
Ax=0\mathbf 1.
$$

All three probabilities are strictly positive. The support is the full set, so there are no outside strategies to check. The
payoff is $u=0$.

## Why It Is an ESS

Put $L=-A$. A tangent direction preserves the sum of the probabilities, so it has the form

$$
z=(p,q,-p-q).
$$

Direct substitution gives

$$
z^TLz=
\begin{pmatrix}p&q\end{pmatrix}
\frac{1}{3}
\begin{pmatrix}
25+\delta & 5+2\delta\\
5+2\delta & 1+4\delta
\end{pmatrix}
\begin{pmatrix}p\\q\end{pmatrix}.
$$

The first leading entry is positive, and the determinant of this two-dimensional reduced matrix is

$$
\frac{(25+\delta)(1+4\delta)-(5+2\delta)^2}{9}
=9\delta
>0.
$$

Thus $L=-A$ is positive definite on the tangent space. Equivalently, $A$ is negative definite there. Consequently the exact
full-support candidate is an ESS.

## Why the Input Checks Allow `fast`

The matrix entries are finite, normal binary64 values of ordinary size. No exact nonzero becomes zero, no value becomes infinity
or NaN, the exponent range is small, and no two distinct exact entries become the same double. The distinct entry values are
separated by roughly one third or more, far above the $10^{-12}$ adjacent-entry check.

The dangerous small number appears only after elimination through cancellation. It is not present as a small matrix entry or as
the difference between two entries. The matrix-wide checks therefore do not select `safe` fallback.

## Where `fast` Rejects It

For the full support, `fast` performs partial-pivot Gaussian elimination on

$$
\begin{pmatrix}
A & -\mathbf 1\\
\mathbf 1^T & 0
\end{pmatrix}
\begin{pmatrix}x\\u\end{pmatrix}
=
\begin{pmatrix}\mathbf 0\\1\end{pmatrix}.
$$

The first pivot has magnitude $3$. After that elimination, partial pivoting selects the normalization row with pivot $4/3$.
The next pivot has exact magnitude

$$
\frac{3\delta}{4}=7.5\times10^{-13}.
$$

The historical fast cutoff is $10^{-12}$. Because

$$
7.5\times10^{-13}<10^{-12},
$$

the historical `find_candidate_fast::find()` returns `false` before back-substitution. The support is discarded and never reaches
the exact candidate solver. Current fast returns candidate-possible here so exact arithmetic decides the support.

This is a cutoff failure, not a floating-point conversion failure. The exact bordered system is nonsingular, but it is close
enough to singular that the fixed fast cutoff rejects it.

## Counterexample 2: Positive Probability

This example still fails after every pivot below $10^{-12}$ is sent to the safe solver. Its exact upper-triangle input is

```text
3#1/50000000,4001/200000000000,5000000001/100000000,1/50000000,10000000002001/200000000000,0
```

The exact full-support result is

```text
1;1/10000000000,1/2,4999999999/10000000000;7;3;7;3;;1;T_reduced_hessian_nd;25000000010002500001/1000000000000000000;25.000000
```

| Pivot-fallback method | Candidates | ESS |
|---|---:|---:|
| historical `fast` | 0 | 0 |
| current `fast` | 1 | 1 |
| safe | 1 | 1 |

Let

$$
\varepsilon=10^{-10},\qquad \delta=10^{-13},\qquad M=100.
$$

Using strategy 3 as the reference coordinate, the exact reduced matrix is

$$
H=
\begin{pmatrix}
A_{11}-2A_{13}+A_{33} & A_{12}-A_{13}-A_{23}+A_{33}\\
A_{12}-A_{13}-A_{23}+A_{33} & A_{22}-2A_{23}+A_{33}
\end{pmatrix}
=-M
\begin{pmatrix}
1&1\\
1&1+\delta
\end{pmatrix}.
$$

The vector

$$
x=(\varepsilon,1/2,1/2-\varepsilon)
$$

solves the exact equal-payoff equations. Its entries are strictly positive. Moreover,

$$
H_{11}=-M<0,
\qquad
\det(H)=M^2\delta>0,
$$

so $H$ is negative definite. Because this is the full support, the exact candidate is an ESS.

In an isolated Release build, both small-pivot branches were changed only from rejection to safe fallback. That build still
returned zero candidates. Instrumenting the next branch showed

```text
probability i=0 value=-0.00031170391640434274
```

Thus the bordered binary64 solve turns the exact value $x_1=10^{-10}$ into a negative approximation far below the
$-10^{-10}$ rejection margin. No pivot fallback occurred; otherwise the exact solver would have retained the support. Returning
safe fallback from the probability branch alone restores the exact ESS.

## Counterexample 3: Outside Payoff

This four-strategy example also fails after small pivots are sent to the safe solver. Its exact upper-triangle input is

```text
4#1/5,40000000001/200000000000,501/10,6040000000001/400000000000,1/5,10020000000001/200000000000,10048000000001/400000000000,0,10040000000001/400000000000,0
```

The safe solver finds two ESS. The one lost by the pivot-fallback fast variant is

```text
1;1/1000,1/2,499/1000,0;7;3;7;3;;1;T_reduced_hessian_nd;10040040000001/400000000000;25.100100
```

| Pivot-fallback method | Candidates | ESS |
|---|---:|---:|
| historical `fast` | 1 | 1 |
| current `fast` | 2 | 2 |
| safe | 2 | 2 |

Write its payoff as

$$
u=\frac{10040040000001}{400000000000}
$$

and its support probabilities as

$$
x=(1/1000,1/2,499/1000).
$$

The reduced matrix of the first three strategies is again

$$
H=-100
\begin{pmatrix}
1&1\\
1&1+10^{-13}
\end{pmatrix},
$$

so the support equations have the stated strictly positive solution and $H$ is negative definite. The fourth strategy's three
payoffs are

$$
\left(u-\frac{100001}{10000},\;u+\frac{199}{10000},\;u-\frac{1}{10000}\right).
$$

Their exact expected payoff is

$$
\begin{aligned}
&\frac{1}{1000}\left(u-\frac{100001}{10000}\right)
+\frac{1}{2}\left(u+\frac{199}{10000}\right)
+\frac{499}{1000}\left(u-\frac{1}{10000}\right)\\
&=u-\frac{1}{10000}<u.
\end{aligned}
$$

Therefore the fourth strategy is not an outside best reply, the extended support remains 7, and negative definiteness proves
that this exact candidate is an ESS.

The pivot-fallback fast variant omits support 7 and reports only the other ESS on support 14. Instrumenting the outside-payoff
branch showed, for support 7,

```text
outside row=3 rowsum=25.1100200000025 threshold=25.100500000002498 payoff=25.100100000002499
```

The exact outside payoff is $u-10^{-4}$, but the binary64 candidate error makes the computed outside payoff about $0.00992$ above
the computed equilibrium payoff and about $0.00952$ above the rejection threshold. The probability test precedes this branch and
did not fire for this support. No pivot fallback occurred; otherwise the exact solver would have retained it. Returning safe
fallback from the outside-payoff branch restores the missing ESS.

## Current Reproduction

From the repository root after a Release build:

```bash
cpp/build/fracessa fast '3#-3,1,2,-1000000000001/3000000000000,-1999999999999/3000000000000,-4000000000001/3000000000000' --candidates
cpp/build/fracessa safe '3#-3,1,2,-1000000000001/3000000000000,-1999999999999/3000000000000,-4000000000001/3000000000000' --candidates

cpp/build/fracessa fast '3#1/50000000,4001/200000000000,5000000001/100000000,1/50000000,10000000002001/200000000000,0' --fullsupport --candidates
cpp/build/fracessa safe '3#1/50000000,4001/200000000000,5000000001/100000000,1/50000000,10000000002001/200000000000,0' --fullsupport --candidates

cpp/build/fracessa fast '4#1/5,40000000001/200000000000,501/10,6040000000001/400000000000,1/5,10020000000001/200000000000,10048000000001/400000000000,0,10040000000001/400000000000,0' --candidates
cpp/build/fracessa safe '4#1/5,40000000001/200000000000,501/10,6040000000001/400000000000,1/5,10020000000001/200000000000,10048000000001/400000000000,0,10040000000001/400000000000,0' --candidates
```

Current fast and safe now agree on these commands. Reproducing the historical mismatches requires revision `8697ebaf`. The first
two historical mismatches also occur with `--fullsupport`. The outside-payoff example must enumerate support 7, so it should be
run without that option.

## Relation to Earlier Regression Matrices

- Former database IDs 38 and 39 exercised exact-to-double input problems and triggered matrix-wide safe fallback. Both were later
  removed as strategically equivalent duplicates.
- IDs 45, 46, and 47 were false rejections caused by the retired normalized heuristic. The current raw fast method solves them.
- The first matrix above exposes the former per-support pivot rejection.
- The second and third matrices show that pivot fallback alone leaves independent false rejections in both remaining per-support
  inequalities; current fast's precision-span gate sends both complete matrices to safe search.

These counterexamples are stored in the canonical SQLite matrix set as IDs 2207, 2208, and 2209 respectively.

## Precision-Span and Small-Pivot Mitigation

Current `fast` makes two changes relative to historical fast revision `8697ebaf`:

1. a pivot below $10^{-12}$ returns candidate-possible so exact arithmetic decides that support; and
2. one exact integer precision-span check replaces the six conversion checks, and the whole matrix uses safe search when
   $P\geq10^9$.

For the precision span, let $d$ be the least common positive denominator and let $Z=dA$ be the integer game. Fast uses

$$
M=\max\left(d,\max_{i,j}|Z_{ij}|\right),
$$

$$
m=\min\left(d,\min_{Z_{ij}\neq0}|Z_{ij}|,\min_{Z_{ij}\neq Z_{kl}}|Z_{ij}-Z_{kl}|\right),
$$

and

$$
P=\frac{M}{m}.
$$

The three counterexamples above now give fast/safe ESS counts of $1/1$, $1/1$, and $2/2$. The pivot example is recovered by the
per-support pivot fallback. The probability and outside-payoff examples exceed the precision-span cutoff and therefore use
matrix-wide safe search. This remains a heuristic: it is not a proof that every matrix with $P<10^9$ is safe for double rejection.
