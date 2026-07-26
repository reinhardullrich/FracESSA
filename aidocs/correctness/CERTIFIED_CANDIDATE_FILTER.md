# Certified Candidate Filter Design Note

Status: design discussion only. Nothing in this document is implemented.

## Purpose

FracESSA may inspect up to `2^n` supports. Running exact rational Gaussian
elimination for every support is potentially too slow, but the current double
filter is not allowed to discard supports safely: verification matrices 38 and
39 contain valid ESS supports that the default path rejects before exact
arithmetic.

The target is a one-sided filter:

```text
CERTIFIED_REJECT  -> safely skip exact candidate solving
EXACT_REQUIRED    -> make no numerical decision; run the exact solver
```

The certificate is allowed to save work only when it proves rejection. It is
never allowed to turn uncertainty into rejection. In particular:

```text
certificate succeeds and proves a violated condition -> reject
certificate fails, is too wide, or overlaps a boundary -> exact fallback
```

`EXACT_REQUIRED` is not an error. It is the normal and correct result for a
singular, nearly singular, or insufficiently separated support system.

## Source Of The Certificate

The proposed scalar error bound is a standard a posteriori verification result
for linear systems. The primary reference is:

- Siegfried M. Rump, *Verification methods: rigorous results using
  floating-point arithmetic*, Acta Numerica 19 (2010), 287-449,
  [DOI 10.1017/S096249291000005X](https://doi.org/10.1017/S096249291000005X),
  Section 10.3, Theorem 10.2 and equation (10.9). An
  [author-hosted PDF](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf)
  is available.

Rump's theorem states that, for a square real matrix `C`, a right-hand side
`b`, an approximate solution `z_hat`, and any matrix `R`, the condition

```text
||I - R C|| < 1
```

proves that `C` is nonsingular and gives

```text
||C^-1 b - z_hat||
    <= ||R (b - C z_hat)|| / (1 - ||I - R C||).
```

The theorem is exact mathematics. It does not assume that `R` really is the
inverse, that `z_hat` is accurate, or that `C` is well conditioned. A poor
`R` merely makes the inequality fail, in which case the theorem says nothing.

For a practical implementation using only round-to-nearest floating-point
arithmetic, the most directly relevant reference is:

- Takeshi Ogita, Siegfried M. Rump, and Shin'ichi Oishi, *Verified Solution of
  Linear Systems Without Directed Rounding*, Technical Report 2005-04,
  [repository record](https://tore.tuhh.de/entities/publication/8ae7c0ad-441d-45da-be13-d6421d5f1c01),
  [PDF](https://www.tuhh.de/ti3/paper/rump/OgRuOi05z.pdf), especially Theorems
  3.1-3.2 and Section 4. That paper derives computable upper bounds, including
  rounding and underflow, without repeatedly changing the processor's rounding
  mode.

The stronger interval fixed-point alternative belongs to the Krawczyk family:

- Rudolf Krawczyk, *Newton-Algorithmen zur Bestimmung von Nullstellen mit
  Fehlerschranken*, Computing 4 (1969), 187-201,
  [DOI 10.1007/BF02234767](https://doi.org/10.1007/BF02234767).
- Siegfried M. Rump, *Validated Solution of Large Linear Systems*, Computing
  Supplementum 9 (1993), 191-212,
  [PDF](https://www.tuhh.de/ti3/paper/rump/Ru93.pdf), Theorem 1.1.

The first FracESSA experiment should use the simpler scalar norm bound. A full
Krawczyk interval iteration should be added only if measurements show that the
scalar bound sends too many otherwise easy rejections to exact arithmetic.

## The Exact FracESSA System

For a support `S` of size `k`, FracESSA solves a bordered system of dimension
`m = k + 1`:

```text
C_S z = b,

      [ A_S  -1 ]          [ x ]          [ 0 ]
C_S = [           ],   z = [   ],     b = [   ].
      [ 1^T    0 ]          [ u ]          [ 1 ]
```

Here:

- `A_S` is the payoff submatrix restricted to `S`;
- `x` contains the `k` support probabilities;
- `u` is the common equilibrium payoff on `S`;
- the first `k` equations say `A_S x = u * 1`;
- the final equation says `sum(x_i) = 1`.

The exact candidate conditions are:

```text
x_j > 0                                      for every j in S
sum(j in S) A(i,j) x_j - u <= 0              for every i outside S
```

An outside value equal to zero is not a rejection. It means that the outside
strategy is another best response and belongs to the candidate's extended
support. Exact equality therefore remains important even if a numerical filter
can prove many strict rejections.

## Exact Affine Normalization

The normalization constant `c` is not an equilibrium payoff. It is an
arbitrary exact scalar taken directly from the input matrix, preferably

```text
c = A(0,0).
```

No game or support has to be solved to obtain it. Define, using rational
arithmetic,

```text
d(i,j) = A(i,j) - c
s      = max(i,j) |d(i,j)|
A'     = d / s,                    when s > 0.
```

If `s = 0`, every matrix entry is equal and the numerical filter should be
bypassed in the first implementation.

For any probability vector `x`, `sum(x_j) = 1`, so every pure-strategy payoff
transforms as

```text
(A' x)_i = ((A x)_i - c) / s.
```

Consequently, for two strategies `x` and `y`,

```text
x^T A' y = (x^T A y - c) / s.
```

Subtracting the common scalar changes every compared payoff by the same amount,
and division by the positive `s` preserves every inequality. Therefore this
normalization preserves best responses, candidate supports, extended supports,
and ESS conditions. If the original equilibrium payoff is `u`, the normalized
payoff is

```text
u' = (u - c) / s.
```

That payoff is produced later by the normalized support solve. It is not needed
to construct the normalization.

The exact normalization is done once per input game in `O(n^2)` rational
operations. The resulting entries lie in `[-1,1]`. The original rational matrix
remains available for exact fallback and final output.

This removes the avoidable scale and common-offset failures:

```text
matrix 38: [0, 10^-20; 10^-20, 0] -> [0, 1; 1, 0]
matrix 39: [10^20, 10^20+1; ...]  -> [0, 1; 1, 0]
```

Normalization improves the numerical problem, but it is not a certificate.
Matrices containing only small integers can still produce singular or
ill-conditioned bordered systems.

## Derivation Of The Bound

Let the exact normalized bordered system be

```text
C z = b.
```

Let `z_hat` be the ordinary double approximation, and let `R` be any double
matrix intended to approximate `C^-1`. Mathematically, every finite double is
an exact real number, so `z_hat` and `R` may be treated as exact real data in the
proof. Define

```text
E = I - R C
r = b - C z_hat.
```

If a submultiplicative matrix norm proves `||E|| < 1`, the Neumann series

```text
(I - E)^-1 = I + E + E^2 + E^3 + ...
```

converges, and

```text
||(I - E)^-1|| <= 1 / (1 - ||E||).
```

Because `R C = I - E`, this also proves that the square matrices `R` and `C`
are nonsingular. Subtracting `C z_hat` from `C z = b` and multiplying by `R`
gives

```text
C (z - z_hat)       = r
R C (z - z_hat)     = R r
(I - E)(z - z_hat)  = R r
z - z_hat           = (I - E)^-1 R r.
```

Taking norms yields

```text
||z - z_hat||
    <= ||R r|| / (1 - ||E||)
     = ||R (b - C z_hat)|| / (1 - ||I - R C||).
```

For implementation, compute rigorous floating-point upper bounds

```text
alpha >= ||I - R C||_inf
beta  >= ||R (b - C z_hat)||_inf.
```

If `alpha < 1`, then

```text
e = beta / (1 - alpha)
```

is a uniform componentwise enclosure because

```text
|z_i - z_hat_i| <= ||z - z_hat||_inf <= e.
```

Thus every exact solution component lies in

```text
Z_i = [z_hat_i - e, z_hat_i + e].
```

The infinity norm is useful here because it immediately gives componentwise
bounds and needs only maximum absolute row sums.

## What Makes The Computation Rigorous

The algebra above is not enough by itself. Ordinary evaluations of `alpha`,
`beta`, and `e` in nearest-double arithmetic can round downward and invalidate
the claimed upper bounds. The implementation must establish all of the
following:

1. The floating representation encloses every original exact rational entry.
2. Every computed interval or analytic error bound encloses all rounding error.
3. The final `alpha`, `beta`, and `e` are rounded upward.
4. Overflow, underflow, non-finite values, and unsupported compiler behavior
   cause `EXACT_REQUIRED`, not rejection.

### Enclosing Exact Rational Inputs

Let `q` be one exact normalized rational entry. The filter needs binary64
endpoints `q_lo` and `q_hi` satisfying

```text
q_lo <= q <= q_hi.
```

Two reasonable implementation routes are:

1. Convert through MPFR once with downward and upward rounding.
2. Convert once, compare the resulting binary64 value back against `q` using
   exact rational arithmetic, and move the required endpoint with
   `nextafter`. One step is enough only when the conversion routine guarantees
   correctly rounded nearest conversion; otherwise continue until enclosure is
   proved.

This cost is paid once per matrix after exact normalization, not once per
support.

### Route A: Outward-Rounded Interval Arithmetic

Construct an interval matrix `C_box` containing the exact `C`. With `R` and
`z_hat` treated as point intervals, evaluate

```text
E_box = I - R C_box
v_box = R (b - C_box z_hat)
```

using outward-rounded interval operations. For an interval value
`X = [lo,hi]`, define

```text
mag(X) = max(|lo|, |hi|).
```

Then valid bounds are

```text
alpha = max(i) sum(j) mag(E_box(i,j))
beta  = max(i) mag(v_box(i)).
```

IEEE 1788.1-2017 specifies interval arithmetic over IEEE 754 binary64
endpoints; its [official standard page](https://standards.ieee.org/ieee/1788.1/6074/)
is the relevant interoperability reference. A standard-conforming interval
type is sufficient, but not required, if FracESSA can prove the same inclusion
property itself.

### Route B: Round-To-Nearest Error Bounds

Ogita, Rump, and Oishi derive formulas that bound the matrix products and
residual while leaving the processor in round-to-nearest mode. This can avoid
rounding-mode switches and may suit FracESSA's millions of small systems better.
Their bounds explicitly account for underflow; copying only the headline
formula while omitting those correction terms would not be a certificate.

The two routes establish the same mathematical premises. They should not both
be implemented initially. Prototype and benchmark the smallest one that is
rigorous on FracESSA's supported compilers.

### Compiler Requirements

The certificate code cannot be compiled under transformations that discard
IEEE floating-point semantics, such as unrestricted `-ffast-math`. Fused
multiply-add contraction is acceptable only when the error analysis accounts
for the actual operation. If a platform cannot provide the arithmetic behavior
assumed by the bound, the filter must be disabled there and exact solving used.

## Turning The Solution Bound Into A Rejection Proof

Assume `alpha < 1` and let

```text
X_j = [z_hat_j - e, z_hat_j + e]    for j in S
U   = [z_hat_k - e, z_hat_k + e]    for the payoff.
```

### Invalid Support Probability

An interior support candidate requires `x_j > 0`. Therefore

```text
upper(X_j) <= 0
```

proves that this support cannot be a candidate and permits
`CERTIFIED_REJECT`.

The converse is deliberately not used. `lower(X_j) > 0` proves positivity, but
FracESSA still needs exact candidate values and exact extended-support equality
for surviving candidates.

### Profitable Outside Strategy

For every `i` outside `S`, evaluate with rigorous enclosures

```text
G_i = sum(j in S) A'(i,j) X_j - U.
```

The exact candidate condition is `G_i <= 0`. Therefore

```text
lower(G_i) > 0
```

proves a profitable outside deviation and permits `CERTIFIED_REJECT`.

These cases do not prove rejection:

```text
G_i overlaps zero
upper(G_i) == 0
lower(G_i) <= 0 < upper(G_i)
```

They all require the exact solver. In particular, exact equality determines
the extended support and cannot be inferred from an interval that merely
contains zero.

## Proposed Per-Support Flow

The simplest complete design is:

```text
prepare once per game:
    construct exact A' by translation and positive scaling
    construct a rigorous binary64 enclosure of A'

for each support S:
    build midpoint bordered system C_mid
    run fast double LU to obtain z_hat

    if LU cannot produce a finite complete z_hat:
        return EXACT_REQUIRED

    evaluate the same cheap double candidate conditions as today

    if double does not propose rejection:
        return EXACT_REQUIRED

    construct/apply approximate inverse R
    rigorously bound alpha and beta

    if alpha >= 1 or any bound fails:
        return EXACT_REQUIRED

    e = upward_bound(beta / (1 - alpha))

    if some support probability has certified upper bound <= 0:
        return CERTIFIED_REJECT

    if some outside gain has certified lower bound > 0:
        return CERTIFIED_REJECT

    return EXACT_REQUIRED
```

This means the current double solver cannot treat a small pivot or an
approximately negative probability as a final Boolean rejection. The future
numerical stage must separate these concepts:

```text
numerical solve status
approximate candidate classification
certified rejection status
```

A small pivot produces no certificate and therefore goes exact. An approximate
negative probability may trigger an attempted certificate, but rejection occurs
only if the exact probability enclosure lies entirely at or below zero.

The current Gaussian elimination also overwrites its bordered matrix. A
certificate later needs the original normalized `C` for the residual and
`I - R C`. The implementation must either retain a small unmodified copy or
reconstruct `C` from the normalized game and support. This is a performance
choice, not a mathematical one.

## Why A Small Residual Alone Is Insufficient

The residual

```text
r = b - C z_hat
```

may be tiny even when `z_hat` is far from the true solution if `C` is badly
conditioned. The factor involving `R` and `1 - alpha` is what converts a
residual into a proof of forward error. Therefore none of these is a
certificate by itself:

```text
small residual
small pivot tolerance not triggered
large estimated reciprocal condition number
agreement between two ordinary floating-point solves
```

They can guide whether certification is likely to succeed, but they cannot
justify rejection.

## Why Direct Interval Gaussian Elimination Is Not The First Choice

Replacing every scalar in Gaussian elimination with a generic interval can
produce severe dependency overestimation: the same uncertain quantity appears
multiple times and is treated as if each occurrence varied independently.
Rump's 2010 review discusses this failure in Section 10.1, and the 1993 paper
gives concrete growth examples.

The proposed method instead:

1. obtains a fast ordinary approximation;
2. preconditions with `R`;
3. performs a separate, one-sided verification step.

That is usually both tighter and closer to the current high-speed solver.

## Cost Model

For a support of size `k`, let `m = k + 1`.

```text
ordinary LU solve                    O(m^3)
explicit approximate inverse R       O(m^3)
bound for I - R C                    O(m^3)
residual bound                        O(m^2)
outside-strategy checks               O((n-k)k)
```

The exact constant factors matter more than the asymptotic notation because
FracESSA performs millions or billions of operations on small matrices.

An explicit inverse is the easiest statement of the theorem, but it may make
the verified filter several times as expensive as the present double solve.
Possible later optimization is to derive or adopt a factorization-based
verified bound that reuses the LU factors without materializing all of `R`.
That is not the first implementation: it adds proof and code complexity and is
justified only if the simple certified prototype is correct but too slow.

The verifier should run only when the ordinary double calculation proposes a
rejection. Candidate-like supports already require exact arithmetic for exact
output and equality classification, so certifying their acceptance provides no
immediate benefit.

## Failure Semantics

Every exceptional numerical case maps to exact fallback:

| Situation | Result |
| --- | --- |
| Exact normalization has `s = 0` | `EXACT_REQUIRED` |
| Midpoint LU fails or produces a non-finite value | `EXACT_REQUIRED` |
| Input conversion cannot be enclosed | `EXACT_REQUIRED` |
| `alpha >= 1` | `EXACT_REQUIRED` |
| `beta` or `e` overflows | `EXACT_REQUIRED` |
| An interval overlaps the decision boundary | `EXACT_REQUIRED` |
| A probability interval has upper endpoint `<= 0` | `CERTIFIED_REJECT` |
| An outside-gain interval has lower endpoint `> 0` | `CERTIFIED_REJECT` |

Most importantly:

```text
alpha >= 1 does not prove that C is singular.
```

It says only that this `R`, this precision, and this bound did not certify the
system.

## Why There Is No Reliable Global "Good Matrix" Rule

An absolute rule such as "no entries below `10^-20`" cannot characterize
safety:

1. Positive scaling and common translation do not change ESS.
2. Small exact integers can cancel to a singular bordered support system.
3. Conditioning belongs to each support system, not just the complete game
   matrix.
4. The correctness decision depends on distance to a boundary, not conditioning
   alone. A moderately conditioned solve can still have a probability or
   outside gain extremely close to zero.
5. Proving that every support is safe beforehand would require examining
   essentially the same exponential support family.

Useful diagnostics are still possible:

- distinct exact normalized values collapse to one binary64 value;
- a nonzero normalized value rounds to zero;
- the estimated reciprocal condition number is small;
- certification gives `alpha` near or above one;
- the proof radius `e` is large relative to the decision margin.

The first two are input-level warnings. The last three are support-level
diagnostics. None replaces exact fallback. The meaningful safety statement is
not "this input looks normal," but "this particular rejection was certified."

## Validation Plan

Correctness must be established before speed is considered.

### Required Correctness Checks

1. Compare every result against the current `--exact` path.
2. During development, independently run the exact candidate solver for every
   `CERTIFIED_REJECT` and assert that it really rejects.
3. Keep verification matrices 38 and 39 as scale/translation regressions.
4. Add or generate exact rational cases with:
   - singular bordered systems;
   - nearly singular but nonsingular bordered systems;
   - zero, tiny positive, and tiny negative support weights;
   - outside gains that are exactly zero and immediately on either side;
   - exact rational values around binary64 rounding boundaries;
   - all-equal payoff matrices.
5. Apply exact transformations `a A + c 11^T` with `a > 0` and verify that the
   candidate and ESS results remain identical.
6. Run the checks on every supported compiler/runner because rounding,
   contraction, and underflow behavior are part of the certificate machinery.

The central test invariant is:

```text
No support accepted by the exact candidate solver may ever be
CERTIFIED_REJECT.
```

### Required Performance Measurements

Measure at least:

1. current default runtime as an incorrect speed reference;
2. current `--exact` runtime as the correctness reference;
3. exact normalization with the existing double filter, for isolated overhead;
4. normalized certified filtering;
5. total and per-support-size counts of:
   - proposed double rejections;
   - certified probability rejections;
   - certified outside-payoff rejections;
   - `alpha >= 1` failures;
   - boundary-overlap fallbacks;
   - exact solves.

Use repeated in-process measurements for small matrices so process startup does
not dominate. Compare total workload time as well as median per-matrix time;
FracESSA's actual cost is the sum across an exponential support search.

## Minimal First Experiment

Ponytail's minimal implementation boundary is:

1. exact affine normalization once per matrix;
2. one scalar infinity-norm certificate based on Rump's Theorem 10.2;
3. only `CERTIFIED_REJECT` and `EXACT_REQUIRED` outcomes;
4. no certified acceptance;
5. no full Krawczyk iteration;
6. no new public option until correctness and speed are measured.

If this version certifies enough rejections and remains faster than full exact
solving, stop. Add a tighter componentwise/Krawczyk enclosure or a
factorization-based verifier only when measurements identify the simple bound
as the limiting factor.

## Open Implementation Decisions

These decisions require experiments rather than assumptions:

1. MPFR-directed endpoint conversion versus exact comparison plus `nextafter`.
2. Outward-rounded interval products versus the Ogita-Rump-Oishi
   round-to-nearest bounds.
3. Explicit `R` versus a factorization-based verified solve.
4. Retaining a small unmodified bordered matrix versus reconstructing it after
   the destructive LU solve.

No choice changes the mathematical contract: rejection is permitted only after
rigorously proving the relevant exact inequality.

## References

1. S. M. Rump, *Verification methods: rigorous results using floating-point
   arithmetic*, Acta Numerica 19 (2010), 287-449,
   [DOI](https://doi.org/10.1017/S096249291000005X),
   [PDF](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf). Theorem 10.2
   is the direct source of the scalar residual/error certificate.
2. T. Ogita, S. M. Rump, and S. Oishi, *Verified Solution of Linear Systems
   Without Directed Rounding*, Technical Report 2005-04,
   [record](https://tore.tuhh.de/entities/publication/8ae7c0ad-441d-45da-be13-d6421d5f1c01),
   [PDF](https://www.tuhh.de/ti3/paper/rump/OgRuOi05z.pdf). Theorems 3.1-3.2
   restate the bound; Section 4 derives rigorous round-to-nearest estimates.
3. S. M. Rump, *Validated Solution of Large Linear Systems*, Computing
   Supplementum 9 (1993), 191-212,
   [PDF](https://www.tuhh.de/ti3/paper/rump/Ru93.pdf). Theorem 1.1 gives the
   interval fixed-point inclusion underlying stronger Krawczyk-style variants.
4. R. Krawczyk, *Newton-Algorithmen zur Bestimmung von Nullstellen mit
   Fehlerschranken*, Computing 4 (1969), 187-201,
   [DOI](https://doi.org/10.1007/BF02234767).
5. [IEEE 1788.1-2017](https://standards.ieee.org/ieee/1788.1/6074/), *IEEE
   Standard for Interval Arithmetic (Simplified)*, for interval operations with
   IEEE 754 binary64 endpoints.
