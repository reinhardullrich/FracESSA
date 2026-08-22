# Model A4: Exact Invader-Interval Pruning

Status: implemented and retained as an isolated research model. Production is unchanged.

## Exact Rule

Let `S` be a support whose exact internal solution `x` has positive probabilities and common support payoff `u`. Define the exact
invader set

\[
V(x)=\{j\notin S:(Ax)_j>u\}.
\]

If `V(x)` is nonempty, every ESS support containing `S` must contain at least one member of `V(x)`. Equivalently, A4 may reject any
later support `T` satisfying

\[
S\subsetneq T\subseteq N\setminus V(x).
\]

The rule is exact. It uses the fraction-free candidate solution and exact integer payoff comparisons; floating-point arithmetic does
not create or apply a pruning rule.

## Retained Implementation

A4 copies A2 and changes only its isolated safe path:

1. The exact candidate solver returns the complete invader set after an otherwise internally valid singleton or pair candidate fails
   the outside-payoff test.
2. The non-circular support generator stores the clause `(S,V)` until the next cardinality.
3. Fixed-order descending-bit generation indexes the clause by the last required bit of `S`. The clause becomes active only after a
   branch contains all of `S`, and the branch is rejected if no selected or future bit can satisfy `V`.
4. Only singleton and pair clauses are retained. Later clauses measured narrower and cost more to inspect than the exact solves they
   avoided.
5. Circular games use A2's unchanged bracelet generator. Circular symmetry expansion and representative-only clauses both measured
   slower.

The implementation supports both the one-word and multiword support paths. Omitting a clause only forgoes optional pruning; there is
no fallback and no change to the exact answer.

## Discarded Designs

The direct fixed-order index is the smallest design that survived measurement. Three broader implementations were removed:

- a dynamic binary decision diagram made representative cases up to roughly 1,163 times slower;
- a general two-watched-literal index made matrix 89 roughly 320 times slower;
- expanding circular clauses over rotations, reflections, and affine symmetries made matrix 34 roughly 1,927 times slower.

No SAT solver, generalized constraint engine, hash index, or materialized support frontier was added.

## Verification and Timing

The maintained ordinary, circular, affine-circular, multiword, and database verification cases all matched A2's ESS count and
structure.

The full comparison selected the 1,071 stored matrices with `0 < safe_calibration_ns < 100,000,000,000`, excluding dimension 2. It
used a 0.5-second target, a 30-second timeout, a scheduler on CPU 1, and persistent workers on CPUs 2 through 9. Unchanged A2 rows
were reused from the stored timing database; A4 alone was rerun.

| Group | Paired matrices | Median A4/A2 | Difference |
|---|---:|---:|---:|
| All | 1,069 | 1.026 | +2.57% |
| Non-circular | 714 | 1.071 | +7.11% |
| Circular | 355 | 1.000 | 0.00% |
| Small | 172 | 1.071 | +7.11% |
| Medium | 359 | 1.024 | +2.40% |
| Large | 295 | 1.008 | +0.79% |
| Super-large | 243 | 1.026 | +2.62% |

A4 completed 1,070 matrices correctly, with no mismatch and one timeout. A2 completed 1,069 correctly and timed out twice. The most
important improvement was matrix 2453: A2 exceeded 30 seconds, while A4 finished in 0.169 ms. Matrices 2424, 52, and 2430 also fell
from seconds to fractions of a millisecond. The largest absolute regression was matrix 211 at +0.330 seconds; the largest ratio
regression was matrix 2169, from 0.374 ms to 29.196 ms.

## Decision

Retain A4 as a high-variance research model. The complete corpus shows a 2.57% median cost, but the exact rule removes one timeout and
collapses several multi-second searches. Do not promote it to production without a separate production decision and a canonical
single-CPU benchmark.
