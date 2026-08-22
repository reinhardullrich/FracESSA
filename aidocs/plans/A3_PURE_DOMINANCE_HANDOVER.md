# A3 Pure-Dominance Experiment Handover

Status: active experiment. The singleton-boundary search-order change is implemented and verified. Production generator behavior is
unchanged; the shared generators now expose one read-only boundary query used by A3.

## Scope

Model A3 lives under `models/a3/` and extends A2. It is a safe-only experimental FracESSA implementation. A3 combines:

- A2's optional out-of-order full-support inertia calculation;
- exact ascending support enumeration;
- exact upper-cone pruning after a curvature failure or an exact Nash equilibrium; and
- exact pure-by-pure weak-dominance pruning.

The experiment is intentionally isolated from production, packaging, releases, and the canonical database.

## Implemented dominance rule

Strategy \(i\) is weakly dominated by pure strategy \(j\ne i\) when

\[
A_{jk}\ge A_{ik}\qquad\text{for every strategy }k.
\]

A3 compares the complete rows of the parser's exact integer matrix. If this condition holds, no ESS support can contain strategy
\(i\), so the singleton support \(\{i\}\) becomes a forbidden root. The existing generator then avoids every support containing it.

For a circular game, one dominated singleton and all of its rotations cover the complete singleton orbit. A3 therefore registers one
circular forbidden orbit instead of inserting every rotated singleton separately.

## Current A3 order

The current implementation in `models/a3/src/fracessa.cpp` performs these steps:

1. Analyze every singleton support exactly. Circular games analyze the one singleton representative and retain its orbit multiplier.
2. Register every exact singleton equilibrium as a pending forbidden root through the normal support generator.
3. Ask the generator whether at least two strategies remain outside those pending singleton roots. If not, finish immediately because
   no support of size two or larger can survive.
4. Otherwise run the binary64 full-support probe, even when a pending singleton root already forbids the full support.
5. If the probe requests exact work, factor the complete reduced Hessian exactly. A full-support ESS ends the search; otherwise exact
   negative inertia may establish a global cardinality ceiling.
6. Scan the full game for weak pure-by-pure dominance.
7. Register dominated singleton roots.
8. Start cardinality two. The generator activates all pending roots and skips their supersets.

`full_support_checked` means that the full support was checked exactly. A binary64 probe alone does not set it. If no earlier pruning
prevents the generator from eventually emitting the full support, its callback performs the exact check when it is still outstanding.

## Verification already completed

The current A3 binary was built after its latest source change. The following checks passed:

- `python3 models/a3/verify.py`: all 23 direct and database cases matched production ESS counts and ESS support-size structures;
- a non-circular dimension-65 multiword dominance case matched production;
- a circular dimension-65 zero game exercised the multiword circular dominance orbit and matched production with 65 singleton
  candidates and no ESS; and
- no ESS-count mismatch occurred in the 1,071-matrix benchmark described below.

Candidate output is not a compatibility target for A1, A2, or A3 because exact curvature and dominance can reject unstable Nash
candidates before their probability and outside-payoff calculations. ESS counts and structures remain the correctness target.

## Original pre-fix A2 versus A3 benchmark

Results are stored in `experiments/model_timings.sqlite3` under session `a2_a3_dominance_20260822`.

Protocol:

- exact 1,071-matrix selection reused from session `timing_20260821T232438.846587Z`;
- every selected matrix had a positive production-safe calibration below 100 seconds;
- dimension two was excluded by the maintained runner;
- persistent Pybind workers;
- parent process pinned to CPU 1;
- workers pinned to CPUs 2 through 9;
- 0.5-second calibration target per matrix; and
- 30-second timeout per matrix.

The current A2 binary had changed since the previous stored A2 run, so both current binaries were measured afresh. Their hashes were:

- A2: `360b4554d169e577ce1535e515d51a048b48b29f46e5bb0ae0c880ea76e5bbaa`;
- A3: `ef6edc3c412642d870218b2c39c591d2c927e6a812ae6b893f217d458480faaf`.

Correctness and completion:

| Model | Correct | Mismatch | Timeout | Error |
|---|---:|---:|---:|---:|
| A2 safe | 1,069 | 0 | 2 | 0 |
| A3 safe | 1,063 | 0 | 8 | 0 |

On the 1,062 matrices completed by both models:

- median A3/A2 runtime was `1.225598`, so A3 took 22.56% more time at the median;
- A3 was faster on 222 matrices, slower on 838, and tied on 2;
- the sum of successful paired native medians was 9.14% lower for A3, but this excludes all timeout rows and must not be interpreted as
  an overall win; and
- replacing every timeout by only its 30-second lower bound gives 118.50 seconds for A2 and 293.16 seconds for A3.

An independent exact scan of the same matrix strings found weak pure-by-pure dominance in 311 matrices:

| Category | Matrices | Shared completions | Median A3/A2 | Sum A3/A2 | A3 wins/losses | A2/A3 timeouts |
|---|---:|---:|---:|---:|---:|---:|
| Dominance present | 311 | 310 | 1.094198 | 0.133944 | 121 / 188 | 1 / 0 |
| No dominance | 760 | 752 | 1.506762 | 0.977875 | 101 / 650 | 1 / 8 |

The dominance rule therefore gives enormous wins on a small expensive subset but adds overhead on most matrices. Examples include:

- ID 2424, Deriv2 order 21: 3.163 seconds to 52.25 microseconds, a 60,534x speedup;
- ID 2405, Deriv2 order 19: 962.5 milliseconds to 46.54 microseconds, a 20,682x speedup;
- ID 2074, order 18: 11.02 milliseconds to 16.63 microseconds, a 663x speedup; and
- ID 2453, Deriv2 order 24: A2 timed out after 30 seconds while A3 completed in 105.9 microseconds.

## Corrected A3 rerun against stored A2

Results are stored under session `a2_a3_dominance_fixed_20260822`. A2 was not rerun: its 1,071 rows were copied unchanged from the
original session after confirming that the current A2 binary still has the exact stored hash. Corrected A3 alone was rerun on the same
matrix selection with the same persistent-worker, CPU-affinity, 0.5-second target, and 30-second timeout protocol.

- reused A2: `360b4554d169e577ce1535e515d51a048b48b29f46e5bb0ae0c880ea76e5bbaa`;
- corrected A3: `2c5ae0907a71d9b76fea9d2b3fde51ddf690197863c2fe15ce9b9f318792de7b`.

| Model | Correct | Mismatch | Timeout | Error |
|---|---:|---:|---:|---:|
| Stored A2 safe | 1,069 | 0 | 2 | 0 |
| Corrected A3 safe | 1,070 | 0 | 1 | 0 |

On all 1,069 shared completions:

- median A3/A2 runtime was `1.108577`, so A3 took 10.86% more time at the median;
- A3 was faster on 247 matrices, slower on 821, and tied on 1;
- summed paired native medians were 2.76% lower for A3 (`0.972415` A3/A2); and
- replacing each timeout by only its 30-second lower bound gives 118.50 seconds for A2 and 86.89 seconds for A3.

The same exact dominance classification gives:

| Category | Matrices | Shared completions | Median A3/A2 | Sum A3/A2 | A3 wins/losses | A2/A3 timeouts |
|---|---:|---:|---:|---:|---:|---:|
| Dominance present | 311 | 310 | 1.029922 | 0.137906 | 132 / 178 | 1 / 0 |
| No dominance | 760 | 759 | 1.127826 | 1.047042 | 115 / 643 | 1 / 1 |

Because the A2 rows are deliberately reused rather than contemporaneously rerun, these ratios compare identical binaries and protocol
on the same machine but do not remove run-to-run system variation.

## Resolved full-support timeout regression

Before the singleton-boundary fix, A3 introduced seven timeout rows that A2 completed:

| IDs | Order | Family | Stored candidate and ESS structure | A2 time | A3 time |
|---|---:|---|---|---:|---:|
| 110, 111, 112 | 48 | Pothen mesh | `{"1":48}` | about 23 microseconds | over 30 seconds |
| 2635 | 54 | Lehmer | `{"1":54}` | 46.3 microseconds | over 30 seconds |
| 2641 | 55 | Lehmer | `{"1":55}` | 49.9 microseconds | over 30 seconds |
| 2658 | 57 | Lehmer | `{"1":57}` | 55.3 microseconds | over 30 seconds |
| 2663 | 58 | Kms | `{"1":58}` | 35.2 microseconds | over 30 seconds |

ID 2493 timed out in both models. A2 additionally timed out on ID 2453, which A3 solved through dominance pruning.

The seven new timeouts are not caused by the row-dominance scan. Every strategy is already a singleton equilibrium and ESS. Their
singleton forbidden roots cover every nonempty support, so the normal generator would finish before cardinality two. A2 avoids the
full-support calculation and finishes in microseconds. A3 nevertheless runs its binary64 full-support probe first; that probe requests
an exact full reduced-Hessian factorization. Fraction-free intermediate integers become sufficiently expensive that the exact
calculation exceeds 30 seconds.

After the fix, all seven rows again match their stored ESS counts and structures. The corrected maintained rerun measured their A3
native medians from 22.8 to 52.7 microseconds.

The exact full-support calculation can still be useful when the full support itself is forbidden but other supports remain. Its
inertia ceiling may then remove many of those remaining supports. The original bug was narrower: A3 paid for the full support even
when no support remained at all.

## Implemented completion behavior

The implementation retains this order and distinction:

1. Analyze the singleton layer and register singleton-equilibrium pruning.
2. Before the full-support probe, determine whether singleton-equilibrium pruning has already left no support to inspect.
3. If the search is empty, finish immediately. Do not run the full-support probe and do not run the later dominance scan because the
   complete answer is already known.
4. If some support remains, run the full-support binary64 probe and any exact calculation it requests, even when the full support
   itself is already forbidden.
5. Only after that full-support work, run the exact pure-dominance scan and register its roots.
6. Continue normal enumeration from cardinality two.

The user suggested the simple completion condition "the number of singleton ESS equals the game order." Internally, the generator's
actual pruning rule is based on exact Nash-equilibrium supports, not ESS status. The precise condition is that fewer than two
strategies remain outside the pending singleton roots: with zero or one available strategy, no size-two or larger support exists.

The implementation does not rely on an analyzer-side count. Each generator variant answers through the read-only
`has_supports_after_singletons()` query, which does not consume or mutate generator state.

## Rejected proposal

Do not implement the previously proposed shortcut of moving dominance before the full-support probe and then waiting for the first
cardinality-two callback. The user explicitly rejected it because dominance must remain after the full-support calculation. Using
dominance roots to decide whether the pre-probe search is already empty changes the requested order.

The implemented query resolves the earlier open design without consuming the first size-two support, repeating generation, or moving
dominance. Focused generator tests cover one-word, circular, and multiword behavior. The normal Release build, all eight CTest entries,
all 76 Python tests, and all 23 A3 verification cases pass.
