# Choice 1 Candidate Filter Reimplementation Plan

Status: future-phase planning only; do not implement until the temporary
unsafe-default phase is complete, reviewed, and accepted.

Branch: `choice-one-candidate-filter`

Base: `main` at `32f61679da64beb30f36870e190538f9d80e5970`

Read-only reference: committed Choice 1 implementation at
`9a33dab8b06865571ef4ec6ca7f0b2aa6a6af6b2`. The implementation will be
re-created manually from `main`; it will not merge or cherry-pick that branch.

Preceding phase: `UNSAFE_CANDIDATE_FILTER.md` specifies the temporary
unsafe-default heuristic and the single validating matrix parser. Implement
that phase first. This document must then be re-audited against the accepted
source before its file list is treated as implementation scope.

When Choice 1 later lands, it becomes the no-flag numerical mode. Explicit
`--unsafe` continues to select the already-implemented heuristic, and `--exact`
continues to bypass both numerical filters. Unsafe code must not enter or weaken
the strict proof kernel described here.

## Objective

Add a rigorous one-sided Choice 1 certificate as the future no-flag filter while
preserving the explicit unsafe heuristic, existing exact solver, and search
flow.

The numerical filter has exactly two outcomes:

```text
CERTIFIED_REJECT -> skip the exact candidate solve
EXACT_REQUIRED   -> run the unchanged exact rational candidate solve
```

It never accepts a candidate numerically. Singular systems, unsupported
floating-point behavior, overflow, non-finite values, failed bounds, and every
interval touching a decision boundary must return `EXACT_REQUIRED`.

Priorities remain:

1. Correctness: no exact candidate may be rejected numerically.
2. Speed: the filter runs millions or billions of times on small systems.
3. Readability: implementation names and formulas should correspond directly.

## Mathematical Contract

For support `S` of size `k`, with bordered dimension `m = k + 1`, solve

```text
        [ A_S  -1 ] [ x ]   [ 0 ]
C z  = [           ] [   ] = [   ].
        [ 1^T   0  ] [ u ]   [ 1 ]
```

Production guarantees `1 <= k <= n <= 63`, hence `2 <= m <= 64`. Internal
proof helpers may rely on this caller contract; do not add repeated dimension
checks to the per-support hot path.

The exact candidate conditions are

```text
x_j > 0                                         for every j in S
sum(j in S) A(i,j)*x_j - u <= 0                 for every i outside S.
```

An outside gain equal to zero is not a rejection. The filter may return
`CERTIFIED_REJECT` only after proving either

```text
x_j <= 0
```

for a support probability, or

```text
sum(j in S) A(i,j)*x_j - u > 0
```

for an outside strategy.

## Choice 1 Algorithm

### 1. Lazy exact normalization

On the first non-exact candidate check for a game, use exact rational arithmetic
to compute

```text
c       = A(0,0)
s       = max(i,j) |A(i,j)-c|
A'(i,j) = (A(i,j)-c)/s.
```

Translation by `c` and positive scaling by `s` preserve candidate supports,
extended supports, payoff comparisons, and ESS decisions. If `s == 0`, disable
the filter for that analyzer and use exact arithmetic.

Initialization is lazy so `--exact` never allocates or prepares binary64 filter
state.

### 2. Rigorous rational enclosures

Include `<mpfr.h>` before the FLINT declaration of `fmpq_get_mpfr`. Convert each
exact normalized entry `q` independently in the same direction at both stages:

```text
lo  = double_RNDD(mpfr_RNDD(q))
mid = double_RNDN(mpfr_RNDN(q))
hi  = double_RNDU(mpfr_RNDU(q)).
```

Require finite `lo <= mid <= hi`. If `lo == hi`, store `mid = lo` and zero
radius. Otherwise store

```text
radius    = max(upward(mid-lo), upward(hi-mid))
magnitude = max(|lo|,|hi|).
```

Store only three `n*n` binary64 matrices:

```text
A0    midpoint
RA    radius satisfying |A'-A0| <= RA
Aabs  magnitude satisfying |A'| <= Aabs.
```

MPFR is already a required dependency. No interval library or new dependency is
needed. Reuse one explicitly initialized 53-bit MPFR temporary for the complete
game conversion rather than initializing one temporary per entry.

### 3. Strict binary64 proof arithmetic

Compile the proof source separately from the fast project target:

```text
IPO/LTO:      disabled for the proof object
GCC/Clang:    -fno-fast-math -ffp-contract=off
MSVC:         /fp:strict /GL-
```

Outward operations use one explicitly grouped nearest-rounded operation followed
by one adjacent binary64 step. Use the previously tested `double` to `uint64_t`
bit representation approach from the reference implementation rather than
`std::nextafter` in the hot path. Handle signed zero explicitly.

This is the intended C++17 stepping code; keep it direct rather than wrapping it
in an interval type:

```cpp
bool round_down(double value, double& result) noexcept
{
    if (!std::isfinite(value)) return false;
    if (value == 0.0) {
        result = -std::numeric_limits<double>::denorm_min();
        return true;
    }

    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    if (value < 0.0) ++bits;
    else --bits;
    std::memcpy(&result, &bits, sizeof(result));
    return std::isfinite(result);
}

bool round_up(double value, double& result) noexcept
{
    if (!std::isfinite(value)) return false;
    if (value == 0.0) {
        result = std::numeric_limits<double>::denorm_min();
        return true;
    }

    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    if (value > 0.0) ++bits;
    else --bits;
    std::memcpy(&result, &bits, sizeof(result));
    return std::isfinite(result);
}
```

Check the binary64 format, round-to-nearest mode, and subnormal preservation once
per analyzer during lazy initialization. Do not put environment checks in the
per-support arithmetic loops. On an unrecognized compiler, build a conservative
stub that always requires the exact solver instead of compiling an unverified
floating-point proof.

The analyzer executes synchronously on one thread. After successful lazy
initialization, its caller must not change that thread's rounding or subnormal
environment during the analyzer lifetime. Each later pybind call constructs a
new analyzer and checks its own calling thread.

### 4. Compact retained LU

Extract support and complement indices once into fixed 64-element stack arrays.
Build the compact `m*m` midpoint bordered matrix and factor it in place with
partial pivoting:

```text
upper triangle        U
strict lower triangle L multipliers
permutation[i]        original row now occupying factor row i.
```

Keep the factors and permutation. Solve the permuted right-hand side with one
fixed stack vector. Do not allocate a solution matrix, construct a dense inverse,
or use the old absolute pivot thresholds. A pivot swap must exchange complete
rows, including multipliers retained from earlier columns.

The reusable compact LU matrix may resize once when support size changes. There
must be no heap allocation for each support.

### 5. Cheap rejection proposals

Before proof work, identify only possible violations from the midpoint solution:

```text
solution[j] <= 0       possible support-probability rejection
midpoint gain > 0      possible outside-strategy rejection.
```

These signs are proposals only. They never reject a support. If no proposal
exists, use exact arithmetic immediately instead of calculating the proof.

### 6. Residual enclosure

Treat each finite binary64 solution component as an exact real number. Let `C0`
be the midpoint bordered matrix and `R_C` its entrywise input radius, whose only
nonzero block is

```text
R_C[p,q] = RA(support_indices[p],support_indices[q])   for p,q < k.
```

First enclose

```text
r0[p] = sum(q<k) A0(support_indices[p],support_indices[q])*solution[q]
        - solution[k]                                      for p < k
r0[k] = sum(q<k) solution[q]-1

r_input = R_C*|solution|
r_lo    = downward(r0_lo-r_input)
r_hi    = upward(r0_hi+r_input).
```

Apply the row permutation only after constructing the unpermuted intervals, so
`[r_lo,r_hi]` contains `P*(C*solution-b)`. Reconstruct active bordered
coefficients directly from the stored game matrices; do not retain a second
bordered matrix or materialize a per-support radius matrix.

### 7. Choice 1 LU error proof

Use the underflow-aware Oishi-Rump LU factorization-defect bound to construct

```text
u       = 2^-53
eta     = 2^-1074
rho     = 1-m*u
gamma   = m*u/rho
delta   = m/rho

d_lu    = gamma*|L|*(|U|*1)
          + delta*(m*1+diag(|U|))*eta
d_input = P*(R_C*1)
d        = upward(d_lu+d_input)
         >= |P*C-L*U|*1.
```

The mathematical source is Oishi and Rump, *Fast Verification of Solutions of
Matrix Equations* (2002), Theorem 4.1. The existing
`../correctness/CERTIFIED_CANDIDATE_FILTER.md` contains the broader certificate
derivation and references.

Apply upward absolute triangular recurrences to bound

```text
q >= ||(L*U)^-1 * (P*C-L*U)||_infinity.
```

For each nonnegative input vector `v`, the recurrence is exactly

```text
y[i] = upward(v[i] + sum(j<i) |L(i,j)|*y[j])
w[i] = upward((y[i] + sum(j>i) |U(i,j)|*w[j]) / |U(i,i)|),
       with the second line evaluated for i = m-1,...,0.
```

Thus `w >= |U^-1|*|L^-1|*v`; no inverse or generic triangular-solve abstraction
is required.

Only `q < 1` proves nonsingularity. Apply the same recurrence to the enclosed
residual and compute

```text
r_abs[i] = max(|residual_lower[i]|, |residual_upper[i]|)
beta     = max(abs_lu_apply(r_abs))
denom    = downward(1-q)
error    = upward(beta/denom)
         >= ||z-solution||_infinity.
```

Any failed operation or `q >= 1` means `EXACT_REQUIRED`.

### 8. Certified rejection

For midpoint support-probability proposals, reject only when

```text
upward(solution[j] + error) <= 0.
```

For midpoint outside-gain proposals, calculate only the needed lower bound:

```text
g0_lo = downward(sum(j in S) A0(i,j)*solution[j]-solution[k])
gain_radius = upward(
    sum(j in S) (RA(i,j)*|solution[j]| + Aabs(i,j)*error) + error
)
gain_lo = downward(g0_lo-gain_radius).
```

The final `+error` covers the payoff component. Reject only when `gain_lo > 0`.
Equality and overlap with zero require the exact solver.

## Future Production Source Scope

The strict mathematics and proof-kernel boundaries below are complete, but this
source list was reviewed against the stated `main` base. The preceding unsafe
phase changes several of the same files. Re-run caller and dead-code searches
before requesting source approval; do not mechanically apply the old removals.

The expected C++ production delta is:

| File | Minimal change |
| --- | --- |
| `cpp/CMakeLists.txt` | Add one strict proof object and one test-only exact-rejection oracle option. |
| `cpp/include/fracessa/fracessa.hpp` | Add only the numerical mode state and certified-filter declaration needed beside the existing unsafe and exact paths. |
| `cpp/include/fracessa/matrix_server.hpp` | Add lazy Choice 1 state while preserving the accepted unsafe matrix and scratch state. Remove no symbol merely because it was unused on the old `main` base. |
| `cpp/src/findeq.cpp` | Preserve the unsafe filter and add the certified-rejection call plus compile-time exact oracle. |
| `cpp/src/fracessa.cpp` | Select Choice 1 by default, unsafe only when explicitly requested, and exact when requested, with one predictable branch per support. |
| `cpp/src/certified_candidate_filter.cpp` | Add the isolated Choice 1 implementation. |
| `cpp/src/certified_candidate_filter_detail.hpp` | Declare only the private helper contracts needed by focused production-helper tests. |

The matching C++ test scope is:

| File | Minimal change |
| --- | --- |
| `cpp/tests/CMakeLists.txt` | Register one focused certificate test executable. |
| `cpp/tests/test_certified_candidate_filter.cpp` | Test the actual strict production helpers and MatrixServer integration. |

The future public-mode handoff is expected to require the existing CLI, pybind,
and wrapper boundary files as well: pass explicit unsafe selection into the
core and add a numerical `unsafe` Python argument only when Choice 1 creates
that distinction. Matrix parsing remains one validated path. Re-audit and list
those exact files before seeking implementation approval.

No Choice 1 change is expected in `cpp/include/linalg/linear_solver.hpp`; the
temporary unsafe phase owns removal of the obsolete generic double solver, and
the exact solver remains unchanged.

If another `.cpp`, `.hpp`, or `.py` file becomes necessary after the required
re-audit, implementation must stop and request renewed approval before changing
it.

## Explicitly Excluded

The reimplementation will not include:

- Choice 2 or a compile-time Choice 1/Choice 2 selector;
- unsafe heuristics or parser routing inside the strict proof kernel; the
  separate public mode is specified in `UNSAFE_CANDIDATE_FILTER.md`;
- certificate diagnostics, counters, or destructor output;
- multiprocessing, sink, or output-schema changes;
- GitHub workflow or release changes;
- copied benchmark dumps, generated sources, or experimental search programs;
- a generic interval type, inverse matrix, virtual dispatch, or public
  certificate API.

The temporary phase retains `find_candidate_dbl()` for the unsafe heuristic.
Choice 1 should add one narrowly named certified-filter entry instead of
renaming or duplicating the existing unsafe implementation. Both numerical
filters keep the same control-flow shape:

```text
false -> the selected filter rejects the support
true  -> exact candidate solving is required
```

Only Choice 1 may describe its `false` result as certified. Unsafe retains the
same boolean shape but not the proof guarantee.

## Comment Requirements

Comments in new source must explain only non-obvious mathematics and machine
contracts:

- the one-sided rejection/fallback invariant;
- why exact affine normalization preserves the game decisions;
- why both MPFR conversion stages use the same directed rounding;
- how adjacent binary64 stepping creates outward bounds;
- retained LU storage and permutation semantics;
- the residual enclosure;
- the Oishi-Rump defect inequality, absolute triangular recurrence, and `q < 1`;
- the final probability and outside-gain rejection inequalities.

Do not comment simple assignments, ordinary loops, or self-explanatory branches.
Keep existing variable names such as `support_indices`, `non_support_indices`,
`support_size`, `solution`, `payoff`, and `rowsum` wherever their meaning remains
the same.

## Verification Gates

### Baseline

Before source edits, record the accepted temporary unsafe-default build and test
state. Verification IDs 1-44 should already pass there; IDs 38 and 39 are no
longer expected failures after unsafe normalization.

### Focused correctness

The focused test must call the actual strict production helpers and cover:

- adjacent binary64 stepping and unsupported environments;
- rigorous rational conversion boundaries;
- retained LU reconstruction, solve, and multiple row swaps;
- singular, near-singular, non-finite, and maximum `m = 64` fallback;
- exact upper bounds for the LU defect and absolute triangular recurrence;
- exact containment by residual and outside-gain bounds;
- strict and boundary support-probability decisions;
- strict and boundary outside-gain decisions;
- exact affine transformations and constant-matrix fallback;
- deterministic exact-oracle cases that must never be falsely rejected.

### End-to-end correctness

Build with `FRACESSA_CERTIFICATE_ORACLE=ON`. Every certified rejection is then
cross-checked with the exact candidate solver and must fail immediately on any
disagreement. The option must compile completely out of normal builds.

Run:

```text
focused certificate and exact-solver CTests
complete core/CLI CTests
complete wrapper tests
all verification-matrix correctness tests
ASan/UBSan core and CLI checks
```

Do not run verification matrix 33 or 34 with `--exact` or the exact-rejection
oracle without Reinhard's separate approval for that run.

### Strict-build inspection

Inspect the generated compile command and verify that the proof source has
strict floating-point flags, contraction disabled, and IPO disabled. The normal
project sources retain their existing speed flags.

### Performance

Use the existing persistent-process three-second benchmark protocol on the same
matrix set and machine conditions as the saved comparison. Exclude matrix 34.
Do not rerun historical binaries or matrix 33 exact. Compare per-matrix medians
with the saved raw-double, exact-only, and previous Choice 1 measurements.

No diagnostics implementation is planned. If timing differs unexpectedly,
instrumentation requires a separate measured reason and separate approval.

## Acceptance Conditions

The implementation is acceptable only when:

1. No exact candidate is ever rejected by the oracle build.
2. All included verification matrices match their exact candidate baselines.
3. `--exact` does not initialize or allocate certificate state.
4. The normal support path performs no per-support heap allocation except the
   existing scratch resize when support size changes.
5. No-flag behavior becomes Choice 1, explicit unsafe remains available, and
   exact, parser, result, and output-schema behavior otherwise remain unchanged.
6. The strict proof source remains isolated from fast-math, contraction, and IPO.
7. The final source diff contains only the approved files and documented code.

## Required Regression Data

`main` contains verification IDs 1-44. Commit
`2be0207242d585aaa14f9c87bbaed1b068c40de0` on `certified-filter` contains the
validated data for IDs 45-47:

- ID 45 is the preserved dimension-20 heuristic counterexample;
- ID 46 reaches the Choice 1 LU error proof and requires boundary fallback;
- ID 47 reaches residual construction but fails the Choice 1 LU proof and
  requires exact fallback.

Transfer only their committed rows during the future Choice 1 implementation
into
`python/verification/verification_matrices.json`,
`python/verification/baseline_candidates.csv`, and
`python/verification/baseline_result.json`. Do not activate them while unsafe is
still the temporary default. The candidate CSV must retain the
mathematically correct exact candidates; the result JSON intentionally remains
the historical wrong-filter speed baseline. Do not copy their experimental
builds, logs, source snapshots, or Choice 2/unsafe code.

This is required test-data scope, not an expansion of the approved C++ source
scope.
