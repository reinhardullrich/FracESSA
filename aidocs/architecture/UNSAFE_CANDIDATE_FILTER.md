# Temporary Unsafe-Default Candidate Filter

Status: historical; local correctness and performance verification is recorded below.

Document role: failed-approach history. The normalized heuristic fixed raw-double failures 38-39 but introduced failures 45-47;
it was not a monotonic safety improvement and was removed. At the end of this phase, production required either `fast`, the
historical unnormalized raw-double heuristic with matrix-wide safe fallback checks, or `safe`, the exact-only path. The current
method contract is in [`../KNOWLEDGE.md`](../KNOWLEDGE.md); reusable numerical counterexamples are in
[`../correctness/FAST_CANDIDATE_FALSE_REJECTION.md`](../correctness/FAST_CANDIDATE_FALSE_REJECTION.md). The remainder of this
document is historical.

Branch: `choice-one-candidate-filter`

Base: `main` at `32f61679da64beb30f36870e190538f9d80e5970`

Read-only reference: the implemented heuristic at
`2be0207242d585aaa14f9c87bbaed1b068c40de0` on `certified-filter`. Re-create
only the behavior specified here; do not merge or copy Choice 1, Choice 2,
diagnostics, experiments, generated results, or unrelated branch changes.

## Temporary Decision

The current phase has exactly two numerical modes:

```text
unsafe   fast heuristic candidate rejection with exact fallback
exact    no floating-point candidate rejection
```

There is no third numerical mode called `default`. While this temporary phase
is active, omitting a numerical flag selects `unsafe`. Writing `--unsafe`
selects the same behavior explicitly. `--exact` selects the existing exact-only
path.

Choice 1 is a later phase. Do not add its proof kernel, strict floating-point
target, state, routing, tests, or regression data in this phase.

Correctness warning: the unsafe filter is not a rejection proof. Its cheap danger
test catches many suspicious floating-point rejections and sends them to exact
arithmetic, but it can still reject a real candidate. Preserved reference cases
45-47 demonstrate this limitation.

## Public Option Contract

### CLI

```text
no numerical flag       unsafe numerical mode (temporary)
-u, --unsafe            unsafe numerical mode, selected explicitly
-e, --exact             exact numerical mode
```

Reject `--exact --unsafe` because the user explicitly selected both numerical
modes. The short flag `-u` belongs to the numerical `--unsafe` option. Matrix
input always uses the single validating parser.

Every non-exact CLI run should identify the numerical mode as unsafe on
`stderr`, whether selected implicitly or explicitly. This warning is outside
the timed analyzer construction and does not enter the per-support hot path. It
must state both that exact candidates and ESS results can be missed and that
suspicious or unusable floating cases fall back to exact work and can therefore
be much slower.

### Pybind And Python

Keep the existing numerical selector:

```text
exact=False   -> unsafe numerical mode
exact=True    -> exact numerical mode
```

Do not add a parser selector or a second numerical `unsafe` boolean yet. The
latter would be redundant while unsafe is the only non-exact mode. Add that
numerical distinction only when Choice 1 becomes the no-flag default.

## Existing Search Routing

The core already expresses the complete temporary mode selection:

```cpp
if (!conf_exact_)
    if (!find_candidate_dbl(support, support_size))
        return;

if (!find_candidate_frc(support, support_size))
    return;
```

Therefore this phase does not add `conf_unsafe_`, a new constructor parameter,
a new filter-selection branch, or `find_candidate_unsafe()`. Improve the body
of the existing `find_candidate_dbl()` and keep its contract:

```text
false -> unsafe rejection; skip this support
true  -> run the unchanged exact candidate solver
```

The unsafe path may reject a support or request exact arithmetic. It never
accepts a candidate and never declares an ESS. Exact stability checking remains
unchanged.

## One-Time Conditional Exact Normalization

Before support enumeration, and only when `conf_exact_` is false, compute from
the exact rational matrix:

```text
c       = A(0,0)
s       = max(i,j) |A(i,j)-c|
A'(i,j) = (A(i,j)-c)/s.
```

Translation by `c` and positive scaling by `s` preserve support probabilities,
best-response comparisons, candidates, extended supports, and ESS decisions.
Every normalized entry is in `[-1,1]`.

If `s == 0`, return false from initialization so the analyzer uses exact
candidate solving for every support. Convert each normalized value once with
the existing `fraction::to_dbl()`. A non-finite conversion, or a nonzero
rational value that converts to zero or a subnormal double, returns false in
the same way. The subnormal check is required because the relative-error model
behind the epsilon guard is not valid below the smallest normal binary64 value.

Implement `MatrixServer::initialize_game_matrix_dbl()` and call it from the
analyzer constructor before the first support. The method only prepares the
normalized matrix and reports whether it is usable. If it returns false, the
analyzer sets its existing `conf_exact_` state to true so the unchanged support
loop uses exact candidate solving thereafter.

This is conditional initialization, not a first-support getter check. It keeps
`--exact` free of double allocation while avoiding initialized/available state,
pointer checks, and an extra branch on every support. It also establishes a
clear precondition: every double matrix reference is obtained only after
successful initialization, so no reference can be invalidated by later matrix
allocation.

Keep normalization in `MatrixServer` with the other matrix-preparation code.
The server owns and returns the exact game, the finished normalized double
game, and reusable matrix scratch, but contains no filter or fallback decision.

Reuse the existing `game_dbl_` and `sub_bordered_dbl_` storage. The uncalled
`get_bee_matrix_dbl()` and `bee_dbl_` remain removed. Do not expand this into
cleanup of unrelated matrix helpers or tests.

FLINT's `fmpq_get_d()` conversion used by `fraction::to_dbl()` rounds toward
zero. Using binary64 epsilon `2^-52`, rather than only unit roundoff `2^-53`, in
the heuristic guard is a conservative heuristic allowance for that conversion
direction; it is not a proof bound.

## Fused Unsafe Kernel

For support size `k`, set bordered dimension `m = k+1` and build

```text
        [ A'_S  -1 ] [ x ]   [ 0 ]
        [ 1^T     0 ] [ u ] = [ 1 ].
```

Implement the historical small-matrix solve directly in
`find_candidate_dbl()`:

1. Extract support and non-support indices once into fixed 64-byte stack arrays.
2. Ask `MatrixServer` only for compact `m*(m+1)` scratch storage, keyed by
   `support_size`; build the complete bordered system in `find_candidate_dbl()`
   from the already extracted indices.
3. Perform destructive Gaussian elimination with physical partial-pivot swaps.
4. Complete back substitution into `double solution[64]`.
5. Scan outside rows directly through the extracted indices.

The scratch accessor must not extract the support or populate matrix entries;
doing so would repeat work already performed by `extract_support_partition()`.
Every one of the `m*(m+1)` entries is overwritten before elimination, so stale
scratch values cannot affect the next support.

Do not allocate a solution matrix, materialize a full-length double strategy,
or add a generic solver API. After this single caller contains the fused solve,
remove the unused `solve_linear_dbl()` and only its tests and now-unused
includes. Leave the exact solver unchanged.

Retain the historical thresholds only as rejection proposals:

```text
selected pivot < 1e-12                 -> exact fallback
x[j] < -1e-10, j < k                   -> probability proposal
outside_gain > 1e-4 * game_dimension   -> outside proposal
```

The final variable `solution[k]` is payoff and receives no positivity test.
Small, zero, or non-finite pivots mean uncertainty, so they request exact
arithmetic. Complete the solve before classifying a proposal so the danger test
has the retained pivot and solution scales it needs.

## Cheap Danger Veto

Retain during the solve:

```text
minimum_pivot  = minimum absolute selected pivot
solution_scale = max(1, |solution[0]|, ..., |solution[m-1]|)
epsilon        = binary64 epsilon = 2^-52
H              = 64
```

For a probability proposal:

```text
margin         = -solution[j]
decision_scale = solution_scale.
```

For an outside proposal:

```text
margin         = outside_gain
decision_scale = m * solution_scale.
```

Use the measured reference comparison directly:

```cpp
constexpr double kGuardMultiplier = 64.0;
constexpr double kBinary64Epsilon = std::numeric_limits<double>::epsilon();
const double base_risk =
    kGuardMultiplier * static_cast<double>(m) * kBinary64Epsilon;

const auto is_suspicious = [&](double margin, double decision_scale) noexcept {
    const double risk = base_risk * decision_scale;
    const double separation = minimum_pivot * margin;
    return !std::isfinite(risk) || !std::isfinite(separation) ||
           separation <= 0.0 || risk >= separation;
};
```

For each proposal:

```text
suspicious      -> return true  -> exact candidate solve
well separated  -> return false -> unsafe rejection
```

Stop at the first proposal. If it is suspicious, go exact immediately. If no
proposal exists, return `true` and run the exact solver. Do not add residuals,
matrix norms, inverses, condition estimators, a second solve, or a configurable
guard multiplier.

`H = 64` is a private heuristic, not a mathematical error bound. This veto
improves the old decision but does not make it correct for every matrix.

The arithmetic assumes an IEEE-754 binary64 `double`. Enforce that once at
compile time with `std::numeric_limits<double>` (`is_iec559`, radix 2, and 53
significand bits). Do not add floating-environment checks to the support loop;
this unsafe mode is not a proof kernel.

## Implementation Scope

The implementation is limited to these source files.

### C++ Core And Boundary

| File | Minimal change |
| --- | --- |
| `cpp/include/fracessa/matrix_server.hpp` | Prepare, store, and return the normalized double game and reusable scratch; contain no filter or fallback decision. |
| `cpp/include/linalg/linear_solver.hpp` | Remove only `solve_linear_dbl()` after its single caller is fused; preserve the exact solver. |
| `cpp/src/findeq.cpp` | Own the fused unsafe solve, proposals, and danger veto while consuming the normalized matrix from `MatrixServer`. |
| `cpp/src/fracessa.cpp` | Request normalized unsafe state once before enumeration; on preparation failure, reuse `conf_exact_` for permanent exact fallback. |
| `cpp/src/main.cpp` | Add numerical `--unsafe`, reject exact+unsafe, use the single parser, and print the unsafe warning. |
| `cpp/src/pybind_module.cpp` | Use the single parser and keep `exact` as the numerical selector. |
| `cpp/tests/test_matrix.cpp` | Check that nonzero normalization values converting to zero or subnormal binary64 make double-matrix preparation fail. |
| `cpp/tests/test_linear_solver.cpp` | Remove only tests for the deleted double solver. |
| `cpp/tests/cli_parser_blackbox.py` | Verify parser routing, numerical routing, warning text, normalization cases, and incompatible CLI options. |

### PyFracESSA

| File | Minimal change |
| --- | --- |
| `python/fracessa/types.py` | Expose no parser selector in `RunConfig`. |
| `python/fracessa/core.py` | Call pybind without a parser selector. |
| `python/fracessa/tests/test_core_unit.py` | Assert no parser selector reaches pybind. |

### Documentation

The implementation task updated `README.md`, `aidocs/KNOWLEDGE.md`, this architecture document, `aidocs/pyfracessa/README.md`, and
`aidocs/CHANGES.md`. The former review finding was later resolved; the durable counterexamples now live in
`aidocs/correctness/FAST_CANDIDATE_FALSE_REJECTION.md`. This phase improved an explicitly unsafe mode and did not make the default
mathematically correct.

The CLI and Python documentation must say plainly that no flag and
`exact=False` use a heuristic filter that can miss candidates and ESS
results, and that conservative fallback can have exact-mode runtime.

No change is planned for `cpp/include/fracessa/fracessa.hpp`, CMake, the
constructor signature, output schemas, multiprocessing, sinks, or dependencies.
If another `.cpp`, `.hpp`, or `.py` file becomes necessary, stop and request
renewed source approval before editing it.

## Explicitly Excluded

- Choice 1 and Choice 2.
- A third numerical mode or generalized filter interface.
- `conf_unsafe_`, a numerical pybind/Python `unsafe` boolean, or a new core
  constructor parameter.
- A candidate-rejector-double source file, strict floating-point object target, MPFR
  proof arithmetic.
- Claiming that the danger veto proves rejection.
- Per-support heap allocation, diagnostics, counters, or logging.
- Per-support normalization-state or availability checks.
- A public `H` setting or automatic calibration.
- Compatibility aliases for the unreleased name `--unsafe-parsing`.
- Changes to exact solving, exact stability, candidate output, or support order.
- Activating reference verification IDs 45-47 in the normal suite while unsafe
  remains the temporary default.

## Verification Plan

### Interface And Routing

1. No numerical flag and explicit `--unsafe` produce identical analyzer output.
2. Both non-exact CLI forms print the unsafe warning.
3. The CLI test proves that `--exact` prints no unsafe warning. Source and diff
   inspection must separately confirm that this path never calls
   `MatrixServer::initialize_game_matrix_dbl()` and leaves both double matrices
   unallocated;
   do not add allocator instrumentation for this one structural property.
4. Every numerical mode uses the same validating parser.
5. `--exact --unsafe` fails before parsing or analysis.
6. Pybind and `RunConfig` expose `exact`; no parser selector or redundant
   numerical boolean is added.
7. Explicit `--unsafe` still rejects dimension 64; the analyzer contract
   remains `1 <= n < 64`.

### Numerical Behavior

1. Verification IDs 1-44 must match their accepted candidate and ESS baselines
   in temporary default unsafe mode.
2. IDs 38 and 39 specifically verify that exact affine normalization fixes the
   old scale and translation failures.
3. Compare allowed unsafe outputs with `--exact`, but never run IDs 33 or 34
   exact without Reinhard's separate approval.
4. Reuse the existing edge fixtures instead of adding another test harness:
   ID 7 checks zero-payoff candidates, ID 41 checks a negative payoff, ID 44
   checks a zero-payoff mixed ESS, while IDs 40 and former ID 43 exercise constant-matrix exact
   fallback, and ID 42 preserves an outside-payoff equality boundary.
5. Keep reference IDs 45-47 outside the active verification files for now.
   They are known negative controls where unsafe loses an exact ESS, not tests
   that the temporary default can pass.
6. The focused MatrixServer test must prove that a nonzero normalized rational
   converting to zero or a subnormal double makes initialization return false.
7. Run normal C++, CLI, wrapper, matrix, and ASan/UBSan coverage.

### Performance

Use the established optimized persistent-process three-second protocol:

```text
Release + native CPU flags + IPO/LTO where supported
CPU affinity and candidates enabled
verification IDs 1-33 and 35-44
matrix 34 excluded
no exact rerun of matrix 33
per-matrix medians
```

Compare against the saved historical raw-double and reference unsafe results in
commit `2be0207242d585aaa14f9c87bbaed1b068c40de0` at
`experiments/speed_comparison_2026-07-26/results/unsafe_filter_vs_choice1_3s_2026-07-27.{json,csv}`.
Read those files from the Git ref; do not copy generated experiment data into
this branch.

The reference unsafe implementation was 6.7% slower in aggregate than the old
incorrect raw filter on their shared subset, but that aggregate hides important
variation: IDs 1-12 were roughly 5.8% to 32% slower, while veto-triggered exact
fallback made IDs 19, 29, and 33 about 3.64x, 8.21x, and 1.77x slower. It was
still materially faster than Choice 1 on every measured matrix.

Report medians separately for IDs 1-12 and for the complete benchmark set.
Because this project runs millions of small matrices, the all-matrix sum cannot
hide a small-matrix regression. The saved reference was measured in an earlier
run, so a difference above 5% is a screening trigger, not by itself proof of a
code regression. It blocks acceptance until the build, machine conditions, and
implementation difference explain it; only a separately approved same-session
interleaved rerun can support a precise regression claim. Do not add a more
elaborate heuristic merely to meet the threshold.

Do not create a no-veto mode or permanent diagnostics merely to isolate a few
scalar operations. Add temporary instrumentation only if the complete result
materially disagrees with the saved reference.

## Acceptance Conditions

This temporary phase is acceptable when:

1. Unsafe and exact are the only numerical modes, and unsafe is selected when
   no numerical flag is supplied.
2. All entry surfaces use one validating matrix parser.
3. Every failed or suspicious unsafe check delegates to the exact candidate
   solver.
4. IDs 1-44 match their accepted exact candidate and ESS output.
5. IDs 45-47 remain documented evidence that unsafe is not proof-backed and are
   not activated as default-mode regressions yet.
6. Exact solving, exact stability, output schema, and support ordering remain
   unchanged.
7. Exact mode never calls unsafe initialization and leaves the double game and
   scratch matrices empty.
8. The unsafe support path allocates only when its reusable scratch shape
   changes, never once for every support, and remains near the saved unsafe
   benchmark in both the small-matrix and complete groups.
9. The source diff stays within the exact files listed above.

## Implementation Result

Verified locally on 2026-07-28:

- Release build succeeded with GCC 16.1.1, native CPU flags, and IPO/LTO.
- All 10 core/CLI CTests, all 24 wrapper tests, and all 44 active matrix
  correctness tests passed. IDs 38-39 now return their accepted normalized
  results.
- Default and `--exact` candidate output matched byte-for-byte for the 42 active
  matrices allowed to run exact; exact IDs 33 and 34 were not executed.
- A separate ASan/UBSan build passed the same 10, 24, and 44 test groups.
- Source inspection confirms `--exact` skips
  `MatrixServer::initialize_game_matrix_dbl()`; the server starts with empty
  double matrices, and only the unsafe path can allocate the double game or
  bordered scratch.
- The existing three-second persistent-process harness measured IDs 1-33 and
  35-44 on CPU 2, with matrix 34 excluded. Every ESS and candidate count matched
  the saved unsafe run. Summed medians changed by `+0.052%` for IDs 1-12 and
  `+0.629%` over the complete set.
- ID 38 alone crossed the 5% historical screening threshold: `3.121 us` versus
  `2.896 us` (`+7.8%`, `+0.226 us`). Both binaries use the same compiler and
  runtime libraries, but the saved measurement was paired/interleaved with
  Choice 1 and its binary contained that additional code, while this standalone
  LTO binary does not. Aggregate groups are stable, so this is recorded as an
  isolated historical-run/code-layout difference, not a paired regression
  claim.
- The source diff remained within the approved files. Exact candidate solving,
  exact stability, output schema, support order, multiprocessing, sinks, CMake,
  and dependencies were not changed.

The known IDs 45-47 limitation remains. Passing IDs 1-44 does not make the
unsafe danger veto a correctness proof.

## Historical Outcome

The later verified proof was implemented and then removed. At the end of this historical phase, production exposed only the
required `fast` and `safe` methods described at the top.
