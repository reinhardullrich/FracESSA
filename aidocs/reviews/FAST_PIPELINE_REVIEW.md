# Fast Candidate Pipeline Review

Last verified: 2026-08-02

Scope: the complete binary64 candidate path from integer game preparation through
`find_candidate_fast`, plus the exact handoff that protects candidates which the fast path cannot reject. This is an actionable
checklist, not a record of every experiment. Correctness comes first; speed changes are retained only when they help a balanced
panel of small and large, circular and non-circular games while keeping the code readable.

Status:

- `[ ]` is open.
- `[x]` was reviewed and needs no change.

## Overall conclusion

The current algebra is coherent: it clears denominators once, applies one common binary scale, equilibrates the complete game,
builds the reduced symmetric candidate system, factors it with Bunch-Kaufman LDLT, recovers the probabilities and payoff, and
hands every surviving support to exact arithmetic.

The factorization and solve paths reviewed here appear internally correct, including their 1-by-1 and 2-by-2 pivots and
permutations. The current regression matrices and random comparisons agree with safe search.

Fast search is nevertheless a heuristic prefilter, not a proof. The precision-span cutoff and small-pivot fallback recover the
known counterexamples, but the remaining binary64 probability and outside-payoff rejections have no general error bound. A false
binary64 rejection can therefore still lose an exact candidate. The known historical examples and this boundary are documented
in [`FAST_CANDIDATE_FALSE_REJECTION.md`](../correctness/FAST_CANDIDATE_FALSE_REJECTION.md).

## Required order

1. Close the two correctness items below.
2. Apply the bit-identical simplifications one at a time.
3. Prototype rounding-changing optimizations only in `find_candidate_test`.
4. Promote a prototype only after exact-result comparison and balanced timing evidence.
5. Consider support-generation work only after the local candidate kernel is clean.

## Correctness first

### FP-C01: Non-finite accumulated values can reject instead of falling back

- [x] Add conservative finite checks after complete accumulations in `find_candidate_fast::find()` and the mirrored test path.

The solver already checks each recovered non-reference probability, but it does not check every value used by the later
inequalities:

- `reference_probability` can become non-finite after repeated subtraction;
- the accumulated and rescaled `payoff` can become non-finite;
- `threshold` can overflow;
- an accumulated and rescaled outside `rowsum` can become non-finite.

For example, `reference_probability == -infinity`, `payoff == -infinity`, or `rowsum == +infinity` can currently satisfy a normal
comparison in the wrong direction and return `false`, thereby discarding the support. The minimal safe behavior is to return
`true` whenever one of these completed values is non-finite, so exact arithmetic decides the support.

No current matrix reproduces this under the dimension, precision-span, equilibration, and pivot guards. This is still a real
logical gap because a fallback path must not interpret a non-finite result as mathematical evidence.

Required check: add one focused regression that reaches each new guard, or document why a guard cannot be reached if a reliable
construction is impossible. Then compare fast and safe candidate contracts over the retained database.

Completed 2026-08-02 for compilers preserving IEEE special-value behavior: both independent paths now return candidate-possible
after a non-finite reference probability, payoff, threshold, or outside payoff. No exact input is known that reliably reaches
these guards through the existing precision-span, equilibration, pivot, and solved-coordinate checks, so no artificial test hook
was added merely to inject an impossible internal state. The complete native and Python suites pass with the guards active. MSVC
Release now requests `/fp:precise`; Windows verification is deferred as recorded under FP-C03.

### FP-C02: Fast-mode regression coverage is incomplete and one comment is stale

- [x] Exercise `fast` directly on every stored historical false-rejection matrix.
- [x] Update the stale `cli_parser_blackbox.py` comment which says fast rejects the precision-span example.

The current black-box block invokes only `test` for the precision-span example even though the current fast implementation now
contains the same protection. The regression suite should explicitly cover:

- the small-pivot counterexample;
- the positive-probability counterexample;
- the outside-payoff counterexample;
- the collapsed-value and tiny-value matrix-wide fallbacks.

For each case, safe gives the reference candidate/ESS result and fast must match it. Tests should assert results, not merely that
the command exits successfully.

Completed 2026-08-02: the CLI black-box test now compares complete candidate output from `fast` and `test` byte for byte with
`safe` for the small-pivot, positive-probability, and outside-payoff counterexamples. It also checks their expected ESS counts of
1, 1, and 2. The focused black-box test and all 10 CTest tests pass; all 56 Python tests pass as well.

### FP-C03: MSVC Release must preserve special-value behavior

- [x] Replace the MSVC Release `/fp:fast` option with `/fp:precise`.

The fast and test paths now rely on standard NaN and infinity propagation plus `std::isfinite()` to select exact fallback.
`cpp/CMakeLists.txt` formerly applied `/fp:fast` to every MSVC Release source, even though Microsoft documents that special values
may not propagate or behave according to IEEE 754 under that option. MSVC Release now uses `/fp:precise`; Linux and macOS already
avoid an equivalent fast-math mode.

Closed 2026-08-02 as a configuration change. No Windows build or regression was run in this task, as requested. Compile the next
Windows release and run the CLI counterexample regressions before publishing it; that deferred platform verification does not keep
this source-review point open.

## Bit-identical or structural simplifications

These changes do not intentionally alter the arithmetic. They should still be tested because they touch the hot solver.

### FP-S01: Reuse the 1-by-1 pivot multiplier — tested, do not promote

- [x] Benchmark reusing the multiplier in `find_candidate_test`; leave production fast unchanged.

`factor_bunch_kaufman()` computes

```cpp
const double multiplier = system(column, k) * inverse_pivot;
```

for the Schur update, then traverses the same column again and repeats the same multiplication to store `L`. Store `multiplier`
into `system(column, k)` after that column's update instead. This removes one multiplication and the second full column pass for
every 1-by-1 pivot. The multiplication itself and its input values are unchanged.

Measured 2026-08-02 against the saved `equilibrated-fast` run. Both builds used the same canonical Release/native/LTO settings,
CPU 2, one persistent Pybind process, native nanosecond medians, the one-second target, and the same 81 matrices of dimensions
3 through 25. The experiment stored the native `test` method under timing mode `fast` only because the historical SQLite schema
does not admit `test`; build label `fp-s01-test`, its comment, and its binary hash identify the experimental code.

| Panel | Matrices | Median change | Mean change |
|:---|---:|---:|---:|
| all | 81 | +0.3064% | +0.2850% |
| non-circular | 51 | +0.2566% | +0.3038% |
| circular | 30 | +0.3854% | +0.2530% |
| dimensions 3-8 | 20 | +1.0169% | +1.1063% |
| dimensions 9-16 | 28 | +0.4301% | +0.3786% |
| dimensions 17-25 | 33 | -0.0810% | -0.2922% |

All 81 ESS counts matched. The variant won 27 matrices, lost 51, and tied 3. The two long one-sample matrices improved, but that
does not outweigh the balanced result: ID 45 changed by -3.48% and ID 47 by -5.49%, while the overall panel regressed and the
large-matrix median was effectively neutral. Decision: do not copy FP-S01 into production fast. Reset test to the production fast
implementation before starting the next independent experiment.

### FP-S02: Extract the complement only when it is needed — promoted

- [x] Keep the complement mask at function entry and enumerate its bits only after the support probabilities pass.

Before promotion, `find_candidate_fast::find()` filled `non_support_indices` before factorization. Many supports return during
factorization or probability validation, so this work and its stack array are unused. Outside strategies can be visited directly
from the existing complement mask with `ctz` plus clear-lowest-bit iteration. Preserve ascending index order.

Measured 2026-08-02 against the saved `equilibrated-fast` run. Both builds used the same canonical Release/native/LTO settings,
CPU 2, one persistent Pybind process, native nanosecond medians, the one-second target, and the same 81 matrices of dimensions
3 through 25. The experiment stored the native `test` method under timing mode `fast` only because the historical SQLite schema
does not admit `test`; build label `fp-s02-test`, its comment, and its binary hash identify the experimental code.

| Panel | Matrices | Median change | Mean change |
|:---|---:|---:|---:|
| all | 81 | -0.8086% | -0.0231% |
| non-circular | 51 | -0.8638% | +0.0989% |
| circular | 30 | -0.6616% | -0.2305% |
| dimensions 3-8 | 20 | +3.6021% | +3.9998% |
| dimensions 9-16 | 28 | -1.5120% | -1.5646% |
| dimensions 17-25 | 33 | -0.9184% | -1.1533% |

All 81 ESS counts matched. The variant won 60 matrices, lost 18, and tied 3. Its small-panel regression is concentrated in
sub-microsecond and few-microsecond cases; the medium and large panels both improve. The change also deletes one 64-byte stack
array and one eager scan while preserving ascending outside-strategy order. Promoted to fast on 2026-08-02; fast and test retain
the same implementation.

### FP-S03: Convert and equilibrate one symmetric triangle — promoted

- [x] Convert each exact symmetric entry once and mirror it.
- [x] Apply each equilibration product once and mirror it.

Before promotion, `prepare_normalized_double_game()` converted both `A(i,j)` and `A(j,i)`. The equilibration pass likewise computed both
halves. The exact source is symmetric and the scale product is symmetric, so one triangular loop plus the mirrored assignment
removes duplicate conversion and multiplication without changing the value placed in either half.

This is once-per-game work, so expect a modest gain rather than a large end-to-end change.

Measured 2026-08-02 with FP-S02 retained underneath it. The canonical Release/native/LTO, CPU-2, persistent-Pybind,
native-nanosecond-median, one-second protocol and the same 81 dimension-3-through-25 matrices were used again. The new session is
`timing_20260802_fp_s02_s03_triangle`, with build label `fp-s02-s03-test`.

| Panel | Matrices | Combined median vs fast | Combined mean vs fast | FP-S03 median vs FP-S02 | FP-S03 mean vs FP-S02 |
|:---|---:|---:|---:|---:|---:|
| all | 81 | -1.9117% | -1.8817% | -0.8333% | -1.8005% |
| non-circular | 51 | -1.6725% | -1.0522% | -0.6160% | -1.0752% |
| circular | 30 | -2.4419% | -3.2919% | -1.6412% | -3.0335% |
| dimensions 3-8 | 20 | -3.2960% | -1.9803% | -6.4912% | -5.7310% |
| dimensions 9-16 | 28 | -3.2173% | -3.3846% | -1.0477% | -1.8376% |
| dimensions 17-25 | 33 | -1.4652% | -0.5469% | -0.0740% | +0.6132% |

All 81 ESS counts matched. Combined FP-S02+FP-S03 won 62 matrices, lost 14, and tied 5 against production fast. FP-S03 alone
won 60 and lost 21 against the saved FP-S02 run. Its direct paired median contribution is -0.8333%; subtracting the two rounded
cumulative medians suggests about -1.10 percentage points, but medians of per-matrix ratios are not additive and the runs contain
ordinary timing noise. As expected for once-per-game setup work, the clear gain is in small and medium games, while the large-game
median is neutral. Promoted together with FP-S02 on 2026-08-02. The temporary shared-helper flag was removed; conversion now always
uses one triangle, while fast and test retain independent but source-identical equilibration and candidate kernels.

### FP-S04: Fold precision-span extrema into the collection pass

- [ ] Track the minimum and maximum nonzero absolute entry while collecting references for sorting.
- [ ] After sorting, scan only adjacent distinct gaps.
- [ ] Remove the unused one-argument `precision_span_at_least()` and its `include_game_denominator` branch if no caller has appeared.

Sorting is still needed to find the smallest exact nonzero pairwise gap. The separate post-sort pass over every entry for extrema
is not needed. Fast preparation always uses the denominator-free definition, and the more general public overload currently has
no production caller.

The reference vector and its one allocation should remain: it is the simple way to sort without copying arbitrary-precision
integers.

### FP-S05: Delete unused `matrix_dbl` operations

- [ ] Reconfirm and remove the mutable `data()` accessor and `swap_rows()` if they still have no caller.

This is a readability cleanup, not a meaningful speed optimization. It belongs after the functional items and should not grow a
separate abstraction or test suite.

## Exact setup allocation

### FP-A01: Avoid the temporary full FLINT rational matrix

- [ ] Prototype direct denominator clearing in `matrix_int::set_from_fraction_matrix()`.

The current wrapper allocates an `fmpq_mat`, copies every entry from `matrix_frc`, asks FLINT to clear denominators, and destroys
the temporary. The same result can be produced in two direct passes over `matrix_frc`:

1. compute the least common positive denominator;
2. write every numerator multiplied by the corresponding exact scale into `matrix_int`.

This removes one full matrix allocation and copy per game. It is not per-support work, and it is a larger exact-arithmetic change,
so retain it only if the implementation stays short and a safe-mode benchmark shows a repeatable benefit.

## Rounding-changing experiments

The following formulas are mathematically equivalent but change binary64 operation order. Implement only one at a time in
`find_candidate_test`; compare its candidate decisions against safe search before measuring speed.

### FP-E01: Store inverse pivot blocks for the solve

- [ ] Test storing the inverse 1-by-1 or 2-by-2 `D` block after its factor update.

Factorization already computes the inverse 1-by-1 pivot. For a 2-by-2 pivot it also computes the diagonal ratios, determinant
factor, and inverse block scale. `solve_bunch_kaufman()` currently recomputes these quantities and performs several divisions.
After the trailing update no longer needs the original block, the factorization could store the inverse block coefficients and the
solve could use multiplications and additions.

This may save divisions, but it changes representation and rounding in the core numerical method. It is an experiment, not an
automatic cleanup.

### FP-E02: Reuse the reduced-system row offset

- [x] Tested and reverted after changing one known fast decision; no general stability advantage for the retained formula was found.

The current formulas are

```text
r_i   = (d_i * B_mm - B_im) / D_m
H_ij  = B_ij - d_i * B_mj - d_j * B_im + d_i * d_j * B_mm
```

Define once per row

```text
row_offset = d_i * B_mm - B_im
r_i        = row_offset / D_m
H_ij       = B_ij - d_i * B_mj + d_j * row_offset
```

This removes two multiplications and one addition per lower-triangle entry. It is algebraically identical, but reassociation can
change cancellation and therefore rejection decisions.

The test-only implementation passed all 10 ordinary Release tests and all 10 ASan/UBSan tests. An initial 81-matrix comparison was
discarded because its temporary helper omitted each matrix dimension and therefore compared identical parser errors rather than
computations. The canonical native/LTO benchmark used valid inputs and exposed a changed decision before timing completed: matrix
46 retained its exact full-support ESS in `fast` but returned zero candidates in `test`. A focused high-precision diagnostic found
only two of 55 lower-triangle entries changed, each by one binary64 ulp (`1.1102230246251565e-16`). For one entry the factored
formula exactly matched the higher-precision evaluation of the double inputs while the expanded formula was one ulp low; for the
other, the two values lay equally far on opposite sides. The factored construction is therefore not generally less stable. The
reduced matrix has 2-norm condition number about `2.984e10`, so it strongly amplifies the perturbation: one computed probability
changes from `3.0581423239e-9` to `-6.5326723740e-9`, while its exact value is `1e-8`. The retained formula happens to err in the
favorable direction for this matrix. The experiment was reverted and no timing result was retained because adopting it would
knowingly change an existing fast regression result; that is a behavioral trade-off, not a numerical-stability proof for the
retained expression.

The exact `find_candidate_safe` counterpart has no rounding trade-off because its reduced system uses denominator-cleared FLINT
integers. It now stores `A_mm-A_im` directly in the existing right-hand-side entry and reuses that non-owning entry reference in
`A_ij-A_mj+(A_mm-A_im)`. This removes one exact subtraction and one repeated matrix lookup per lower-triangle entry without an
allocation or integer copy. Release and ASan/UBSan passed all 10 tests, and complete exact candidate-output hashes matched for all
76 quick-panel matrices from dimensions 3 through 25.

The first whole-panel timing comparison was discarded: running the entire baseline before the changed build produced an implausible
40-50% difference that disappeared when the baseline was rerun after the processor had warmed. In the controlled warm reverse pair,
the CPU-2 persistent-Pybind median improved by `6.02%` overall, `7.30%` for dimensions 9-16, and `6.76%` for dimensions 17-25.
Non-circular and circular medians improved by `6.95%` and `5.15%`; dimensions 3-8 had a zero median change. The changed build won
63 matrices, tied 11, and lost two. The exact optimization is retained in the working tree.

Hardware counters confirm that this is executed-work reduction rather than frequency bias. Three Perf repeats on non-circular
dimension-23 matrix 53 measured `7.48%` fewer instructions, `8.24%` fewer branches, and `6.58%` fewer cycles. On circular
dimension-25 matrix 70, the corresponding reductions were `6.38%`, `7.30%`, and `5.53%`. One source-level `fmpz_sub` is not one
machine subtraction: FLINT dispatches arbitrary-precision representation, sign, allocation, and normalization cases. Removing that
call from every lower-triangle entry of every visited support eliminates billions of instructions over a complete search.

### FP-E03: Store the reduced scratch in the factorization's traversal order

- [x] Tested and rejected; production and experimental search retain the original row-major lower-triangle storage.

The Bunch-Kaufman loops mostly walk fixed columns while `matrix_dbl` is row-major. Storing the logical lower triangle transposed
would make those hot walks contiguous. Keep this local to the fast/test implementation; do not create a new general matrix class
for one kernel.

Required check: cover both 1-by-1 and 2-by-2 pivot paths, compare exact candidate contracts, and use a balanced timing panel. Keep
it only for a clear gain across matrix types and dimensions.

The test-only transposed accessor passed all 10 Release tests and all 10 ASan/UBSan tests. Its initial claimed 81-matrix complete
comparison used the same invalid dimensionless helper described under FP-E02 and is discarded. The valid CPU-2 persistent-Pybind
benchmark checked matching ESS counts and fallback classifications while covering 76 dimension-3-to-25 matrices in both call
orders. Transposition was slower: median changes were `+1.0605%` with `fast` called first and `+1.2794%` with `test` called first.
For dimensions 17-25 the corresponding medians were `+1.9864%` and `+2.1579%`; both circular and non-circular subsets regressed.
On dimension-23 matrix 53, three 100-run Perf repeats measured about `+4.55%` instructions, `+4.22%` branches, and `+2.37%`
cycles, with essentially unchanged branch misses. A complete candidate-contract comparison was unnecessary after this decisive
performance rejection.

The scratch matrices are at most 24-by-24 in this panel and fit comfortably in L1 cache. Transposition therefore saves little
cache-miss cost, while it reverses the efficient row-major traversal used to construct the reduced system and produces more
instructions. An earlier test-only row-oriented update without transposition was also rejected: its forward/reverse median
regressions were `+0.7609%` and `+0.9746%`. No further storage-order variant is justified without a materially different kernel.

### FP-E04: Fuse the forward right-hand-side transformation into factorization

- [ ] Defer unless profiling still identifies the triangular solve as a material cost after the simpler items.

The factorization could update the single right-hand side while applying each pivot, leaving only the backward solve. This removes
one triangular traversal but couples factorization and solve and makes the code harder to verify. The expected gain is smaller than
the first three experiments, so it is deliberately last.

## Larger target after the local kernel

### FP-G01: Support generation and forbidden-support checks

- [ ] Reprofile only after the local checklist above is resolved.
- [ ] Benchmark the unused alternative circular generator or delete it; do not retain an untested second implementation.

Profiles show that support generation can dominate complete runs even when the candidate kernel is efficient. It accounted for
about 33% on the reviewed non-circular dimension-20 case and about 48% on the circular dimension-25 case. This is the plausible
place for a large whole-program gain, but it is also algorithmic work and should not be mixed into the local solver cleanups.

## Reviewed and deliberately retained

- [x] Fast and safe both build a reduced system. Safe must reconstruct and prove every surviving support exactly; this is not
  redundant work that can be deleted.
- [x] Fast and safe both extract support indices for surviving supports. Passing a tiny stack array through the orchestration would
  add plumbing for little expected gain.
- [x] Probability validation remains separate from payoff validation. It preserves the early exit on a negative probability.
- [x] Per-support numerical scratch (`solution`, scale ratios, pivots, and index arrays) is fixed-size stack storage; there is no
  hidden per-support heap allocation in the fast candidate routine.
- [x] `reduced_system_` is reallocated only when support cardinality changes, not once per support.
- [x] The precision-span vector is allocated once per game and is needed for sorting exact entry references.
- [x] The unused upper half of `reduced_system_` remains zero-initialized. Avoiding it would require packed or uninitialized custom
  storage for only a few cardinality transitions.
- [x] Exact outside-payoff validation after the fast check is intentional proof, not duplicated work.
- [x] The rational game owned by `fracessa` remains owned storage. Borrowing or moving it would complicate lifetime and API rules
  for a once-per-game copy.
- [x] Both fast and test helpers may remain constructed as members. They contain fixed arrays and default matrices, not recurring
  heap allocations; a variant or optional wrapper would add complexity.
- [x] `game_scales_` is value-initialized and then overwritten during equilibration. Avoiding that trivial fixed-array
  initialization is not worth custom initialization rules.
- [x] Extracting `support_count` and asserting that it equals the supplied `support_size` documents a useful invariant. Trusting
  only one value would not remove the index extraction that candidate construction needs.
- [x] Fast and test remain separate source files while test is the independent experiment path. Do not abstract their shared code
  until that independence is no longer required.
- [x] The common denominator is cleared once, equilibration runs once per game, exact inertia is reused, and exact candidate-vector
  materialization is deferred until a successful candidate needs output or logging.

## Evidence already collected

Correctness checks at the reviewed tree:

- CTest: 10 of 10 passed.
- Fixed-seed random comparison: 7,750 integer games matched between fast and safe candidate results.
- Constructed full-support comparison: 2,441 positive candidates, including zero-diagonal indefinite systems, produced no fast
  miss.
- Retained database comparison: all 90 matrices matched the exact candidate contracts.
- Earlier promotion checks: ASan/UBSan 10 of 10 and Python tests 56 of 56.

These checks provide useful evidence but do not turn the remaining binary64 inequality tests into a proof.

One diagnostic profile divided time as follows; these are profile shares, not benchmark medians:

| Matrix | Shape | Support generation | Fast body | Factorization | Solve |
|---:|:---|---:|---:|---:|---:|
| 53 | n=23, non-circular | 0.7% | 25.1% | 46.1% | 15.6% |
| 61 | n=20, non-circular | 33.1% | 26.0% | 25.9% | 12.2% |
| 70 | n=25, circular | 48.3% | 10.0% | 28.2% | 8.6% |

Earlier setup experiments also found:

- removing the precision-span decision changed the median by about +0.1%, effectively noise;
- equilibration was about 1.13% slower overall on the balanced 81-matrix panel, with mixed effects by size and symmetry;
- equilibration is retained for numerical behavior, not because it is a universal speed improvement.

## Point-by-point completion rule

For every open item:

1. make the smallest readable change in `find_candidate_test` when rounding can change, otherwise use the shared production
   location directly;
2. run the focused unit or black-box regression;
3. compare candidate and ESS results with safe mode;
4. for a speed claim, measure native C++ elapsed nanoseconds and report the per-matrix median change on a balanced panel;
5. copy a successful numerical experiment into fast only after explicit approval;
6. mark the item complete here and record the retained evidence or rejection reason.
