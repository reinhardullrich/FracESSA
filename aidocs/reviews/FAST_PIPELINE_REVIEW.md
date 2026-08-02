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
was added merely to inject an impossible internal state. The complete native and Python suites pass with the guards active. The
MSVC Release configuration remains open as FP-C03 below.

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

### FP-C03: MSVC Release disables reliable special-value behavior

- [ ] Remove `/fp:fast` from the MSVC Release options or replace it with `/fp:precise`.

The fast and test paths now rely on standard NaN and infinity propagation plus `std::isfinite()` to select exact fallback.
`cpp/CMakeLists.txt` still applies `/fp:fast` to every MSVC Release source. Microsoft documents that special values may not
propagate or behave according to IEEE 754 under that option. The Linux and macOS configurations do not enable an equivalent
fast-math mode, so this is a Windows Release correctness gap rather than a defect in the fallback branches themselves.

Required outcome: compile MSVC Release with `/fp:precise`, run the CLI counterexample regressions on Windows, and only then treat
the non-finite fallback as a cross-platform guarantee.

## Bit-identical or structural simplifications

These changes do not intentionally alter the arithmetic. They should still be tested because they touch the hot solver.

### FP-S01: Reuse the 1-by-1 pivot multiplier

- [ ] Remove the second scaling pass over the accepted 1-by-1 pivot column.

`factor_bunch_kaufman()` computes

```cpp
const double multiplier = system(column, k) * inverse_pivot;
```

for the Schur update, then traverses the same column again and repeats the same multiplication to store `L`. Store `multiplier`
into `system(column, k)` after that column's update instead. This removes one multiplication and the second full column pass for
every 1-by-1 pivot. The multiplication itself and its input values are unchanged.

### FP-S02: Extract the complement only when it is needed

- [ ] Keep the complement mask at function entry and enumerate its bits only after the support probabilities pass.

`find_candidate_fast::find()` currently fills `non_support_indices` before factorization. Many supports return during
factorization or probability validation, so this work and its stack array are unused. Outside strategies can be visited directly
from the existing complement mask with `ctz` plus clear-lowest-bit iteration. Preserve ascending index order.

### FP-S03: Convert and equilibrate one symmetric triangle

- [ ] Convert each exact symmetric entry once and mirror it.
- [ ] Apply each equilibration product once and mirror it.

`prepare_normalized_double_game()` currently converts both `A(i,j)` and `A(j,i)`. The equilibration pass likewise computes both
halves. The exact source is symmetric and the scale product is symmetric, so one triangular loop plus the mirrored assignment
removes duplicate conversion and multiplication without changing the value placed in either half.

This is once-per-game work, so expect a modest gain rather than a large end-to-end change.

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

- [ ] Test factoring the repeated row term out of reduced-system construction.

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

### FP-E03: Store the reduced scratch in the factorization's traversal order

- [ ] Prototype a local transposed accessor over the existing row-major `matrix_dbl` storage.

The Bunch-Kaufman loops mostly walk fixed columns while `matrix_dbl` is row-major. Storing the logical lower triangle transposed
would make those hot walks contiguous. Keep this local to the fast/test implementation; do not create a new general matrix class
for one kernel.

Required check: cover both 1-by-1 and 2-by-2 pivot paths, compare exact candidate contracts, and use a balanced timing panel. Keep
it only for a clear gain across matrix types and dimensions.

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
