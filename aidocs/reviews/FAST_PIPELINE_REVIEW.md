# Fast Candidate Pipeline Review

Last verified: 2026-08-03

Scope: the complete binary64 candidate path from integer game preparation through `find_candidate_fast`, plus the exact handoff
that protects candidates which fast cannot reject. This is an actionable checklist, not a transcript of every experiment.
Correctness comes first; speed changes are retained only when they help a balanced panel of small and large, circular and
non-circular games while keeping the code readable.

Status:

- `[ ]` is open.
- `[x]` is completed, rejected after measurement, or deliberately retained.

## Overall conclusion

The current algebra is coherent: it clears denominators once, applies one common binary scale, equilibrates the complete game,
builds the reduced symmetric candidate system, factors it with Bunch-Kaufman LDLT, recovers probabilities and payoff, and hands
every surviving support to exact arithmetic.

The reviewed factorization and solve paths appear internally correct, including 1-by-1 and 2-by-2 pivots and permutations. The
retained regression matrices and random comparisons agree with safe search.

Fast remains a heuristic prefilter, not a proof. The precision-span cutoff and small-pivot fallback recover the known
counterexamples, but the remaining binary64 probability and outside-payoff rejections have no general error bound. A false
binary64 rejection can therefore still lose an exact candidate. The boundary and concrete examples are documented in
[`FAST_CANDIDATE_FALSE_REJECTION.md`](../correctness/FAST_CANDIDATE_FALSE_REJECTION.md).

## Open findings

| Priority | ID | Work | Expected impact |
|:---:|:---|:---|:---|
| 1 | FP-E01 | Store inverse pivot blocks for the double solve | The only remaining local hot-path experiment with a plausible measurable gain |
| 2 | FP-S04 | Fold precision-span extrema into its collection pass | Small once-per-game simplification |
| 3 | FP-A01 | Avoid the temporary full FLINT rational matrix | Safe setup only; retain only if a benchmark justifies the extra exact code |
| 4 | FP-S05 | Delete unused `matrix_dbl` operations | Readability only |
| 5 | FP-E04 | Fuse the right-hand-side transformation into factorization | Deferred: more coupling for a smaller expected gain |

Prototype FP-E01 and FP-E04 only in `find_candidate_test` because they change binary64 representation or operation order. Promote
either only after exact-result comparison and a balanced benchmark. FP-S04, FP-A01, and FP-S05 do not justify another abstraction.

### FP-E01: Store inverse pivot blocks for the solve

- [ ] Test storing the inverse 1-by-1 or 2-by-2 `D` block after its factor update.

Factorization already computes the inverse 1-by-1 pivot. For a 2-by-2 pivot it also computes the diagonal ratios, determinant
factor, and inverse block scale. `solve_bunch_kaufman()` recomputes these values and performs several divisions. After the trailing
update no longer needs the original block, factorization could store the inverse block coefficients and solve could use
multiplications and additions.

This may save divisions, but it changes representation and rounding in the numerical core. It is an experiment, not an automatic
cleanup.

### FP-S04: Fold precision-span extrema into the collection pass

- [ ] Track the minimum and maximum nonzero absolute entry while collecting references for sorting.
- [ ] After sorting, scan only adjacent distinct gaps.
- [ ] Remove the unused one-argument `precision_span_at_least()` and its `include_game_denominator` branch.

The 2026-08-03 recheck confirms that the one-argument overload still has no caller. Sorting remains necessary to find the smallest
exact nonzero pairwise gap, but the separate post-sort scan over all entries for extrema is unnecessary. Active fast preparation
always uses the denominator-free definition.

The reference vector and its one allocation should remain: sorting references is simpler than copying arbitrary-precision
integers. This is once-per-game work, so the likely end-to-end gain is small.

### FP-A01: Avoid the temporary full FLINT rational matrix

- [ ] Prototype direct denominator clearing in `matrix_int::set_from_fraction_matrix()`.

The 2026-08-03 recheck confirms that the wrapper still allocates an `fmpq_mat`, copies every entry from `matrix_frc`, asks FLINT to
clear denominators, and destroys the temporary. A direct version needs two passes: compute the least common positive denominator,
then write each numerator multiplied by its exact scale.

This removes one full matrix allocation and copy per game, not per support. It is a larger exact-arithmetic change and may lose to
FLINT's optimized bulk routine. Retain it only if the implementation stays short and a safe benchmark shows a repeatable gain.

### FP-S05: Delete unused `matrix_dbl` operations

- [ ] Remove mutable `matrix_dbl::data()` and `matrix_dbl::swap_rows()`.

The 2026-08-03 caller recheck found no use of either operation. This is readability cleanup, not a speed optimization; it needs no
new test or abstraction.

### FP-E04: Fuse the forward right-hand-side transformation into factorization

- [ ] Defer unless profiling identifies the triangular solve as material after FP-E01.

Factorization could update the single right-hand side while applying each pivot, leaving only the backward solve. This removes one
triangular traversal but couples factorization and solve and makes verification harder. Its expected gain is smaller than FP-E01.

## Completed, rejected, and deliberately retained

### Correctness

#### FP-C01: Non-finite completed values fall back

- [x] Fast and test return candidate-possible after a non-finite reference probability, payoff, threshold, or outside payoff.

No retained exact input reliably reaches these guards through the existing precision-span, equilibration, pivot, and coordinate
checks, so no artificial production hook was added merely to inject an internal state.

#### FP-C02: Historical false-rejection regressions are covered

- [x] CLI tests compare complete `fast` and `test` candidate output with `safe` for the small-pivot, positive-probability, and
  outside-payoff counterexamples and check their ESS counts of 1, 1, and 2.
- [x] The stale comment claiming that current fast rejects the precision-span example was corrected.

#### FP-C03: MSVC Release preserves special values

- [x] MSVC Release uses `/fp:precise`, not `/fp:fast`.

Linux and macOS already avoid equivalent fast-math behavior. The next Windows release still needs its normal platform build and
counterexample run before publication; that deferred platform execution does not make the source finding open.

### Structural and numerical experiments

#### FP-S01: Reuse the 1-by-1 pivot multiplier — rejected

- [x] Tested in `find_candidate_test`; production fast remains unchanged.

The balanced historical 81-matrix panel regressed by 0.3064% at the median. It won 27, lost 51, and tied 3. The removed
multiplication did not compensate for the changed storage/traversal pattern.

#### FP-S02: Extract the complement only when needed — promoted

- [x] Fast and test keep the complement mask at entry and enumerate outside bits only after probability checks pass.

The historical 81-matrix panel improved by 0.8086% at the median, winning 60, losing 18, and tying 3. The change also removed one
eager scan and one 64-byte stack array.

#### FP-S03: Convert and equilibrate one symmetric triangle — promoted

- [x] Exact conversion and equilibration process one triangle and mirror it.

Added on top of FP-S02, its direct paired median contribution was 0.8333%; combined FP-S02 and FP-S03 improved 1.9117% against the
original fast panel. All 81 ESS counts matched.

#### FP-E02: Reuse the reduced-system row offset

- [x] The binary64 reassociation was tested and reverted.

Matrix 46 changed two reduced entries by one binary64 ulp. Its ill-conditioned system amplified that difference enough to change
a probability from positive to negative and lose a known ESS. Neither algebraic form is generally more stable, so production keeps
the form that preserves the retained regression.

- [x] Exact arithmetic reuses the row offset safely.

The exact solver first reused `A_mm-A_im` directly. The controlled warm panel improved 6.02% at the median, and Perf confirmed
fewer instructions, branches, and cycles. It subsequently added a dense lazy cache for exact reduced entries; the balanced
seven-matrix speed panel improved 4.99% at the median, while all 81 quick matrices matched for correctness.

The fraction-free solve now starts back substitution with the absolute final pivot, directly returning a positive common
denominator and avoiding a whole-vector negation. The full 81-matrix quick benchmark matched every ESS count: 64 improved, 11 tied,
and 6 regressed, with a 1.42% median improvement. Excluding dimension 2 under the benchmark rule gives 64 improved, 8 tied, and 5
regressed across 77 matrices, with 1.50% median and 1.37% mean improvements.

#### FP-E03: Store reduced scratch in factorization traversal order — rejected

- [x] Tested transposed storage and a row-oriented update; both were slower.

The valid forward/reverse transposed benchmarks regressed by 1.0605% and 1.2794% at the median; dimensions 17-25 regressed by
1.9864% and 2.1579%. The earlier row-oriented update regressed by 0.7609% and 0.9746%. The small scratch matrices fit in L1, so
transposition saves little and makes construction less efficient.

#### FP-G01: Direct circular support generation — completed with V3

- [x] V2 was benchmarked against V1 and retained test-only.
- [x] Direct fixed-density bracelet generator V3 was implemented, verified, benchmarked, and made production.
- [x] V1 and V2 remain test-only at the user's explicit request.

V3 and V1 emitted identical ordered representatives through dimension 24, matched a dynamic forbidden-support run, and passed the
dimension-63 boundary. Complete quick-set candidate output matched in both fast and safe modes, with all C++, sanitizer, and Python
tests passing.

On the 31 circular quick matrices of dimension at least 3, V3 improved the median by 23.68% when V1 ran first and 19.90% in reverse
order. The conservative reverse order had 30 wins and one tie; dimensions 19 and above improved 34.70% at the median. Full evidence
is in [`PRODUCTION_V3_COMPARISON_2026-08-03.md`](../../experiments/direct_bracelet_generation_2026-07-29/PRODUCTION_V3_COMPARISON_2026-08-03.md).

### Reviewed and deliberately retained

- [x] Fast and safe both build a reduced system. Safe must reconstruct and prove every surviving support exactly.
- [x] Both paths extract support indices. Passing a tiny stack array through orchestration would add plumbing for little gain.
- [x] Probability validation remains before payoff validation, preserving the early exit on a negative probability.
- [x] Per-support fast scratch uses fixed-size stack storage; there is no hidden per-support heap allocation.
- [x] `reduced_system_` is reallocated only when support cardinality changes.
- [x] The precision-span reference vector is allocated once per game and is needed for sorting.
- [x] The unused upper half of `reduced_system_` stays zero-initialized; packed custom storage is not justified.
- [x] Exact outside-payoff validation after fast is proof work, not duplicated work.
- [x] The rational game owned by `fracessa` remains owned storage; borrowing it would complicate lifetime rules for one copy.
- [x] Fast and test helpers remain members. They contain fixed arrays and default matrices, not recurring allocations.
- [x] `game_scales_` value initialization remains; custom initialization rules are not justified.
- [x] Extracting and checking `support_count` documents an invariant and does not duplicate removable work.
- [x] Fast and test remain separate source files so test can be changed independently.
- [x] The common denominator is cleared once, equilibration runs once per game, exact inertia is reused, and public exact vectors are
  materialized only for successful candidates that need output or logging.

## Retained evidence

- Current native C++/CLI suite: 10 of 10 passed.
- Current Python suite: 63 of 63 passed.
- Current sanitizer suite: 10 of 10 passed.
- Fixed-seed historical comparison: 7,750 integer games matched between fast and safe candidate results.
- Constructed historical comparison: 2,441 positive full-support candidates produced no fast miss.
- The latest 81-matrix exact quick comparison and timing run matched every expected ESS count.

These checks are evidence, not a proof for the remaining binary64 inequality decisions.

## Completion rule for an open item

1. Make the smallest readable change in `find_candidate_test` when rounding can change; otherwise use the owning production
   location directly.
2. Run the focused unit or black-box regression.
3. Compare candidate and ESS results with safe mode.
4. For a speed claim, report per-matrix native-nanosecond changes on a balanced panel.
5. Promote a successful numerical experiment only after explicit approval.
6. Move the item into the completed section with its retained or rejected evidence.
