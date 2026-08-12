# Fast Candidate Pipeline Experiments

Status: historical decision record completed on 2026-08-03. This file contains no open tasks and is not the source of current
architecture or validation facts. Use `../PROJECT.md` for the live pipeline and
`../correctness/FAST_CANDIDATE_FALSE_REJECTION.md` for its known correctness boundary.

This record preserves completed and rejected binary64 candidate-path experiments whose measurements prevent unproductive work from
being repeated without new evidence.

## Conclusion At Completion

The production algebra at completion was coherent: it cleared denominators once, applied one common binary scale, equilibrated the complete game,
builds the reduced symmetric candidate system, factors it with Bunch-Kaufman LDLT, recovers probabilities and payoff, and hands
every surviving support to exact arithmetic.

The reviewed factorization and solve paths appear internally correct, including 1-by-1 and 2-by-2 pivots and permutations. The
retained regression matrices and random comparisons agree with safe search.

Fast remains a heuristic prefilter, not a proof. The precision-span cutoff and small-pivot fallback recover the known
counterexamples, but the remaining binary64 probability and outside-payoff rejections have no general error bound. A false
binary64 rejection can therefore still lose an exact candidate. The boundary and concrete examples are documented in
[`FAST_CANDIDATE_FALSE_REJECTION.md`](../correctness/FAST_CANDIDATE_FALSE_REJECTION.md).

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

#### FP-E04: Fuse the forward right-hand-side transformation into factorization — promoted

- [x] Fast and test apply each pivot and forward right-hand-side update during factorization, leaving only the backward solve.

This removes one triangular traversal but couples factorization and solve. All 10 Release tests, all 63 Python tests, all 10
ASan/UBSan tests, and complete candidate output for all 81 quick matrices matched production fast before promotion.

The requested single-call benchmark improved 61, tied 1, and regressed 15 of the 77 performance-eligible matrices. Its median,
mean, and summed-runtime changes were `-3.61%`, `-2.74%`, and `-3.02%`. The median changes were `-4.01%` for dimensions 9--16 and
`-4.40%` for dimensions 17 and above. The isolated single-run regressions include noisy short cases and matrix 31 at `+18.58%`;
the broad larger-matrix result justified promotion to production fast.

#### FP-A01: Avoid the temporary full FLINT rational matrix — rejected

- [x] Tested direct two-pass denominator clearing; retained FLINT's bulk matrix conversion.

The direct version needed no second full matrix: it first accumulated the least common positive denominator, then wrote each
scaled numerator directly into the existing `matrix_int` destination. It passed all 10 Release tests, all 63 Python tests, all 10
ASan/UBSan tests, and matched every expected ESS count in the complete 81-matrix quick benchmark.

It was nevertheless slower. Across the 77 performance-eligible matrices of dimension at least 3, it improved 14, tied 2, and
regressed 61. The median, mean, and summed-runtime changes were `+1.18%`, `+0.93%`, and `+1.15%`. The removed allocation occurs only
once per game, while the replacement performs entry-by-entry LCM, exact division, and multiplication. Retain the temporary
`fmpq_mat` because FLINT's optimized bulk conversion is measurably faster end to end.

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

- [x] V2 was benchmarked against V1 and retained as a historical failed experiment.
- [x] Direct fixed-density bracelet generator V3 was implemented, verified, benchmarked, and made production.
- [x] V1 remains for comparisons. V2 was slower than V1 and V3, never entered production, and remains only to record the rejected
  design.

V3 and V1 emitted identical ordered representatives through dimension 24, matched a dynamic forbidden-support run, and passed the
dimension-63 boundary. Complete quick-set candidate output matched in both fast and safe modes, with all C++, sanitizer, and Python
tests passing.

On the 31 circular quick matrices of dimension at least 3, V3 improved the median by 23.68% when V1 ran first and 19.90% in reverse
order. The conservative reverse order had 30 wins and one tie; dimensions 19 and above improved 34.70% at the median. The retained
tracked summary is in [`SUPPORT_GENERATORS.md`](../architecture/SUPPORT_GENERATORS.md); the full raw report remains local under
`experiments/direct_bracelet_generation_2026-07-29/PRODUCTION_V3_COMPARISON_2026-08-03.md`.

#### FP-E01: Reuse the inverse 2-by-2 pivot-block scale — rejected

- [x] Tested only in `find_candidate_test`; production fast remained unchanged and test was restored to match it.

After FP-E04 fused the forward right-hand-side solve into factorization, the original proposal no longer had a separate solve in
which to store inverse pivot blocks. The remaining experiment reused factorization's existing `inverse_block_scale` for the
2-by-2 right-hand side, replacing four divisions with two multiplications. The forms are algebraically equal but round differently
in binary64.

All 81 quick matrices produced identical complete fast and test candidate output and their expected ESS counts, and all 10 Release
C++ tests passed. In the requested one-call CPU-2 comparison, the 77 performance-eligible matrices of dimension at least 3 split
almost exactly evenly: 39 improved and 38 regressed. The median change was `-0.0006%`, effectively zero, while summed runtime
regressed by `0.1201%`. The arithmetic mean was not informative because the test call for matrix 2207 paid the process's one-time
exact-fallback initialization. The experiment showed no measurable gain and does not justify changing the numerical core.

#### FP-S04: Fold precision-span extrema into the collection pass — rejected

- [x] Tested in the shared precision-span helper and restored the original separate scan.

The experiment collected the minimum and maximum nonzero absolute entries while gathering references for sorting, then scanned
only adjacent distinct gaps after the sort. It also removed the unused denominator-aware overload and branch. This eliminated one
short cached traversal per game, but interleaved FLINT comparisons with the otherwise simple reference-collection loop.

All 81 quick matrices produced identical complete candidate output, fallback classifications, and expected ESS counts, and all 10
Release C++ tests passed. The canonical CPU-2 comparison used each stored calibration to target 0.5 seconds per build on the 77
performance-eligible matrices: 17 improved, 15 tied, and 45 regressed. The median, mean, and summed-runtime changes were
`+0.2177%`, `+0.0425%`, and `+0.4995%`. The all-zero dimension-50 matrix improved by 16.07% because the experiment returned before
sorting its zero references, but that special-case gain did not justify the broad regression or changing the clearer general path.

#### FP-S05: Delete unused `matrix_dbl` operations — completed

- [x] Removed mutable `matrix_dbl::data()` and `matrix_dbl::swap_rows()` after a complete active-caller audit found no use.

The similarly named calls elsewhere belong to `std::vector`, the exact integer and fraction types, or the retained rational LU
class. Removing the two dead operations also removed the unused `<utility>` include and corrected the class comment. This is a
readability reduction with no intended runtime effect; both Release configurations built successfully and all 10 C++ tests passed.

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

## Validation At Completion

- Native C++/CLI suite: 10 of 10 passed.
- Python suite: 63 of 63 passed.
- Sanitizer suite: 10 of 10 passed.
- Fixed-seed historical comparison: 7,750 integer games matched between fast and safe candidate results.
- Constructed historical comparison: 2,441 positive full-support candidates produced no fast miss.
- The latest 81-matrix exact quick comparison and timing run matched every expected ESS count.

These checks are evidence, not a proof for the remaining binary64 inequality decisions.
