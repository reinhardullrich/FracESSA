# Support Generator Design Handover

Status: implemented on `main`; the certified-filter worktree uses the same
generator architecture unchanged.

This is a durable design handover, not a transcript of the implementation
session. The detailed algorithms, proofs, counterexamples, and benchmarks remain
in [`../plans/SUPPORT_GENERATORS.md`](../plans/SUPPORT_GENERATORS.md).

## Branch Context

The current worktree contains three connected phases:

1. the temporary numerical candidate-search experiments;
2. the support-frontier and circular-bracelet experiments that established the
   replacement design;
3. the production one-support-at-a-time generators and compressed circular
   candidate output.

Production now requires `fast` or `safe` before the matrix. Fast uses the
historical raw-double heuristic before exact confirmation; safe starts with
exact candidate solving. There is no default method.

## Agreed Production Architecture

FracESSA no longer creates the complete set of nonempty supports or even a
complete cardinality layer before solving. A support is one `uint64_t` mask.
Each concrete generator owns the matrix dimension, the cardinality sweep, its
recursive state, and the forbidden-support rules discovered during analysis.

- `NonCircularSupportGenerator` uses fixed-cardinality binary DFS.
- Production `CircularSupportGeneratorV3` uses direct fixed-density bracelet recursion, placing one selected bit per recursive
  step and rejecting reflected duplicates during construction.
- Both emit cardinality ascending first and numeric mask ascending second.
- Both pass exactly one support and its cardinality to the analyzer at a time.
- Both activate newly forbidden supports only when moving to the next
  cardinality. Distinct supports in the same cardinality therefore cannot prune
  one another.

The full-support option remains a special first attempt. If it fails to find an
ESS and normal generation follows, the generator's full-support emission is
ignored because that mask was already analyzed.

## Why Generation Uses A Callback

The production interface is deliberately callback-driven:

```cpp
generator.generate([&](bitset64 support, size_t support_size) {
    // analyze this one support, register an exact candidate, then return
});
```

`generate()` is recursive. When it finds a support, it calls the consumer
synchronously while the DFS or bracelet call stack is still alive. After the callback
returns, recursion continues from exactly the saved position. The callback is a
template defined in the header, so the compiler can inline it; there is no
`std::function`, virtual dispatch, or allocation per support.

A pull interface such as `next_support()` would be logically valid and would
make the caller superficially simpler. The difficulty would move inside the
generator: every return would otherwise destroy the recursive call stack, so
the implementation would need a manual stack, an explicit state machine, or
C++ coroutine machinery. That adds state loads and branches and is not expected
to outperform the inline callback. The callback was therefore retained because
it is the simplest implementation of these recursive algorithms, not because a
pull interface would be incorrect.

The callback receives both the support mask and support size because the
generator owns the cardinality loop. When an exact candidate is found,
`add_forbidden()` is called synchronously before the callback returns. This is
especially important for the test-only V2 generator, whose returned multiplier
belongs to the support most recently emitted to that callback.

## Candidate And Pruning Lifecycle

`fracessa::analyze_support()` performs the optional fast search, safe exact
equilibrium solve, candidate population, and stability classification.
If it reports an exact candidate, the caller registers `candidate_.support` with
the active generator and then calls `finalize_candidate()`.

Every exact equilibrium support becomes forbidden for larger cardinalities,
whether or not stability classifies that candidate as an ESS. A result accepted
only by exact rational analysis may become a pruning rule; neither numerical
procedure creates one by itself.

The explicit sequence remains visible in each of the three call paths: full
support, circular generation, and non-circular generation. Moving generator
registration into `analyze_support()` was considered and rejected. The two
generators have different `add_forbidden()` results, while the full-support path
has no generator at all; hiding those differences would require templates,
overloads, inheritance, or another callback for only a few direct lines.

`finalize_candidate()` remains separate because it performs work that follows
successful analysis and generator registration: assigning the representative
ID, applying the optional orbit multiplier to the ESS count, optionally storing
the candidate, and logging it.

## Why There Is No Generator Base Class

The two production generators share a tiny compile-time calling convention but
not an implementation. Their recursive state and pruning representations are
different, and circular `add_forbidden()` also produces an orbit multiplier.
Two explicit matrix-type branches are shorter and clearer than a virtual base,
adapter, or type-erased wrapper. They also keep the per-support hot path free of
runtime dispatch.

The same reasoning applies to the three retained circular implementations: keeping historical and experimental algorithms does
not justify a runtime hierarchy.

## Circular Output And Multiplier

By default, circular matrices store only the smallest integer support in each complete rotation/reflection orbit. The candidate's
non-null `multiplier` is the number of distinct raw supports represented by that row; it can be one. Non-circular candidates use a
null multiplier.

Every circular matrix automatically detects exact affine index multipliers once from the rational game. The helper filters V3
output before `analyze_support()`. For an exact candidate it canonicalizes each multiplier image back to a bracelet, deduplicates
those bracelets, and calls V3's existing `add_forbidden()` for each one, so the forbidden family remains closed under the verified
affine group. It also reconstructs one candidate row for every distinct dihedral image by exactly permuting the solved vector,
support, and extended support. Each row's multiplier therefore retains its universal meaning—only distinct rotations and
reflections, at most $2n$—and no matrix-specific symmetry information is required to expand the output. V3's recursion is
unchanged, the helper disengages when it finds no extra multiplier, and it never runs for non-circular matrices.

Candidate and ESS totals remain mathematical totals, so circular ESS rows are
counted with their multiplier. The nullable field replaced `shift_reference`
through the CLI, pybind boundary, Python APIs, Arrow and Parquet sinks, CSV
verification baseline, and SQLite fixtures.

Compressed output and pruning are separate concerns. The production circular
generator still registers every distinct rotation and reflection of an exact
candidate support as compact forbidden masks. One canonical mask alone is not
enough because canonicalization does not preserve subset relationships.

## Retained Circular Generators

`CircularSupportGenerator` is V1: fixed-content FKM necklace recursion followed by reflection reduction. It remains available for
regression and performance comparisons but is no longer wired into production.

`CircularSupportGeneratorV3` is the production path. It is the binary specialization of Karim-Alamgir-Husnine `BraceFD`: each
recursive step places one selected bit and skips the zero run before it, while the paper's reversal state rejects mirror branches
during construction. V3 retains V1's callback, expanded forbidden-orbit masks, and independently calculated orbit multiplier.

### Experimental `CircularSupportGeneratorV2`

`cpp/include/fracessa/supports.hpp` also contains
`CircularSupportGeneratorV2`. It is explicitly test-only and is not wired into
production.

V2 stores one forbidden support per bracelet orbit. During recursion it uses
two 64-bit alignment masks to test all rotations and reflections in parallel.
The arbitrary-width left rotation needed by that algorithm lives with the other
bitset operations as `bs64::rot_left()`.

V2 does not enumerate rotations merely to calculate the multiplier. The FKM
leaf already supplies the necklace period, and reflection comparison determines
whether the reflected necklace is the same rotation class. The multiplier is
therefore `period` or `2 * period` and is cached for the synchronous
`add_forbidden()` call.

The compact representation is correct but slower. The focused 2026-08-02 fast-mode benchmark verified identical results on all 81
quick-test matrices and timed the 33 circular cases. V2 was slower by 41.48% at the median and 91.40% by geometric mean; matrix 34
was 6.824 times slower. V1 was faster on 28 cases, equal on three, and V2 won only two sub-microsecond dimension-2/3 measurements.
The expanded forbidden-orbit representation therefore remains in V3; only V2's compact pruning representation was rejected. See
`experiments/circular_support_v2_2026-08-02/README.md` for the complete method and table.

## Validation Record

The production replacement was checked order-independently against all active
matrices so deliberate candidate-ID and ordering changes could not hide a
mathematical difference. Release and ASan/UBSan builds each passed all 62 CTests,
and the wrapper and regenerated baseline/database checks passed.

Independent generator checks also compared both circular implementations with
brute-force enumeration and pruning for manageable dimensions, verified every
reported orbit multiplier, exercised arbitrary forbidden rules, and covered the
dimension-63 boundary. The temporary V2 checker passed in Release and under
ASan/UBSan.

The production V3 promotion compared V1 and V3 generation through dimension 24, covered V3's dimension-63 boundary, and produced
byte-identical fast and safe candidate output on all 81 quick-test matrices. Release and ASan/UBSan CTests plus all 63 Python tests
passed. Two CPU-2 persistent-Pybind panels found conservative reverse-order median runtime reductions of 19.90% overall and 34.70%
for dimensions 19 and above. See `experiments/direct_bracelet_generation_2026-07-29/PRODUCTION_V3_COMPARISON_2026-08-03.md`.

## Remaining Work

- Do not revisit V1 or V2 for production unless a new end-to-end measurement beats V3 while preserving the same pruning contract.
- Do not replace the callback with `next_support()` without a measured reason
  that justifies a resumable generator implementation.
- Do not funnel generator registration into `analyze_support()` unless the
  surrounding architecture changes enough to remove, rather than hide, the
  current type differences.
- After this branch is merged and updated `main` is verified, remove the
  `choice-one` worktree and delete the local and remote feature branch.
