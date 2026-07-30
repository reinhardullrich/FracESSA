# Support Generator Design Handover

Status: implemented on `main`; the certified-filter worktree uses the same
generator architecture unchanged.

This is a durable design handover, not a transcript of the implementation
session. The detailed algorithms, proofs, counterexamples, and benchmarks remain
in [`../plans/SUPPORT_GENERATORS.md`](../plans/SUPPORT_GENERATORS.md).

## Branch Context

The current worktree contains four connected phases:

1. the temporary explicitly unsafe numerical candidate filter and the exact
   fallback mode;
2. the support-frontier and circular-bracelet experiments that established the
   replacement design;
3. the production one-support-at-a-time generators and compressed circular
   candidate output;
4. the rigorously one-sided Choice 1 candidate filter.

The strict bounded-error "Choice 1" numerical filter is the no-flag default. Its
implementation record remains in
[`CANDIDATE_REJECTOR_DOUBLE.md`](CANDIDATE_REJECTOR_DOUBLE.md); explicit
`--unsafe` still selects the separate heuristic.

## Agreed Production Architecture

FracESSA no longer creates the complete set of nonempty supports or even a
complete cardinality layer before solving. A support is one `uint64_t` mask.
Each concrete generator owns the matrix dimension, the cardinality sweep, its
recursive state, and the forbidden-support rules discovered during analysis.

- `NonCircularSupportGenerator` uses fixed-cardinality binary DFS.
- `CircularSupportGenerator` uses fixed-content FKM-style necklace recursion
  followed by reflection reduction, so it emits one bracelet representative.
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
synchronously while the DFS or FKM call stack is still alive. After the callback
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

`fracessa::analyze_support()` performs the selected bounded-error or unsafe filter,
exact rational equilibrium solve, candidate population, and stability classification.
If it reports an exact candidate, the caller registers `candidate_.support` with
the active generator and then calls `finalize_candidate()`.

Every exact equilibrium support becomes forbidden for larger cardinalities,
whether or not stability classifies that candidate as an ESS. A result accepted
only by exact rational analysis may become a pruning rule; neither numerical
filter creates one by itself.

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

The same reasoning applies to the experimental third implementation: it exists
for a later comparison, not as justification for a hierarchy.

## Circular Output And Multiplier

Circular matrices now store only the smallest integer support in each complete
rotation/reflection orbit. The candidate's non-null `multiplier` is the number
of distinct raw supports represented by that row; it can be one. Non-circular
candidates use a null multiplier.

Candidate and ESS totals remain mathematical totals, so circular ESS rows are
counted with their multiplier. The nullable field replaced `shift_reference`
through the CLI, pybind boundary, Python APIs, Arrow and Parquet sinks, CSV
verification baseline, and SQLite fixtures.

Compressed output and pruning are separate concerns. The production circular
generator still registers every distinct rotation and reflection of an exact
candidate support as compact forbidden masks. One canonical mask alone is not
enough because canonicalization does not preserve subset relationships.

## Experimental `CircularSupportGeneratorV2`

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

The compact representation is not assumed to be faster: it saves stored orbit
masks but performs more bit-parallel work in the recursive pruning test. A
future focused comparison test and end-to-end benchmark must decide whether V2
or the expanded production representation should remain. Until then,
`CircularSupportGenerator` is the only production circular path.

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

## Remaining Work

- Add the planned permanent V1-versus-V2 comparison test when the generator test
  framework is ready, then benchmark and keep only the justified representation.
- Do not replace the callback with `next_support()` without a measured reason
  that justifies a resumable generator implementation.
- Do not funnel generator registration into `analyze_support()` unless the
  surrounding architecture changes enough to remove, rather than hide, the
  current type differences.
- After this branch is merged and updated `main` is verified, remove the
  `choice-one` worktree and delete the local and remote feature branch.
