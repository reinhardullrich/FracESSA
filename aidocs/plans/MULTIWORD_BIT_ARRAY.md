# Multiword Bit Array Plan

Status: staged implementation for dimension 65 and above. Dimension 64, the isolated multiword primitive, and the dimension-independent
exact Hadeler checker were completed on 2026-08-08. Candidate search and public analyzer dispatch remain limited to dimension 64.

## Goal

Extend the completed one-word analyzer beyond dimension 64 without slowing its normal support path.

- Dimensions 1 through 64 use exactly one `uint64_t` support word; this boundary is outside the remaining work described here.
- Dimension 65 and above use `ceil(n / 64)` words.
- The choice is made once before the search. A normal one-word set operation must not inspect a variant, pointer, vector, or word count.
- Correctness comes first. The existing one-word performance is the second requirement.

This is a representation extension, not a claim that exhaustive search at dimension 65 is practical. A complete search still has up
to `2^n` supports. Larger dimensions are nevertheless useful for full-support checks, heavily pruned searches, and the exact
copositivity checker.

## Boundary: 64, Not 63

A `uint64_t` represents strategy positions 0 through 63, so it represents every support of a dimension-64 game. The multiword path
must begin at dimension 65.

The old dimension-63 limit came from using `1ULL << n` as an exclusive end value for complete mask enumeration. Production no longer
enumerates supports that way:

- non-circular supports use binary depth-first generation;
- circular supports use direct fixed-density bracelet generation;
- `bs64::two_to_the_power_of()` has no production caller.

Dimension 64 therefore needs only correct boundary handling:

```cpp
set_all_n_bits(64) == ~uint64_t{0}
```

No code may shift a 64-bit word by 64. The implemented result structures have 65 entries so support-size index 64 is valid.

### One-word mask microbenchmark

On 2026-08-08, GCC 16.1.1 compiled all three formulas with `-O3 -DNDEBUG -march=native` and the benchmark ran pinned to CPU 2 on the
ARM64 laptop. Each timing window executed 3 billion materialized inline mask calculations and lasted about 1.35-1.39 seconds. The
balanced comparison executed 120 billion calls in total.

```cpp
// Former implementation, valid for n < 64.
(uint64_t{1} << n) - 1

// Implemented explicit boundary, valid for n <= 64.
n == 64 ? ~uint64_t{0} : (uint64_t{1} << n) - 1

// Branchless, valid for 1 <= n <= 64.
~uint64_t{0} >> (64 - n)
```

For dimensions 1, 7, 31, and 63, the former expression measured about 0.459-0.460 ns, the explicit condition about 0.450 ns, and the
branchless expression about 0.464 ns per calculation. In this generated loop, the explicit condition was about 2.0-2.2% faster than
the former expression; the branchless form was about 0.8-0.9% slower. At dimension 64, the former expression is invalid, while the
explicit and branchless forms measured 0.421 ns and 0.464 ns respectively. The one-word implementation therefore uses the readable
explicit condition directly; no cached-mask member was added.

## Ponytail Audit

### Delete before generalizing

The following `bitset64.hpp` helpers have no production caller. Their remaining callers are tests that test only the dead helper:

- `two_to_the_power_of()`;
- `find_pos_next_set_bit()`;
- `bits_before_pos()`;
- `is_smallest_representation()`.

Delete those helpers and their dedicated tests. Do not port them to the multiword type.

`CircularSupportGeneratorV2` had no caller in production or tests and was removed. Its failed experiment remains preserved in the
historical documentation. Do not restore or port it.

`CircularSupportGenerator` V1 is test-only, but it remains a useful independent oracle for V3. Keep it one-word-only during this
work; do not generalize it. Production uses only `NonCircularSupportGenerator` and `CircularSupportGeneratorV3`.

### Keep

The remaining helpers are active in at least one of support generation, circular symmetry, candidate solving, stability, output, or
logging. In particular:

- `is_set_at_pos()` is indirectly used by runtime bitstring logging;
- `popcount64()` is the one-word primitive that a multiword count also needs;
- `ctz64()` is used by support pruning, candidate indexing, and connected components;
- rotations and reflection are active in circular support generation and symmetry handling.

No standard C++17 type supplies a runtime-sized bit array. `std::bitset` has a compile-time size. Adding Boost only for
`dynamic_bitset` would add a dependency and would not by itself preserve the raw-`uint64_t` hot path. A small word-vector type is the
minimal required implementation.

## Representation

Keep two concrete types:

```cpp
using bitset64 = uint64_t;       // production path for n <= 64
class bitset_multiword;         // production path for n >= 65
```

`bitset_multiword` owns `std::vector<uint64_t> words_`. Bit `i` is stored in:

```text
word index = i / 64
bit offset = i % 64
```

Its constructor receives the dimension, allocates `dimension / 64 + (dimension % 64 != 0)` words once, and records a mask for the
unused high bits of the last word. This overflow-safe expression matters before the parser has a separately approved practical
maximum. Production constructs this type only for dimensions above 64, but the primitive accepts every positive dimension so its
later generators can be compared exhaustively with the one-word implementation at small dimensions. Every mutating operation must
preserve zero in unused high bits.

### Fixed width for one complete matrix run

The dimension is known before support generation starts. Calculate

```cpp
word_count = (dimension + 63) / 64;
```

once and keep it unchanged for the complete matrix run. Every support belonging to that matrix uses exactly `word_count` words. A
singleton support in a dimension-129 game therefore still has three words, just like its full support. Cardinality changes the
number of set bits, not the representation width.

The generator owns one pre-sized mutable working mask. Long-lived masks—such as forbidden supports and retained candidates—each own
their own `word_count` words. They are copied only when they must outlive the synchronous generator callback. Never resize a support
because its current highest set bit is small; that would complicate comparisons, rotations, complements, and unused-bit invariants
for no useful gain.

Required operations are only the operations with current production callers:

- set and clear one bit;
- fill the first `n` bits;
- empty/equality/numeric-order comparison;
- subset, union, and set difference;
- population count;
- lowest set position and extraction of all set positions;
- lowest set bit removal for graph traversal;
- in-place cyclic rotation by one position and reflection;
- diagnostic bitstring output.

Numeric ordering compares words from highest to lowest. This preserves the current meaning of ascending support masks.

Extraction writes into a caller-owned, already reserved `std::vector<size_t>` so repeated candidate solves do not allocate. Decimal
conversion is deliberately absent from this primitive. Stage 6 imports its fixed words into the existing FLINT integer wrapper only
at the CLI/Pybind output boundary.

Do not expose arbitrary shifts, proxy references, iterators, hashing, arithmetic, or a general bitset algebra until a production
caller needs them.

## Concrete Template Boundary

The dimension cannot be a template parameter because it comes from the input at runtime. Do not use `std::bitset<N>` and do not
instantiate one analyzer for every possible dimension. Templates select only the storage representation:

```cpp
template<class SupportMask>
class basic_fracessa;

template<class SupportMask>
struct basic_candidate;

using fracessa = basic_fracessa<bitset64>;
using multiword_fracessa = basic_fracessa<bitset_multiword>;
```

The existing `fracessa` alias preserves the current one-word C++ tests and generated code. CLI and Pybind dispatch once after parsing
the dimension:

```cpp
if (dimension <= 64)
    run_typed<bitset64>(...);
else
    run_typed<bitset_multiword>(...);
```

`find_candidate_fast`, `find_candidate_test`, and `find_candidate_safe` should keep owning their present numerical matrices and
workspaces. Only their support-dependent member functions become function templates over `SupportMask`. `check_stability()` and the
main analyzer orchestration use the same `SupportMask` as `basic_candidate`. Define these templates in implementation files and add
explicit instantiations for the two mask types where practical; do not move unrelated numerical kernels into public headers.

Share small read operations through overloaded inline functions or a two-specialization `support_mask_ops<SupportMask>`:

- set, clear, and test one bit;
- count and extract set positions;
- subset, difference, equality, and numeric order;
- rotate, reflect, and convert to output text.

Do not implement another arbitrary-precision decimal converter inside the bitset. At the CLI/Pybind output boundary, reuse the
existing FLINT integer wrapper to import the fixed words and produce decimal text; Python then constructs its arbitrary-size `int`
from that exact text. The hot support operations remain independent of FLINT conversion.

Do not force the generators behind one identical implementation. Their storage economics are different:

- the current `uint64_t` generators keep their by-value recursion unchanged;
- the multiword generators mutate one fixed-width working mask and undo each bit before returning;
- circular affine symmetry gets a separate multiword specialization because its current table of one-word destination masks is not
  an efficient representation above 64.

## Static Dispatch, Not A Runtime Variant

Do not replace every `bitset64` with one class containing
`std::variant<uint64_t, std::vector<uint64_t>>`. That design makes every small operation branch and makes every small support object
large enough to contain a vector.

Dispatch once from the analyzer boundary:

```cpp
if (dimension <= 64)
    run_one_word_search(...);
else
    run_multiword_search(...);
```

Shared search logic may be compiled for both concrete mask types through small templates or overloads. The generated one-word code
must still operate directly on `uint64_t`. Use templates only where both paths perform the same algorithm; do not build a hierarchy,
virtual interface, allocator framework, or generic bitset library.

## Multiword Generator Rule

The current one-word generators pass masks by value during recursion. Copying a `std::vector<uint64_t>` at each recursive step would
be catastrophically expensive.

The multiword generators must instead own one preallocated working mask:

1. set a bit before descending into the include branch;
2. clear that bit after returning;
3. pass the completed mask to the callback by `const&`;
4. copy the mask only when storing a newly forbidden exact-equilibrium support or a retained candidate.

The callback is synchronous. Its reference becomes invalid as soon as the callback returns and generation continues. Candidate
solvers may read it without copying. A successful exact candidate copies it once into the working candidate; `add_forbidden()` then
copies it once into persistent pruning storage. Candidate output copies it only when candidate output or logging was requested.

The one-word generators should retain their current by-value operations. Do not force the small path through the mutable multiword
implementation merely to make both source paths look identical.

For multiword circular generation, allocate the V3 `positions` and `prefix_density` work arrays once at `dimension + 2`. Both arrays
must use `size_t`; the current `uint8_t` density entries cannot represent dimensions above 255. Circular affine symmetry should
store destination **positions**, not one separately allocated single-bit multiword mask per source position. Transform a support by
iterating its set positions and setting their mapped destination bits in a pre-sized scratch mask.

## File-by-File Implementation Map

| File | Dimension-above-64 change |
|---|---|
| `cpp/include/fracessa/bitset_multiword.hpp` | Add the fixed-width-per-run word vector and only the required operations. |
| `cpp/tests/test_bitset_multiword.cpp` | Compare every operation with an independent `std::vector<bool>` reference at the 63/64, 127/128, and unused-high-bit boundaries. |
| `cpp/tests/CMakeLists.txt` | Build and register the isolated primitive test. |
| `cpp/include/fracessa/supports.hpp` | Leave the one-word generators intact; add separate mutable-work-mask non-circular and V3 circular generators. Do not restore or port V2. |
| `cpp/include/fracessa/circular_affine_symmetry.hpp` | Keep the present one-word implementation and add a multiword specialization with destination-position tables and reusable scratch masks. |
| `cpp/include/fracessa/candidate.hpp` | Make the working candidate depend on `SupportMask`; keep support and extended support in the same fixed-width representation. |
| `cpp/include/fracessa/fracessa.hpp`, `cpp/src/fracessa.cpp`, `cpp/src/checkstab.cpp` | Template the search engine over the mask type, keep one static dispatch before the search, and replace result-structure arrays with dimension-sized vectors. |
| `cpp/include/fracessa/find_candidate_*.hpp`, `cpp/src/find_candidate_*.cpp` | Template only support-dependent methods. Give the large instantiation reusable `size_t` index and numerical scratch vectors; retain all one-word stack arrays. |
| `cpp/include/linalg/copositive_integer.hpp` | Add recursive large-dimension Hadeler enumeration and multiword negative-neighbor components; retain the one-word Hadeler implementation unchanged. |
| `cpp/src/main.cpp`, `cpp/src/pybind_module.cpp` | Dispatch once, serialize arbitrary-width masks, and use dimension-sized result structures. |
| `cpp/include/fracessa/matrix_parser.hpp` | Remove the bitmask-derived maximum only after all paths are connected; add checked triangular-size arithmetic. |
| `python/pyfracessa/sinks_parquet.py` | Reject an out-of-range support explicitly until a separately approved schema replaces `uint64`; CSV and JSON need no schema change. |

The large candidate solvers must not allocate index or numerical scratch storage for every support. Resize their vectors once from the
game dimension and reuse them. In particular, replace the large-path forms of `support_indices`, `non_support_indices`, `solution`,
`scale_ratios`, `pivots`, `game_scales`, and reduced-B component indices. The safe solver's dense `dimension^3` reduced-entry cache
also needs checked size multiplication; its memory growth is a separate practical limit even after support masks become multiword.

## Migration Stages

Each stage is independently reviewable and requires approval for its exact source files.

### Completed one-word prerequisite

Dimension 64 works with raw `uint64_t`, including bit 63, rotation, reflection, extraction, full-support construction, result
structures through support size 64, parser and CLI/Pybind boundaries, and the copositivity component mask. Focused generator tests
prune after the singleton layer and never attempt to enumerate `2^64` supports. The remaining stages do not reimplement that case.

Remove the dead `bitset64.hpp` helpers listed in the Ponytail audit in a separately approved cleanup before templating or duplicating
their callers. V2 has already been removed.

The canonical SQLite candidate table may remain restricted to dimension at most 63 because signed SQLite `INTEGER` cannot store a
mask with bit 63 set. Runtime support and canonical database storage are separate changes.

### Completed Stage 1: Add and prove the multiword primitive

The minimal multiword type and focused tests were added without connecting it to FracESSA. Covered boundaries are:

- dimensions 65, 128, 129, and a non-power-of-two dimension;
- bits 0, 63, 64, 127, and 128;
- unused high-bit masking;
- subset, union, difference, count, extraction, ordering, one-position rotation, reflection, and bitstring round trips.

The tests compare every operation against a simple independent boolean-vector reference. Reflection is checked at every dimension
from 1 through 260, including unused-word boundaries, and is verified as an involution. Release and ASan/UBSan tests pass.

The initial per-bit reflection was replaced with in-place word-order reversal, one bit reversal per word, and one fused alignment
pass. On the ARM64 laptop, representative 65-256-bit kernel measurements were about 7-10x faster than the per-bit implementation
except at dimension 70, which was about 2.5x faster. The shared one-word reversal uses ARM64's single `rbit` instruction and retains
the portable mask-and-shift implementation elsewhere. A matched CPU-2 persistent-Pybind benchmark over ten representative circular
and non-circular matrices reported a median current/baseline ratio of exactly 1.000 in both fast and safe modes. This confirms that
the isolated stage did not regress the existing production path; no multiword analyzer benchmark exists because that path is not
connected yet.

### Completed Stage 2: Extend exact Hadeler copositivity

Hadeler enumeration and the fixed 64-row sign-scan storage now have a separate large-dimension path:

- one-word matrices continue to use the existing fixed `uint64_t` neighbor rows;
- larger matrices recursively enumerate fixed-cardinality principal subsets and carry their selected positions directly in one
  pre-sized `std::vector<size_t>`; no subset mask, index extraction, or per-subset allocation is needed;
- larger matrices use one multiword neighbor mask per row and dynamically sized row sums;
- connected-component traversal uses the same mathematical algorithm and reuses its component, frontier, and discovery masks;
- only graph-discovered components require one index extraction before constructing their dense principal matrices.

Focused tests cover a bad pair across positions 63/64, a cardinality-three rejection spanning positions 0/64/128 after all proper
principal subsets pass, a disconnected component spanning three words, and a connected 65-vertex negative graph. Every case decides
before exponential enumeration becomes impractical. Release, ASan/UBSan, Python, and database-integrity checks pass.

The checker remains in FracESSA because no separate Coposit repository currently owns this implementation.

### Stage 3: Template the candidate and stability core

Introduce `basic_candidate<SupportMask>` and `basic_fracessa<SupportMask>`. Template only the support-dependent methods of the fast,
test, and safe candidate finders. Add the reusable large-path index and numerical scratch vectors and the dimension-sized candidate
and ESS structures.

Connect this stage first through internal full-support tests at dimensions 65, 128, and 129. Full support needs no combinatorial
generator, so it isolates mask extraction, candidate solving, extended-support construction, reduced-B stability, and output-mask
ownership. Do not remove the public parser limit in this stage.

### Stage 4: Add the multiword non-circular generator

Implement the separate mutable-work-mask depth-first generator. Preserve cardinality-first and numeric-within-cardinality order,
activate forbidden supports only between cardinalities, and copy only persistent forbidden supports.

Tests must stop cheaply: compare complete generation only at small dimensions by translating one-word masks, and at dimensions 65
and 129 prune after the singleton or another deliberately tiny layer. Verify forbidden supports spanning bits 63 and 64.

### Stage 5: Add multiword circular generation and symmetry

Implement the multiword V3 generator with dimension-sized `size_t` work arrays, then the multiword affine-symmetry specialization.
Preserve bracelet ordering, dihedral multipliers, reconstructed candidates, and pruning across every stored rotation and reflection.
Do not port V1 or V2 to multiword storage.

Tests should compare the multiword path against the one-word V3 oracle at feasible dimensions, then use immediately pruned
dimension-65 and dimension-129 cases to cover the word boundary without enumerating all bracelets.

### Stage 6: Expose the large path and update output contracts

Add the one-time CLI and Pybind dispatch and remove the parser's bitmask-derived maximum. Keep checked arithmetic for
`n * (n + 1) / 2`, word counts, dense cache dimensions, result counters, and allocation sizes. Allocation failure remains an explicit
error, not a fallback to truncated storage.

Keep the one-word path's stack arrays and `uint8_t` indices unchanged. No per-support allocation is allowed in either generator or
candidate solver. The large path uses fixed-width-per-run masks and reusable `size_t` buffers.

Result structures become dimension-sized vectors. CLI support fields remain decimal integers. Pybind converts a multiword decimal
value to Python's arbitrary-size `int`, preserving the public Python type.

The Parquet candidate schema currently uses `uint64`. Supporting larger masks requires a separately reviewed schema decision. The
minimal general representation is a decimal string; CSV and JSON already accept arbitrary-size integer text. Until that decision,
the Parquet sink must reject an out-of-range support clearly. Do not silently truncate or split the value into undocumented columns.

Update `README.md`, `aidocs/PROJECT.md`, Python documentation, and parser errors. State explicitly that dimensions above 64 are
representable but exhaustive runtime remains exponential.

## Correctness Checks

The minimum acceptance set is:

1. Existing C++/CLI and Python suites remain green.
2. Existing candidate output is byte-identical for all current one-word database matrices.
3. Dimension 64 remains on the completed one-word path; the multiword path is never instantiated for it.
4. Dimensions 65 and 128 use two fixed words for every support in their respective runs; dimension 129 uses three.
5. Supports crossing the 63/64 and 127/128 boundaries survive set, extraction, complement, rotation, reflection, decimal-output,
   and Python-integer round trips.
6. Non-circular and circular generators preserve cardinality-first and numeric-within-cardinality order.
7. Forbidden-support pruning removes every strict superset on both paths.
8. Rotational/reflection multipliers and affine symmetry images remain exact.
9. Copositivity connected components match an independent queue-and-boolean-adjacency reference above 64.
10. ASan/UBSan tests cover dimensions 65, 128, and 129, while the existing one-word suite retains dimension 64.
11. Instrumented generator tests observe no allocation during recursive descent after construction; only persistent forbidden or
    candidate copies may allocate.
12. Parquet either uses an approved arbitrary-width schema or rejects a support above `uint64` explicitly.

Use full-support or immediately pruned matrices for analyzer boundary tests. Never try to verify dimension 65 by enumerating all
`2^65` supports.

## Performance Gate

Benchmark the one-word baseline before Stage 3 and again after Stages 3 and 6 with the canonical CPU-2 persistent-Pybind method.

- The dimension-at-most-63 path must show no repeatable regression on representative small, medium, circular, and non-circular
  matrices.
- Inspect the optimized one-word set/subset/extraction kernels if timing moves; they must not contain a large-mask branch or heap
  access.
- Multiword performance is secondary, but its recursion must allocate no memory per generated partial support.

If the shared template changes one-word code generation or produces a measurable regression, keep separate one-word and multiword
implementations. A little source duplication is cheaper than slowing the operation executed billions of times.

## Next Approval

Remove the remaining dead one-word helpers in one isolated cleanup before Stage 3. Then implement the candidate and stability
templates separately; do not combine that work with either support generator.
