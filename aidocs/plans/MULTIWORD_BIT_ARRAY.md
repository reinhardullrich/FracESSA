# Multiword Bit Array Plan

Status: implementation plan; no source changes have been made from this plan.

## Goal

Remove the accidental dimension-63 representation limit without slowing FracESSA's normal one-word support path.

- Dimensions 1 through 64 use exactly one `uint64_t` support word.
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

No code may shift a 64-bit word by 64. Dimension 64 also requires result-structure storage for index 64; the current
`std::array<size_t, 64>` objects stop at index 63.

### One-word mask microbenchmark

On 2026-08-08, GCC 16.1.1 compiled both formulas with `-O3 -DNDEBUG -march=native` and the benchmark ran pinned to CPU 2 on the
ARM64 laptop. Three independent runs each measured 10.8 billion materialized inline mask calculations, with nine alternating rounds
per formula and dimensions 1, 7, 31, and 63.

```cpp
// Current, valid for n < 64.
(uint64_t{1} << n) - 1

// Branchless, valid for 1 <= n <= 64.
~uint64_t{0} >> (64 - n)
```

The current expression measured about 0.458-0.461 ns per calculation. The branchless expression measured about 0.463-0.464 ns,
roughly 0.4-1.3% slower, or only about 3-5 ms per billion isolated calculations. Generated ARM64 code was:

```text
current:    mov + lsl + sub
branchless: mov + neg + lsr
```

The difference is tiny but consistently does not favor the replacement. Stage 2 should therefore calculate the dimension mask once
per matrix with the readable `n == 64` special case and reuse it in the per-support path. This both supports dimension 64 and removes
the repeated mask calculation instead of making it marginally slower.

## Ponytail Audit

### Delete before generalizing

The following `bitset64.hpp` helpers have no production caller. Their remaining callers are tests that test only the dead helper:

- `two_to_the_power_of()`;
- `find_pos_next_set_bit()`;
- `next_same_cardinality()`;
- `bits_before_pos()`;
- `is_smallest_representation()`.

Delete those helpers and their dedicated tests. Do not port them to the multiword type.

`CircularSupportGeneratorV2` in `supports.hpp` has no caller in production or tests. Its failed experiment is already preserved in
the historical documentation. Delete the implementation instead of porting it.

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
using bitset64 = uint64_t;       // n <= 64
class multiword_bitset;         // n >= 65
```

`multiword_bitset` owns `std::vector<uint64_t> words_`. Bit `i` is stored in:

```text
word index = i / 64
bit offset = i % 64
```

Its constructor receives the dimension, allocates `(dimension + 63) / 64` words once, and records a mask for the unused high bits of
the last word. Every mutating operation must preserve zero in those unused bits.

Required operations are only the operations with current production callers:

- set and clear one bit;
- fill the first `n` bits;
- empty/equality/numeric-order comparison;
- subset and set difference;
- population count;
- lowest set position and extraction of all set positions;
- lowest set bit removal for graph traversal;
- cyclic rotate and reflection;
- decimal and diagnostic bitstring output.

Numeric ordering compares words from highest to lowest. This preserves the current meaning of ascending support masks.

Do not expose arbitrary shifts, proxy references, iterators, hashing, arithmetic, or a general bitset algebra until a production
caller needs them.

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
3. copy the mask only when emitting a support or storing a newly forbidden exact-equilibrium support.

The one-word generators should retain their current by-value operations. Do not force the small path through the mutable multiword
implementation merely to make both source paths look identical.

## Migration Stages

Each stage is independently reviewable and requires approval for its exact source files.

### Stage 1: Remove dead one-word code

Files:

- `cpp/include/fracessa/bitset64.hpp`;
- `cpp/tests/test_bitset64.cpp`;
- `cpp/include/fracessa/supports.hpp` for unused V2 only.

Delete only the dead helpers, their tests, and `CircularSupportGeneratorV2`. Run the complete C++ tests. This stage should have no
runtime behavior change.

### Stage 2: Make dimension 64 a valid one-word boundary

Keep support masks as raw `uint64_t`. Special-case and cache the all-bits mask once per matrix, extend support-size result storage
through index 64, and change the validating parser and boundary tests from 63 to 64.

Audit every shift with dimension 64, especially rotation and reflection. Add quick dimension-64 generator tests that prune after the
singleton layer rather than attempting the complete search.

The canonical SQLite candidate table currently stores supports in signed `INTEGER`, so masks with bit 63 set do not fit. Do not
migrate the database merely to complete this runtime stage. Keep canonical database rows at dimension at most 63 until a separate
database-format change is approved.

### Stage 3: Add and prove the multiword primitive

Add the minimal multiword type and focused tests without connecting it to FracESSA yet. Required boundaries:

- dimensions 65, 128, 129, and a non-power-of-two dimension;
- bits 0, 63, 64, 127, and 128;
- unused high-bit masking;
- subset, difference, count, extraction, ordering, rotation, reflection, and round trips.

No benchmark is needed for an unconnected primitive. The tests must compare operations against a simple independent boolean-vector
reference.

### Stage 4: Extend the copositivity graph first

This is the smallest useful production slice. Replace the fixed 64-row sign-scan storage only on the large path:

- one-word matrices continue to use the existing fixed `uint64_t` neighbor rows;
- larger matrices use one multiword neighbor mask per row and dynamically sized row sums;
- connected-component traversal uses the same mathematical algorithm.

Test disconnected and connected negative-entry graphs across bits 63/64 and across several words. The cone algorithm itself has no
64-dimensional bitset limit; `pending.reserve(64)` is only an initial vector capacity.

### Stage 5: Extend complete analyzer support masks

Compile the active non-circular generator, V3 circular generator, circular affine symmetry, candidate search, candidate storage, and
stability orchestration for the multiword mask. Replace fixed-size support-index scratch arrays only on the large path with reusable
`std::vector<size_t>` storage.

Keep the small path's stack arrays and `uint8_t` indices unchanged. No per-support allocation is allowed in either generator or
candidate solver.

Result structures become dimension-sized vectors. CLI support fields remain decimal integers. Pybind converts a multiword decimal
value to Python's arbitrary-size `int`, preserving the public Python type.

The Parquet candidate schema currently uses `uint64`. Supporting larger masks requires a separately reviewed schema decision. The
minimal general representation is a decimal string; CSV and JSON already accept arbitrary-size integer text. Do not silently
truncate or split the value into undocumented columns.

### Stage 6: Remove the parser limit and update contracts

Remove the fixed maximum only after every runtime path above is connected. Keep overflow checks for `n * (n + 1) / 2` and normal
allocation failures at the validating input boundary.

Update `README.md`, `aidocs/PROJECT.md`, Python documentation, and parser errors. State explicitly that dimensions above 64 are
representable but exhaustive runtime remains exponential.

## Correctness Checks

The minimum acceptance set is:

1. Existing C++/CLI and Python suites remain green.
2. Existing candidate output is byte-identical for all current dimension-at-most-63 database matrices.
3. Dimension 64 uses the one-word path and handles bit 63 without undefined shifts.
4. Dimension 65 uses two words and handles supports crossing the 63/64 boundary.
5. Non-circular and circular generators preserve cardinality-first and numeric-within-cardinality order.
6. Forbidden-support pruning removes every strict superset on both paths.
7. Rotational/reflection multipliers and affine symmetry images remain exact.
8. Copositivity connected components match an independent queue-and-boolean-adjacency reference above 64.
9. ASan/UBSan tests cover dimensions 64, 65, 128, and 129.

Use full-support or immediately pruned matrices for analyzer boundary tests. Never try to verify dimension 65 by enumerating all
`2^65` supports.

## Performance Gate

Benchmark before and after Stages 2 and 5 with the canonical CPU-2 persistent-Pybind method.

- The dimension-at-most-63 path must show no repeatable regression on representative small, medium, circular, and non-circular
  matrices.
- Inspect the optimized one-word set/subset/extraction kernels if timing moves; they must not contain a large-mask branch or heap
  access.
- Multiword performance is secondary, but its recursion must allocate no memory per generated partial support.

If the shared template changes one-word code generation or produces a measurable regression, keep separate one-word and multiword
implementations. A little source duplication is cheaper than slowing the operation executed billions of times.

## Recommended First Approval

Approve only Stage 1 first. It removes code that otherwise would be unnecessarily ported. After its diff and tests are reviewed,
implement Stage 2 and establish dimension 64 before introducing dynamic storage.
