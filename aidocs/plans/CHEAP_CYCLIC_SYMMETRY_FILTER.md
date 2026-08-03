# Cheap Cyclic Symmetry Filter

Status: implemented and always active for circular-symmetric matrices. Correctness verification and the experimental on/off
benchmark are complete.

## Goal

Circular-symmetric games currently solve one support per binary bracelet, so rotations and reflections are already removed by
`CircularSupportGeneratorV3`. Some circular matrices have further exact symmetries because their repeated payoff values are
arranged compatibly with multiplication of cyclic strategy indices.

The first implementation should detect those cheap symmetries once per game and discard equivalent bracelet representatives
after V3 generates them but before `fracessa::analyze_support()` solves them.

The intended pipeline is:

```text
V3 generates one dihedral bracelet representative
    -> exact affine-symmetry filter
        -> skip an equivalent representative, or
        -> solve one retained representative
            -> expand an exact candidate over its full detected orbit
            -> register every orbit member for later superset pruning
            -> output one reconstructed row per dihedral orbit
```

This is an additional orbit reduction around V3, not a replacement for V3.

## Always-On Circular Behavior

The cyclic symmetry filter runs automatically for every matrix already identified as circular symmetric. There is no CLI, Python,
or C++ switch for it.

For a non-circular matrix the helper is never constructed; it must not attempt to infer circular structure from a general matrix.
For a circular matrix with no additional multiplier class, it disengages immediately after the one-time exact detection pass and
the legacy V3 path continues unchanged. The filter does not change the selected `fast`, `safe`, or `test` candidate method.

The former experimental switch existed only to complete the correctness and performance comparison. It was removed after 101 of
120 stored circular matrices were found to have no additional multiplier class and the completed no-extra control measurements
showed median on/off ratios of 1.0000 in both fast and safe mode.

## Scope

The first version will detect affine multiplier symmetries of the cyclic indices:

$$
i\longmapsto ai+b\pmod n,
$$

where $\gcd(a,n)=1$. Translations $b$ are already rotations, and $a=-1$ is already reflection. The useful new cases are valid
multipliers other than $\pm1$.

The first version will not:

- run a general colored-graph automorphism package;
- recognize rook graphs or product-grid groups such as $S_3\times S_8$;
- replace the fixed-density bracelet recursion with a generalized canonical generator;
- assume that repeated payoff values alone create a symmetry;
- precompute or store exponential support tables.

Those exclusions are intentional. The affine filter is exact, small, dependency-free, and directly testable. More general matrix
automorphisms can be considered only after this version has a measured benefit.

## Mathematical Basis

For a circular matrix, index strategies by $\mathbb Z_n$ and write

$$
A_{ij}=c_{j-i\bmod n}.
$$

A unit $a\in(\mathbb Z/n\mathbb Z)^\times$ is a matrix automorphism exactly when

$$
c_{ad\bmod n}=c_d
\qquad\text{for every }d\in\mathbb Z_n.
$$

Equivalently, it is sufficient to compare the exact entries in row zero:

$$
A_{0,d}=A_{0,ad\bmod n}
\qquad\text{for every }d.
$$

Define the detected multiplier group

$$
H=\left\{a\in(\mathbb Z/n\mathbb Z)^\times:A_{0,d}=A_{0,ad}\text{ for every }d\right\}.
$$

Together with translations, these multipliers give the verified affine group

$$
G=\mathbb Z_n\rtimes H.
$$

The matrix comparisons must use `fraction` equality. There is no tolerance and no double conversion. At the supported dimensions,
testing every possible multiplier costs only $O(n\varphi(n))$ exact comparisons once per game.

Because the matrix is symmetric, $a$ and $-a$ produce supports in the same dihedral orbit. Store only one representative of each
pair in $H/\{\pm1\}$. If this quotient contains only the identity class, the helper is inactive and every generated bracelet passes
immediately.

### Why filtering bracelets is correct

Let $D_n$ be the rotation/reflection group already handled by V3, and let

$$
C_D(S)=\min\{gS:g\in D_n\}
$$

be the smallest integer mask in the dihedral orbit of support $S$. V3 emits exactly these masks.

The affine multipliers act on dihedral orbits. For each V3 result $S$, calculate

$$
C_G(S)=\min_{a\in H} C_D(aS).
$$

Solve $S$ exactly when $S=C_G(S)$. Every $G$-orbit therefore contributes one retained bracelet and all other equivalent
bracelets are skipped before any candidate system is built.

This is safe because every tested multiplier is a literal permutation $P$ satisfying

$$
P^TAP=A.
$$

It preserves the candidate equations, positivity and outside-payoff conditions, and ESS stability. A candidate or ESS on $S$
exists exactly when the permuted candidate or ESS exists on $aS$.

## Current Insertion Point

The circular branch in `cpp/src/fracessa.cpp` currently has the exact lifecycle needed by this change:

1. `CircularSupportGeneratorV3::generate()` emits one bracelet and its cardinality.
2. The callback calls `analyze_support()`.
3. An exact candidate is passed to `generator.add_forbidden()`.
4. `add_forbidden()` expands its complete dihedral orbit and returns the number of distinct concrete masks.
5. `finalize_candidate()` stores that number as the output multiplier and uses it for the ESS total.

The new filter belongs between steps 1 and 2. V3's recursion, reversal checks, fixed-density order, and forbidden-mask buckets should
remain unchanged.

## Proposed Design

### 1. Add one small circular-affine helper

Add a focused header under `cpp/include/fracessa/`, provisionally named `circular_affine_symmetry.hpp`. It should own:

- the dimension;
- the verified representatives of $H/\{\pm1\}$;
- the precomputed bit permutation for each retained multiplier;
- dihedral canonicalization of a transformed support;
- the two operations needed by the orchestration code:
  - `is_representative(support)`;
  - iteration over the distinct bracelet images of a support.

Keep this helper independent of candidate solving and stability. It deals only with exact matrix symmetry and support masks.

Precompute destination bits once, so applying $a$ to a support does not perform a modulo operation for every emitted bracelet:

```text
destination[a][i] = 1 << ((a * i) mod n)
```

Transform a support by walking only its set bits with `ctz64()` and OR-ing their precomputed destinations. The dimension is below
64, so fixed stack arrays are sufficient. No hash table or heap allocation is required in the per-bracelet path.

### 2. Detect valid multipliers once

Construct the helper once in the circular branch. For every $a$ from 1 to $n-1$:

1. reject it unless `std::gcd(a, n) == 1`;
2. compare `game_matrix_(0, d)` with `game_matrix_(0, (a * d) % n)` for all offsets $d$;
3. retain it only when every exact comparison succeeds;
4. collapse the equivalent pair $a$ and $n-a$.

The existing `is_cs` contract already says that the input is symmetric circulant. This feature should not add a second full
circulant-matrix validator to the hot program path.

### 3. Filter immediately before solving

At the beginning of the V3 callback:

```cpp
if (!symmetry.is_representative(support))
    return;
```

`is_representative()` should return `true` immediately when no extra multiplier exists. Otherwise, for each retained multiplier it
transforms the support and scans its rotations and reflected rotations. It can reject as soon as any concrete image is numerically
smaller than the already-canonical V3 support. It does not need to allocate or materialize the complete orbit.

V3's `emitted_` state must continue to describe bracelets generated by V3, not bracelets accepted by this outer filter. Therefore,
the filter must stay outside the recursion and must not affect V3's cardinality early-stop rule.

### 4. Reuse V3 for enlarged pruning while preserving dihedral output

When a retained support is an exact candidate, enumerate the distinct values

$$
C_D(aS),\qquad a\in H/\{\pm1\}.
$$

There are at most $\varphi(n)/2$ such bracelet images. Deduplicate them in a fixed local array; a short linear scan is simpler than
a hash table and candidate supports are much rarer than generated supports.

For every distinct bracelet image, call the existing

```cpp
generator.add_forbidden(image)
```

separately. This gives both required results without changing V3:

- every concrete support in the verified affine orbit becomes a future forbidden subset;
- every output row keeps its own dihedral multiplier returned by `add_forbidden()`, which is at most $2n$.

Do not combine multiple dihedral orbits into one output multiplier. A count alone would not identify the matrix-specific affine
transformations needed to reconstruct those candidates. Instead, exactly permute the solved vector, support, and extended support
into each bracelet image and finalize one row per image. This adds no candidate-system solve.

Registering the complete affine orbit is also the invariant that makes the outer representative filter compatible with dynamic
superset pruning. The forbidden family remains closed under $G$: if pruning removes the smallest bracelet in an affine orbit, the
corresponding transformed forbidden support removes every other bracelet in that orbit as well. It is not correct to filter by the
larger group while continuing to register only the old dihedral candidate orbit.

The callback becomes conceptually:

```cpp
if (!symmetry.is_representative(support))
    return;

if (analyze_support(support, support_size)) {
    const candidate solved = candidate_;
    symmetry.for_each_distinct_bracelet_image(candidate_.support, [&](auto image) {
        candidate_ = permute(solved, image);
        finalize_candidate(generator.add_forbidden(image.support));
    });
}
```

The full-support shortcut remains unchanged: the full support has orbit size one under every permutation.

## Expected Example

For the published dimension-24 rook matrix, all units modulo 24 preserve the compact payoff pattern:

$$
H=\{1,5,7,11,13,17,19,23\}.
$$

After identifying $a$ with $-a$, four multiplier classes remain. The filter can therefore merge up to four ordinary bracelet
classes before solving.

It will not recover the complete $S_3\times S_8$ automorphism group. That larger group can merge up to 5,040 dihedral classes in a
generic orbit, but recognizing and generating under it is a separate project. The affine experiment must not claim those savings.

## Implemented Files

- `cpp/include/fracessa/circular_affine_symmetry.hpp`
  - exact multiplier detection;
  - precomputed bit permutations;
  - representative test;
  - distinct bracelet-image iteration.
- `cpp/src/fracessa.cpp`
  - construct the helper in the circular branch and disengage it when only the identity class exists;
  - filter before `analyze_support()`;
  - call `add_forbidden()` for every distinct bracelet image and sum the multipliers.
- `cpp/tests/test_supports.cpp`
  - focused group, canonicalization, orbit, and pruning tests.
- `cpp/tests/cli_parser_blackbox.py`
  - end-to-end automatic filtering and weighted ESS-total coverage.

No matrix-parser, candidate-solver, stability, or database-schema change was needed.

## Correctness Tests

### Detection

- A circular matrix with distinct distance values detects only $\{\pm1\}$ and leaves the helper inactive.
- An $n=8$ matrix with the required equality between distance classes 1 and 3 detects multiplication by 3.
- A matrix with repeated values in positions not preserved by a multiplier does not accept that multiplier.
- The published $n=24$ rook matrix detects all eight units and four classes modulo sign.

### Support orbits

For manageable dimensions, enumerate every nonempty support by brute force and independently form its complete affine orbit.
Verify that:

- V3 plus the filter retains exactly one representative per affine orbit;
- the retained representative is the smallest integer mask in that orbit;
- distinct-bracelet iteration covers the complete orbit exactly once after V3 expands each image;
- every reconstructed row has the same vector, support, and extended support as the legacy dihedral result;
- every individual multiplier is the size of one dihedral orbit and is at most $2n$.

### Forbidden-superset pruning

Add one candidate support through the new path and generate the next cardinalities. Confirm that no emitted support contains any
affine image of the candidate. Include a support whose rotations, reflections, and multiplier images have different smallest set
bits so every V3 forbidden bucket is exercised.

### End-to-end results

- During the experiment, run every affected fixture with the option off and confirm byte-identical legacy output.
- Compare affected circular matrices with the former filter-off output and require exact candidate rows after deterministic
  ordering.
- Run the published dimension-24 matrix and preserve its exact represented totals: 15,120 candidates and
  15,120 ESS.
- Verify candidate vectors, supports, payoffs, ESS flags, and weighted multipliers independently of candidate IDs.
- Run all existing CTests plus the maintained Python verification suite in Release and ASan/UBSan builds.

Candidate rows, IDs, vectors, supports, and individual multipliers must remain identical to the former dihedral-only output. The
filter changes only which systems are solved and which equivalent pruning rules are registered.

## Performance Experiment

Instrument the experimental build, not production, with:

- V3 bracelets emitted;
- bracelets rejected by the affine filter;
- calls to `analyze_support()`;
- exact candidates found;
- distinct affine bracelet images registered;
- concrete forbidden masks registered;
- time spent detecting multipliers;
- time spent in the representative filter;
- total end-to-end time.

Benchmark:

1. circular quick-test matrices with no extra multiplier, to measure the inactive-path overhead;
2. synthetic $n=8$ and other small matrices with known multiplier symmetry;
3. the published dimension-24 rook matrix (matrix 34);
4. any stored circular rook-family matrices for which runtime is practical;
5. both `fast` and `safe`, because the value of avoiding a solve differs greatly between them.

The experimental binary exposed an option-off baseline only for this comparison. Both states used identical iteration counts in
one persistent CPU-pinned process, and the order was reversed between passes. Production no longer exposes either state.

Use the project's persistent-process benchmark method, pin both compared builds to the same CPU, and compare per-matrix medians.
The experiment should also report solver-call reduction separately from runtime so a weak result can be diagnosed.

## Promotion Decision

The filter was promoted after all of the following held:

- exact expanded candidates and ESS are unchanged on every validation matrix;
- every multiplier and forbidden orbit passes independent brute-force checks;
- affected matrices show a material end-to-end gain;
- ordinary circular matrices have no meaningful regression;
- the implementation remains the small outer filter described here and does not complicate V3's recursion.

Promotion required:

1. update the circular multiplier documentation;
2. verify that existing candidate baselines and database rows remain unchanged;
3. refresh current calibration timings for affected matrices;
4. record the measured reduction and the fact that the method detects only the affine subgroup, not the complete automorphism
   group.

## Deferred Extensions

Only measured evidence should justify a second phase. Possible later work is:

- a special constant-off-diagonal case, where $\operatorname{Aut}(A)=S_n$ and one support per cardinality is sufficient;
- an exact recognizer for rectangular rook matrices and row/column support canonicalization;
- a general colored-graph automorphism tool;
- direct generation under a larger group instead of post-bracelet filtering.

These are deliberately not part of the first implementation.
