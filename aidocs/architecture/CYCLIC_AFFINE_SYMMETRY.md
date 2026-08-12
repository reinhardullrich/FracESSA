# Cyclic Affine Symmetry Reduction

Status: implemented and always active for circular-symmetric matrices.

## Purpose

`CircularSupportGenerator` already emits one support per rotation/reflection orbit, or bracelet. Some circular matrices have
additional exact symmetries because multiplication of cyclic strategy indices preserves their repeated payoff pattern.

FracESSA detects those symmetries once per circular game and uses them to:

- skip multiplier-equivalent bracelets before solving a candidate system;
- register every verified image of an exact candidate for later superset pruning; and
- reconstruct the same per-bracelet candidate rows and multipliers that the public output contract uses.

This is an exact outer reduction around the bracelet generator. It does not change candidate solving, stability, or the selected
`fast` or `safe` method.

## Mathematical Basis

Index the strategies of a circular matrix by $\mathbb Z_n$ and write

$$
A_{ij}=c_{j-i\bmod n}.
$$

A unit $a\in(\mathbb Z/n\mathbb Z)^\times$ is a matrix automorphism exactly when

$$
c_{ad\bmod n}=c_d
\qquad\text{for every }d\in\mathbb Z_n.
$$

It is therefore sufficient to compare exact entries in row zero:

$$
A_{0,d}=A_{0,ad\bmod n}
\qquad\text{for every }d.
$$

The detected multiplier group is

$$
H=\left\{a\in(\mathbb Z/n\mathbb Z)^\times:A_{0,d}=A_{0,ad}\text{ for every }d\right\}.
$$

Translations are already rotations. Because the matrix is symmetric, $a$ and $-a$ produce supports in the same dihedral orbit.
The implementation therefore retains one representative of each class in $H/\{\pm1\}$. Together with translations, the verified
affine group is

$$
G=\mathbb Z_n\rtimes H.
$$

All comparisons use exact `fraction` equality. Repeated payoff values alone are insufficient: their positions must satisfy the
multiplier identity above.

### Why One Bracelet Per Affine Orbit Is Enough

Let $D_n$ be the rotation/reflection group already handled by the circular generator, and define the canonical bracelet

$$
C_D(S)=\min\{gS:g\in D_n\}.
$$

For each generated bracelet $S$, define

$$
C_G(S)=\min_{a\in H} C_D(aS).
$$

FracESSA solves $S$ only when $S=C_G(S)$. Every verified affine orbit therefore contributes one solved bracelet.

This is correct because each accepted multiplier is a literal permutation matrix $P$ satisfying

$$
P^TAP=A.
$$

It preserves the candidate equations, positivity, outside-payoff inequalities, and ESS stability. A candidate or ESS on $S$
exists exactly when its permuted counterpart exists on $aS$.

## Implemented Flow

```text
CircularSupportGenerator emits one bracelet
    -> reject it if an exact affine image has a smaller bracelet
    -> analyze the retained support once
    -> for an exact candidate:
         canonicalize every distinct affine image to a bracelet
         register each bracelet's complete dihedral orbit for superset pruning
         reconstruct and output one candidate row per distinct bracelet
```

`cpp/include/fracessa/circular_affine_symmetry.hpp` contains two statically selected implementations:

- `CircularAffineSymmetry` uses precomputed destination bits and fixed arrays for dimensions through 64;
- `CircularAffineSymmetryMultiword` uses precomputed destination positions and reusable multiword scratch storage above 64.

The one-word path walks only set bits. The multiword path extracts positions into reserved storage. Neither path allocates while
filtering an emitted bracelet.

Construction checks exact row-zero entries once, at cost $O(n\varphi(n))$. If only the identity class remains, `fracessa.cpp`
disengages the optional helper and continues through the ordinary circular path.

### Filtering Before Solving

The circular callback rejects noncanonical bracelets before `analyze_support()`:

```cpp
if (symmetry && !symmetry->is_representative(support)) return;
```

The filter stays outside bracelet recursion. Generator ordering, fixed-density recursion, and forbidden-support buckets remain
unchanged.

### Exact Candidate Expansion And Pruning

For an accepted exact candidate $S$, the helper enumerates the distinct bracelets

$$
C_D(aS),\qquad a\in H/\{\pm1\}.
$$

Each bracelet is passed separately to `CircularSupportGenerator::add_forbidden_orbit()`. This preserves both required invariants:

1. every support in the verified affine orbit becomes a future forbidden subset;
2. every public candidate row keeps its ordinary dihedral multiplier, bounded by $2n$.

The implementation does not combine several bracelets into one matrix-specific multiplier. Instead it exactly permutes the solved
probability vector, support, and extended support into each distinct bracelet image and finalizes one row per image. No additional
candidate system is solved.

Registering the complete affine orbit is necessary for dynamic pruning. Filtering by the larger group while registering only the
original dihedral orbit could prune a canonical bracelet without installing the corresponding forbidden support for another affine
image.

## Boundaries

The implementation recognizes only affine index maps

$$
i\mapsto ai+b\pmod n,
\qquad \gcd(a,n)=1.
$$

It does not infer a general colored-graph automorphism group or recognize product-grid groups such as $S_3\times S_8$. For the
published dimension-24 rook matrix, it detects four multiplier classes modulo sign, but not the complete $S_3\times S_8$ group.
More general orbit generation remains an unimplemented research direction in
[`../plans/MAJOR_SINGLE_CORE_PERFORMANCE_OPPORTUNITIES.md`](../plans/MAJOR_SINGLE_CORE_PERFORMANCE_OPPORTUNITIES.md).

## Verification And Evidence

Focused tests cover exact multiplier detection, repeated-value false positives, representative filtering, affine-image expansion,
dihedral multipliers, forbidden-superset pruning, and one-word/multiword agreement across word boundaries.

The promotion audit established:

- byte-identical expanded candidate output with the former dihedral-only path on 200 generated extra-symmetry matrices in both
  `fast` and `safe` mode;
- unchanged stored baselines for all 19 extra-symmetry database matrices;
- 15,120 candidates and 15,120 ESS for the published dimension-24 matrix while reducing solved affine representatives;
- no measurable inactive-path regression: 101 of 120 stored circular matrices had no extra multiplier, and median on/off ratios
  were 1.0000 in both `fast` and `safe` mode; and
- multiword affine decisions and images matching the one-word oracle, including dimensions 65 and 129.

The dated implementation, benchmark, and promotion record remains searchable in `aidocs/CHANGES.md`, entries 292-295 and 378.
