# Generating Fixed-Density Bracelets of Arbitrary Base

**Authors:** S. Karim, Z. Alamgir, and S. M. Husnine  
**Journal:** *International Journal of Computer Mathematics*, Vol. 91, No. 3, 2014, pp. 434–446  
**DOI:** [10.1080/00207160.2013.805753](https://doi.org/10.1080/00207160.2013.805753)  
**Source PDF:** [karim_alamgir_husnine_2014_fixed_density_bracelets.pdf](karim_alamgir_husnine_2014_fixed_density_bracelets.pdf)

> Algorithm-focused extraction and FracESSA implementation notes. This is not a complete transcription. The notation below follows the paper where useful, but the explanations and implementation sketch are restated for this project.

## Result in one sentence

The paper generates fixed-density bracelets directly in lexicographic order and constant amortized time by recursing once per nonzero symbol, carrying the necklace period, and rejecting prefixes as soon as a reversed rotation proves that they cannot be bracelet representatives.

In the FracESSA generator experiments, this 2014 fixed-density algorithm is called **V3**.

For a binary support mask, “density” is exactly support cardinality. The paper is therefore more specialized to FracESSA's circular-support problem than the previously tested fixed-content `BraceletFC` algorithm.

## Objects and notation

A length-​$n$ word over $\{0,\ldots,k-1\}$ has density $d$ when exactly $d$ symbols are nonzero. The paper writes

$$
\mathcal B_k(n,d)
$$

for the set of lexicographically smallest representatives under rotations and reversal.

For FracESSA, $k=2$. A word is a support mask, its ones are supported strategies, and

$$
d=|S|.
$$

Thus the required objects are precisely $\mathcal B_2(n,d)$ for support sizes $d=1,\ldots,n$.

The paper uses the fundamental theorem of necklaces. If a prenecklace $\alpha$ has longest Lyndon prefix length $p$, appending a symbol $b$ preserves the prenecklace property exactly when the repeated-prefix symbol is no larger than $b$. Repeating that symbol preserves $p$; choosing a larger one makes the new word length the next Lyndon-prefix length.

## Why this differs from the current production generator

The production `CircularSupportGenerator` performs a fixed-content FKM recursion one word position at a time. It generates a complete necklace and only then compares it with the smallest rotation of its reflection.

The 2014 algorithm changes both parts:

1. It represents a word by the positions of its nonzero symbols and makes one recursive call per added nonzero symbol.
2. It carries reversal information during recursion and can discard a branch before completing the necklace.

The earlier experimental direct `BraceletFC` generator implements the 2013 **fixed-content** paper. It also performs direct reversal pruning, but it still appends one symbol per recursion. The 2014 fixed-density algorithm is a distinct, more specialized algorithm.

## State used by the paper

The optimized recursion is called `BraceFD(t,p,r,RS)` in the paper.

- $t$: number of nonzero symbols already placed.
- $a_t$: position of the $t$-th nonzero symbol. The positions are strictly increasing.
- $b_{a_t}$: value of that nonzero symbol. In the binary specialization this is always one and need not be stored.
- $p$: density of the longest Lyndon prefix; $a_p$ is its length in word positions.
- $r$: length of the longest prefix that is equal to its reversal and is still relevant to the final bracelet comparison.
- $RS$: result accumulated for the comparison of the suffix after position $r$ with its reversal.
- $s_i$: density of the prefix ending at position $i$. The algorithm only needs this value at nonzero positions.

The complete word array is conceptually initialized to zero. Placing a nonzero symbol changes one position; backtracking restores that position to zero.

For FracESSA, the word itself can be a `uint64_t` mask. The position array still matters because the fixed-density recurrence is expressed in terms of the previous nonzero positions.

## Fixed-density recursion

Assume $t$ nonzero symbols have already been placed. The next recursive call places symbol $t+1$.

The latest position at which it can be placed while leaving enough positions for all remaining symbols is

$$
\operatorname{tail}=n-(d-t)+1.
$$

The repeated-prefix continuation required by the necklace recurrence is

$$
\operatorname{max}=a_{t+1-p}+a_p.
$$

If $\operatorname{max}\leq\operatorname{tail}$, the algorithm first tries that periodic continuation. For an arbitrary alphabet, it tries the repeated symbol and then larger nonzero symbols. In the binary case there is only symbol one, so this becomes one branch.

It then tries every remaining position

$$
j=\operatorname{max}-1,\operatorname{max}-2,\ldots,a_t+1,
$$

or starts at `tail` when the periodic continuation lies beyond `tail`. Choosing one of these earlier positions starts a new Lyndon prefix, so the next period-density parameter is $t+1$.

The important computational consequence is that zero runs are skipped rather than created one zero at a time. Recursion depth is $d$, not $n$.

The final nonzero symbol is not generated recursively. When $t=d-1$, it is placed directly at position $n$, and a constant-time necklace-period test determines whether the completed word is a necklace.

The fixed-density necklace completion uses

$$
\operatorname{next}
=
\left\lfloor\frac{d}{p}\right\rfloor a_p+a_{d\bmod p}.
$$

If `next < n`, the partial word cannot complete a necklace. Equality determines whether the last symbol repeats the prefix or starts a new period. For binary words the last symbol value is fixed to one, so the arbitrary-alphabet value loop disappears.

## Reversal pruning

The paper applies the bracelet criterion incrementally instead of generating all fixed-density necklaces and reflecting every completed one.

Every nonconstant binary necklace begins with a longest zero run. Only reverse rotations beginning with the same zero run can be lexicographically smaller. The algorithm therefore calls its reversal check only when the current prenecklace has the form

$$
0^i\,u\,0^i\,c,
$$

where the first and second zero runs have equal length and the surrounding symbols are nonzero.

The check compares the current prefix with its reversal from the first nonzero position toward the midpoint:

- prefix smaller: keep the branch;
- equal: keep the branch and extend the palindromic-prefix length $r$;
- reversal smaller: discard the branch.

The remaining final comparison concerns the suffix after $r$. In the ordinary bracelet algorithm one new character is appended per recursive call, so one mirrored character can be compared at each step. Fixed-density recursion may jump across several zeros at once. The paper restores constant work by observing that every jump contains only one newly appended nonzero symbol. It compares that symbol with the corresponding mirrored symbol and, if they are equal, compares the two zero-run lengths.

Let

$$
e=n-a_t+r+1.
$$

The zero run beginning at $e+1$ has length

$$
\ell_{e+1}=a_{s_e+1}-a_{s_e}-1.
$$

This lets the algorithm update `RS` without scanning the skipped zeros. The prefix-density metadata $s_i$ is updated only when a nonzero symbol is placed.

## Initialization and edge cases

The paper treats $0<d<n$ in the main algorithm. The all-zero and all-nonzero words are trivial.

The first nonzero position ranges downward from

$$
n-d+1
$$

to

$$
\left\lfloor\frac{n-1}{d}\right\rfloor+1.
$$

The position of the last nonzero symbol is fixed at $a_d=n$. Each initialization then calls the recursive generator with one nonzero symbol already placed.

For its binary CAT proof, the paper assumes $d\leq n/2$. It refers $d>n/2$ to the fixed-content bracelet algorithm. Complementing a binary bracelet gives a bijection between cardinalities $d$ and $n-d$, but the complemented raw word is not necessarily the required canonical representative or in the required order.

This distinction matters to FracESSA. Complement generation is valid for an isolated complete layer after canonicalization, but it is not automatically compatible with dynamic forbidden-subset pruning because complementation reverses inclusion.

## Correctness argument

The algorithm first generates only fixed-density necklaces. It then removes every necklace for which a reversed rotation is lexicographically smaller. Consequently each dihedral orbit contributes exactly one bracelet representative.

The paper states this as: `InitFixed()` lists every member of $\mathcal B_k(n,d)$ exactly once and in lexicographic order.

For a FracESSA implementation, this produces the same raw representative and ascending order required by the existing circular generator, provided the binary position-to-bit convention preserves word lexicographic order.

## Complexity argument

The underlying fixed-density necklace computation tree is constant amortized time. A bracelet contains at most two necklace classes, so restricting that tree to bracelet-producing branches preserves an output-proportional number of recursive nodes.

The only apparently nonconstant operation is the reversal check. The proof maps each equal character comparison injectively to a generated prenecklace. Therefore the total reversal-comparison work over the complete run is proportional to the number of generated bracelets.

The resulting bounds are:

- constant amortized time per bracelet;
- asymptotically linear auxiliary space.

This is an amortized bound for generating a complete fixed-density layer. It does not by itself prove an end-to-end speedup when FracESSA stops early or prunes branches dynamically.

## Binary specialization for FracESSA

The following general-paper machinery disappears for $k=2$:

- no loop over nonzero symbol values;
- no storage for nonzero values, because every nonzero is one;
- no available-symbol linked list;
- no general-content counters;
- no run-length symbol field.

The retained state can be small fixed arrays plus one mask:

- positions of placed ones;
- prefix density at the relevant positions;
- current support mask;
- $t,p,r,RS$.

This is the paper's most promising advantage over the already benchmarked fixed-content `BraceletFC`: it removes both the general alphabet machinery and recursion through explicit zero symbols.

## Required FracESSA adaptation

An isolated generator needs only a callback at each completed bracelet. Production integration additionally requires:

1. support sizes generated in ascending order;
2. representatives emitted in ascending mask order;
3. the callback to receive the support and its cardinality immediately;
4. exact-candidate supports to become pruning rules only at the next cardinality;
5. partial branches to be checked against every stored rotation/reflection of a forbidden support;
6. the full-support configuration and early-empty-layer behavior to remain unchanged;
7. candidate multipliers to retain their current orbit semantics.

Because placed-one positions increase monotonically, a forbidden-subset test can be attached when a new one bit is added. The proof obligation is the same as in the current generator: pruning is valid only once the partial support contains the complete forbidden mask.

## Implementation boundary

The article says that a C implementation was available from the authors on request, but it does not print the program and no downloadable source was found. Joe Sawada's public `necklaces.c` contains the underlying fixed-density **necklace** recurrence, but not this fixed-density bracelet variant.

The project implementation should therefore be treated as an independent adaptation and verified against exhaustive canonicalization, the current FKM-plus-reflection generator, and the previously adapted fixed-content `BraceletFC` generator before any production use.

## Initial assessment

This paper is a credible route to a substantial generator-only speedup because it attacks work that the current and 2013 algorithms both retain: walking through zero positions. It is most attractive for sparse layers, exactly where $d\ll n$.

It is not yet evidence of a full FracESSA speedup. Matrix solving usually dominates, and production pruning changes the recursion tree. The proper next gate is an isolated binary implementation with order-insensitive correctness checks followed by generator-only measurements over complete layers. Production integration should wait for those results.

## V3 experiment status

The isolated binary adaptation is implemented as `V3BraceletGenerator` in `experiments/direct_bracelet_generation_2026-07-29/compare_bracelet_generators.cpp`. It follows the paper's position-of-ones recursion, necklace completion at the final position, incremental reversal test, and delayed zero-run comparison when the mirrored run is not known yet.

V3 matches the existing FKM-plus-reflection generator and the independent fixed-content `BraceletFC` adaptation for every support size through dimension 24. It emits each canonical representative once and in ascending order, and passes optimized plus ASan/UBSan verification. Its generator-only CPU-2 benchmark is recorded in the experiment README and `results/paired_v3_3s.tsv`; production code remains unchanged.
