# Streaming Support Generation With DFS And FKM

Status: deferred architecture and experiment plan. No production implementation
has been approved or made.

## Purpose

FracESSA currently builds support masks eagerly and then physically removes
strict supersets after an exact candidate is found. This document records the
design discussed for replacing that frontier with one-at-a-time generation:

- fixed-cardinality depth-first search (DFS) for general symmetric games;
- fixed-content necklace/bracelet generation for circular-symmetric games;
- candidate pruning inside the generator, so forbidden subtrees are never
  generated;
- no complete `2^n` frontier and, eventually, no complete cardinality layer.

This is a record of the reasoning and experimental evidence, not a decision to
change production. Correctness remains absolute; speed is second; code size and
readability matter after both.

## Short Decision

The promising common architecture is:

```text
for support_size = 1, ..., dimension:
    generate one support or one circular orbit representative
    prune a recursive branch as soon as its partial support contains an
        earlier exact candidate support
    solve each surviving support immediately
    retain only exact candidate masks needed by later cardinalities
```

The two paths differ only in their generator:

- non-circular: a binary fixed-cardinality DFS;
- circular: an FKM-style fixed-content necklace recursion plus bracelet
  reduction, with the same forbidden-subset test inserted into its branches.

A callback from the generator to the existing support solver is likely simpler
than a resumable `get_next_support()` state machine in C++17. Both are streaming;
the implementation should use whichever is shorter and easier to verify.

## Terms

### Support mask

A support is one `uint64_t` mask. Bit `i` means strategy `i` is present. The
analyzer contract is `1 <= dimension < 64`.

### Exact candidate and forbidden support

An exact candidate is a support accepted by `find_candidate_frc()`. Only that
exact result may create a pruning rule. A support accepted only by the double
filter must never prune anything.

Production already calls `remove_supersets()` for every exact candidate,
regardless of whether stability later classifies it as ESS. The invariant to
preserve is:

```text
After exact candidate S is accepted, no strict superset T of S needs to be
searched for an ESS.
```

This does not mean every strict superset is algebraically incapable of being an
equilibrium. It means FracESSA's ESS search and candidate-output semantics
deliberately stop at the earlier exact candidate support. Therefore all exact
candidates, not only ESS supports, become forbidden subsets.

### Rotation class or necklace

For a circular game, cyclic rotations of a support are equivalent. A necklace
class contains all rotations of one support. Its canonical representative is
the smallest binary word, equivalently the smallest mask under the project's
fixed bit ordering.

### Bracelet class

A bracelet class identifies both cyclic rotations and reflections. Its
canonical representative is the smallest member of the full dihedral orbit:
all rotations of the support and all rotations of its reflection.

### DFS

Here DFS means a recursive bitwise construction of a support, not a graph
search over already materialized masks. At each bit the recursion chooses
"exclude" or "include" while tracking how many more selected bits are required.

### FKM

FKM refers to the Fredricksen-Kessler-Maiorana family of necklace-generation
ideas. The experiment discussed here uses the familiar recursive
`Gen(position, period)` form and its fixed-content adaptation: copy the symbol
from the current period, or choose the next larger symbol and start a new
period. For binary supports, the content is fixed by the required number of
ones.

The experimental circular generator is best described precisely as
"fixed-content FKM-style necklace recursion followed by reflection reduction."
It first emits one representative per rotation class, then retains one of the
two reflected rotation classes to obtain bracelets.

## Current Production Behavior

Production `Supports::initialize()` enumerates every nonempty support and stores
it in a vector bucket by cardinality. For circular games, the current code keeps
one rotation representative only when the support cardinality and matrix
dimension are coprime; otherwise it stores every mask.

After an exact candidate is found, `remove_supersets()` scans all larger
cardinality buckets and erases every stored mask containing that candidate.
This has two independent costs:

1. all support masks are built before the first normal support is solved;
2. pruning repeatedly scans and compacts large vectors.

The search order is cardinality first. That property must remain. A candidate
of size `s` can only eliminate strict supersets of sizes `s+1` through `n`; it
cannot eliminate another distinct mask of size `s`.

## Non-Circular Fixed-Cardinality DFS

### Generation order

The successful experiment chooses bits from high to low and visits the exclude
branch before the include branch. This preserves increasing numeric mask order:

```text
generate(bit_count, needed, partial):
    if needed == 0:
        consume(partial)
        return

    bit = bit_count - 1

    if needed < bit_count:
        generate(bit, needed, partial)             // exclude high bit

    with_bit = partial | (1 << bit)
    if with_bit does not complete a forbidden subset:
        generate(bit, needed - 1, with_bit)        // include high bit
```

The actual implementation also needs the usual feasibility bounds: stop when
too few undecided bits remain to reach `needed`, and emit only masks of the
requested cardinality.

### Why one branch test removes a whole subtree

Suppose an earlier exact candidate is mask `S`. While bits are assigned from
high to low, the lowest set bit of `S` is the last required bit of `S` that can
be added. At that moment all higher bits of `S` have already been decided.

If the partial support now contains `S`, every descendant also contains `S`
because recursion only adds bits; it never removes selected bits. The complete
recursive call can therefore be skipped. This is the linked-tree behavior we
wanted: one decision removes a whole implicit subtree without touching every
descendant mask.

### Lowest-bit buckets

The experiment stores each forbidden mask in the bucket indexed by its lowest
set bit. When bit `b` is included, only bucket `b` can contain a forbidden mask
that has become complete for the first time. The hot check is then a contiguous
scan of that bucket:

```text
for S in forbidden_by_lowest[b]:
    if size(S) < current_cardinality and (S & partial) == S:
        prune this branch
```

The `size(S) < current_cardinality` guard states the strict-superset condition.
It also makes it safe to register candidates immediately even though candidates
found in the current layer cannot prune one another.

The one-bit bucket is not merely an arbitrary hash in this traversal: it is the
bit at which containment can first become final. More elaborate indexing may
still reduce scans, but should be added only after measurements show this scan
is a bottleneck.

### Two-bit and deeper indexes

Indexing by two or three selected bits can shrink each bucket, but a query must
inspect every matching pair or triple present in the partial support. It is not
correct to inspect only the two lowest bits of the generated support.

For a rough `n=24`, `k=11` comparison, a one-bit key leaves about `11/24` of
candidates in a bucket, a two-bit key about `C(11,2)/C(24,2)`, and a three-bit
key about `C(11,3)/C(24,3)`. The smaller buckets come with `11`, `55`, or `165`
possible key probes. Contiguous vectors are preferable to linked lists because
the subset scan is simple and cache-sensitive.

Deferred alternatives include a rare-bit anchor, a candidate trie, bit-sliced
candidate IDs, or SIMD scans. None belongs in the first implementation without
profiling evidence.

### No-forbidden fast path

The isolated DFS experiment used a split high-prefix/low-suffix generator while
no forbidden candidates existed, because ordinary combination generation was
then cheaper than recursive branch management. Once the first exact candidate
appeared, it changed to branch-pruning DFS.

This optimization is experimental, not a required part of the architecture. A
first implementation should benchmark the simpler always-DFS form before
retaining a second generator path.

## Gosper Generation And Why It Is Secondary

Gosper's hack generates the next mask with the same number of set bits using a
few arithmetic operations. It naturally supports one-at-a-time generation and
uses constant generator memory.

However, a plain Gosper iterator sees only complete masks. It can skip a
forbidden complete mask or calculate a jump over some contiguous forbidden
range, but it does not expose the binary decision subtree that DFS can abandon
before those masks are constructed.

The experiments established three separate facts:

- streaming Gosper is a strong low-memory baseline;
- candidate-aware Gosper jumps can help strongly pruned matrices;
- a `2^n` byte jump table merely replaces one exponential frontier with
  another and is not the desired final design.

The final non-circular choice must be based on end-to-end benchmarks. DFS is the
cleanest way to integrate pruning into generation; Gosper remains the simpler
fallback if DFS does not produce a material speed advantage.

## Circular Generation: From Necklaces To Bracelets

### Why circular symmetry helps

For a symmetric circulant payoff matrix, rotating all strategy indices only
permutes the selected submatrix. Reflection also preserves it because matrix
symmetry identifies cyclic distance `d` with `n-d`. Therefore one support solve
can represent its complete bracelet orbit.

For `n=8`, support `{0,1,3}` reflects under `i -> -i mod 8` to `{0,7,5}`.
Rotating the reflected set by three positions gives `{0,2,3}`. These masks are
not rotations of each other, but they are reflections and therefore belong to
one bracelet class.

The small generator experiment produced:

| Dimension | Support size | Raw supports | Rotation classes | Bracelet classes |
|---:|---:|---:|---:|---:|
| 8 | 3 | 56 | 7 | 5 |
| 10 | 3 | 120 | 12 | 8 |
| 10 | 5 | 252 | 26 | 16 |

For `n=8`, `k=3`, representative `{0,1,3}` has a dihedral orbit of 16 masks:
eight rotations and eight additional reflected rotations. Orbit sizes can be
smaller when the support itself has rotational or reflection symmetry. The five
bracelet-orbit sizes sum to all 56 raw supports.

### Existing experimental FKM-plus-reflection generator

The retained simple experiment performs these steps for each cardinality `k`:

1. recursively generate fixed-content binary necklaces;
2. at a valid necklace leaf, compare it with the canonical rotation of its
   reflection;
3. emit only the smaller representative, giving one bracelet per dihedral
   orbit;
4. solve that representative once;
5. reconstruct all distinct rotated/reflected candidate rows and ESS counts.

This simple approach was compared with the paper's optimized direct
fixed-content bracelet recursion (`BraceletFC`). Both generated identical orbit
sets through dimension 24 under Release and sanitizers. Direct generation was
faster for the largest generator-only cases but slower for small dimensions;
the end-to-end difference did not justify replacing the simpler
FKM-plus-reflection code.

### Current experimental pruning

The bracelet experiment generates one complete cardinality layer of bracelet
representatives. At each emitted representative it checks whether any stored
forbidden orbit mask is a subset. When an exact candidate is found, all distinct
rotations and reflected rotations are inserted as future forbidden masks.

This already avoids the complete raw support frontier, but it still has two
avoidable pieces:

- the whole bracelet layer is materialized before solving;
- subset pruning occurs at complete leaves rather than inside FKM recursion.

## Integrating Candidate Pruning Into FKM

FKM is already a recursive tree, so the same monotone DFS argument can be used
without generating a complete bracelet first.

The recursion assigns the binary word from its most significant position toward
its least significant position. It tracks at least:

- current position;
- current necklace period;
- number of ones used or still required;
- partial support mask;
- the FKM word buffer needed by the period-copy branch.

Each recursive branch assigns the next symbol. Whenever that symbol is one:

```text
partial |= bit_for_current_position

for S in forbidden_by_lowest[current_bit]:
    if size(S) < target_cardinality and (S & partial) == S:
        do not recurse
```

The test must run both when the FKM period-copy branch copies a one and when an
explicit branch writes a one. The existing fixed-content feasibility checks,
period update, necklace leaf condition, and reflection canonicalization remain
unchanged.

This pruning is correct independently of FKM canonicality. Once a partial mask
contains forbidden `S`, every completion below that recursion node still
contains `S`. FKM may remove additional noncanonical branches, but it cannot
make a selected bit disappear from a descendant.

Candidates found at cardinality `k` cannot prune another distinct mask of the
same cardinality. The forbidden family is therefore logically stable throughout
one FKM layer. The implementation can either activate new masks between
cardinalities or store them immediately and retain the strict-size guard.

## Why One Canonical Forbidden Mask Is Not Enough

The tempting optimization is to store only the canonical bracelet
representative `S` and test canonical generated representative `T` with
`S subset T`. This is incorrect: canonicalization under rotation or reflection
does not preserve subset relations.

### Rotation counterexample

At `n=5`:

```text
S = 00101   canonical mask 5
T = 01011   canonical mask 11
```

The canonical masks do not satisfy `S subset T`. But rotation `01001` of `S`
is a subset of `T`. A canonical-only lookup would solve a support that should
have been pruned.

### Two-necklace counterexample

Storing one canonical rotation for `S` and one for its reflected family is also
insufficient. At `n=6`:

```text
bracelet S              = 001011
stored rotation minima  = 001011 and 001101
canonical T             = 010111
```

Neither stored minimum is a subset of `T`, but reflected rotation `010011` of
`S` is. Thus even two canonical masks do not encode every relevant alignment.

The structural reason is simple: an FKM canonical representative is the first
word in an orbit, but the orientation that embeds into a larger support can be
any other orbit member. The FKM period and reflection state establish
canonicality; they do not provide a subset-preserving map between two orbits.

## Correct Orbit-Pruning Representations

At least one of the following is required.

### Option A: expand orbit masks

Store every distinct rotation and reflected rotation of each exact candidate,
bucketed by lowest set bit. This uses at most `2*n` 64-bit masks per bracelet
candidate and makes each branch check one ordinary mask containment operation.

This is the recommended first experiment because it is the shortest and easiest
to prove. The apparent duplication is modest at the tested dimensions. Matrix
34 has 15,120 output candidates, all at support size 10; storing one raw
`uint64_t` for each is 120,960 bytes, about 118 KiB, before vector capacity and
bucket overhead.

### Option B: test all rotations in parallel from canonical masks

Keep one necklace mask, plus a second one when reflection creates another
rotation class. For a partial support `P`, maintain a bitset of still-valid
cyclic shifts. For every set bit of candidate `S`, rotate or align `P` and
intersect the valid-shift set. A nonzero result means some rotation of `S` is
contained in `P`.

This stores only one or two masks per candidate orbit, but each query costs
roughly `popcount(S)` rotations and intersections per candidate representative.
It may move more work into the hot branch test than explicit orbit expansion.
It should be tested only as a separate benchmark against Option A.

### More elaborate representations

A trie, zero-suppressed decision diagram, compiled automaton, or bit-sliced
candidate index could share orbit structure. These add substantial code and
state. They are not justified unless expanded masks become a measured memory or
scan bottleneck. No exact published algorithm for this specific dynamic
forbidden-subset, fixed-content, dihedral-orbit problem was identified in the
limited literature search.

## Alternative: Complete-Support Reverse Lookup

There is a simpler two-step alternative that does not alter FKM recursion:

1. generate a complete bracelet representative `T`;
2. enumerate subsets of `T` at candidate cardinalities and look them up in an
   exact forbidden set.

For each candidate size `s`, choose the cheaper direction:

```text
min(number of stored candidates of size s, C(popcount(T), s))
```

This can reduce leaf-level subset comparisons, especially when candidate
families are large. It cannot prune internal FKM nodes, so it still visits every
canonical leaf. It is a secondary fallback, not the primary integrated design.

## Retained State In The Target Design

The search should need only:

- current cardinality and recursion state;
- current partial support mask;
- exact forbidden candidate masks and their cardinalities;
- matrix and solver scratch buffers;
- full candidate objects only when candidate output is requested.

The forbidden set and returned candidate rows are separate. ESS-count-only
execution requires support masks for pruning but not stored vectors, payoffs, or
stability strings.

Do not precompute or serialize an exponential support tree. The static generator
state is cheap; the useful forbidden family is matrix-specific and is discovered
during the solve. Loading a `2^n` table cannot remove that dynamic work.

## Correctness And Output Contract

An intentional generator replacement may change deterministic enumeration
metadata:

- candidate row order;
- `candidate_id`;
- `shift_reference`.

Those fields remain useful regression checks while the generator is unchanged,
but they are not mathematical correctness contracts across a deliberate
enumeration redesign. A new deterministic order and regenerated baseline are
acceptable only after an independent order-insensitive comparison proves that
the following remain identical:

- the complete FracESSA candidate set under the existing pruning semantics;
- candidate count and ESS count;
- each support and strategy vector;
- exact payoff and stability result;
- each candidate's ESS classification.

For a circular representative, candidate reconstruction must include every
distinct orbit member exactly once. If an orientation `g(T)` contains candidate
`S`, then canonical `T` contains `g^-1(S)`. Registering the complete dihedral
orbit of `S` therefore makes the pruning test orientation-independent.

## Experimental Evidence Already Available

The full experiment report and raw results are retained on branch
`choice-one-candidate-filter` in
`aidocs/experiments/SUPPORT_FRONTIER_2026-07-29.md` and the corresponding
`experiments/` directories.

### Streaming Gosper

- Complete candidate output was byte-identical for all 52 matrices tested.
- Removing the one-layer vector changed summed medians by `+0.24%` and was
  `1.02%` slower by geometric mean on dimensions 15-24.
- In a dimension-24 no-pruning control, peak RSS fell from `34,576 KiB` to
  `13,632 KiB`, while elapsed time remained `9.37 s` versus `9.39 s`.

Conclusion: streaming itself is a memory improvement, not a demonstrated speed
optimization.

### Non-circular DFS

Pure fixed-cardinality DFS produced byte-identical candidate output on the
tested non-circular corpus. Compared with production, it improved large IDs 49,
51, 53, and 54 by `40.68%`, `67.50%`, `33.19%`, and `14.52%`, but regressed ID
52 by `2.67%`. On the dimension-20 no-pruning full-support control, production,
Gosper, and DFS were effectively equal.

Conclusion: integrated branch pruning can be valuable, but no dimension-only
gate predicts it. The discovered candidate family is matrix-specific.

### Circular bracelet generation

- The standalone generator found 29 nonempty bracelet classes at `n=8`, 77 at
  `n=10`, and 352,697 at `n=24`, representing all `2^n-1` nonempty supports.
- The orbit-aware bracelet analyzer matched mathematical candidate/ESS results
  on all 27 circular fixtures and 104 generated circular games; the 24
  non-circular fixtures remained byte-identical; sanitizers passed.
- It improved 24 of 27 circular fixtures with a `4.28x` geometric-mean speedup.
- Matrix 34 improved from `22.132 s` to `0.650 s`, while peak RSS fell from
  `99,520 KiB` to `13,008 KiB`.
- Direct `BraceletFC` generation was `17.34%` slower by geometric mean at
  dimensions 8-14 but `12.19%` faster at dimensions 19-24, including a `22.31%`
  generator-only win at dimension 24.

Conclusion: circular orbit reduction is the strongest measured opportunity.
The simpler FKM-plus-reflection generator remains the preferred base until an
end-to-end result justifies the more complex direct bracelet recursion.

## Implementation Experiment, When Approved

No production change should be made directly. Use isolated experiment sources
and proceed in small stages so each effect is measurable.

### Stage 1: branch pruning in the existing circular generator

- keep the existing FKM-plus-reflection algorithm;
- keep expanded dihedral forbidden masks;
- insert the lowest-bit subset test into both one-writing recursion branches;
- retain the current bracelet layer temporarily, isolating only the effect of
  moving pruning from leaves into recursion.

### Stage 2: remove the circular layer

- replace layer accumulation with a callback that immediately solves each
  surviving representative;
- retain the same generation, pruning, and orbit reconstruction logic;
- compare memory and end-to-end time with Stage 1.

### Stage 3: unify the non-circular path

- expose the same callback contract from fixed-cardinality binary DFS;
- preserve cardinality-first traversal and increasing numeric order where
  practical;
- retain Gosper as the experimental baseline, not as another permanent mode
  unless measurements require it.

### Optional Stage 4: compact orbit test

- benchmark the one/two-canonical-mask bit-parallel rotation test against
  expanded orbit masks;
- keep it only if end-to-end time or memory improves materially and the proof
  remains understandable.

## Measurements Required

Instrument only the experiment, not production. Record:

- recursive nodes entered;
- branches rejected by cardinality bounds;
- branches rejected by forbidden subsets;
- forbidden-mask comparisons;
- complete necklace leaves;
- reflection rejections;
- bracelet representatives yielded;
- support solves;
- candidate representatives and expanded orbit-mask counts;
- peak RSS and end-to-end median time.

Correctness checks must include:

- all maintained circular and non-circular verification matrices;
- randomized circular games for manageable dimensions, for example `2..14`;
- exact order-insensitive candidate and ESS comparisons;
- explicit orbit uniqueness and complete-orbit reconstruction checks;
- ASan and UBSan runs for the isolated generator/analyzer.

Performance should use the existing persistent-process, pinned-CPU,
approximately three-second median protocol. Include small, medium, and large
cases, matrix 34, strongly pruned non-circular games, and the no-pruning
full-support control.

## Acceptance Rule

Production adoption requires all of the following:

1. zero mathematical result differences;
2. a meaningful end-to-end speed benefit on the intended matrix classes;
3. lower memory without replacing the frontier by another exponential table;
4. a short implementation whose pruning proof is readable in the source;
5. no permanent matrix-property switch unless a measured, reliable predicate
   justifies it.

If integrated branch checks regress weakly pruned games and no simple measured
fast path fixes that, keep the current production path. A sophisticated
generator is not worthwhile for a small or theoretical improvement.

## Literature Boundary

The following sources support the general recursive-generation approach:

- [Ruskey and Sawada, An Efficient Algorithm for Generating Necklaces with
  Fixed Density](https://www.socs.uoguelph.ca/~sawada/papers/fix.pdf): recursive
  fixed-density necklace generation and the suitability of recursion for
  backtracking and subtree pruning.
- [Karim, Sawada, Alamgir, and Husnine, Generating Bracelets with Fixed
  Content](https://www.socs.uoguelph.ca/~sawada/papers/fix-brace.pdf): direct
  fixed-content bracelet recursion and its optimized implementation.
- [Ruskey and Sawada, Generating Necklaces and Strings with Forbidden
  Substrings](https://webhome.cs.uvic.ca/~ruskey/Publications/Forbidden/Forbidden.html):
  evidence that constraints can be carried inside necklace recursion rather
  than applied only after generation.

These papers do not directly solve FracESSA's dynamic condition: forbidden
supports are discovered during computation, the condition is subset
containment rather than a fixed forbidden substring, and circular candidates
must be considered under rotations and reflections. The proposed integration
and its proof therefore remain project-specific and must be tested as such.
