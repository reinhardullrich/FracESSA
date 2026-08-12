# Major Single-Core Performance Opportunities

Status: current research backlog. The integer exact kernel and shared candidate/stability factorization proposed by the original
review have been implemented; the remaining opportunities are unimplemented hypotheses that require measurement before adoption.

## Scope

The target is a material reduction in the native time for one matrix on one CPU
core. Python already provides process-level parallelism, so C++ threads, GPU
execution, and additional process parallelism are explicitly out of scope.
Micro-optimizations that plausibly top out at a few percent are also out of
scope.

The important mathematical observation is that FracESSA accepts only symmetric
payoff matrices. The parser expands either an upper triangle or a
circular-symmetric compact row into a symmetric matrix. Consequently, for

$$
q(x)=\frac12 x^T A x,\qquad
\Delta=\{x\geq0: \mathbf 1^T x=1\},
$$

the exact candidate conditions are precisely the KKT conditions for maximizing
$q$ over $\Delta$:

$$
(Ax)_i=u\quad(i\in S),\qquad
(Ax)_i\leq u\quad(i\notin S).
$$

The current program is therefore an exhaustive active-set method for a quadratic program: it enumerates a support $S$, solves its KKT
equations, and classifies stability, reusing the candidate factorization whenever possible. The largest remaining improvements must
either examine far fewer active sets or avoid material exact work that profiling still shows to be expensive.

## Evidence From The Retained Workload

These are historical measurements from the 87-matrix snapshot used for the original review, not current whole-database results or new
tests:

- Across that 87-matrix timing session, the then-used verified, unsafe, and exact
  totals are respectively `347.148 s`, `59.389 s`, and `1437.953 s`.
- Five structured many-candidate matrices (IDs 66, 65, 90, 64, and 89) account
  for `291.510 s`, or about 84% of the complete verified total.
- Of 86,152 represented exact candidates in the 87-matrix timing corpus, 85,249
  (98.95%) have extended support equal to support. Thus almost every retained
  candidate reaches the common easy stability case with no unused tied best
  reply.
- The 5.8-fold verified-to-unsafe total gap was not permission to use the unsafe result and did not prove that the gap was attainable.
  It quantified how much single-core time the rigorous per-support path added at that time.

That snapshot suggested two different strategies: structure-aware algorithms for the dominant graph-like families and a broadly faster
exact kernel for unstructured games.

## Completed Foundations

The following proposals from the original review are now production behavior:

- The parser clears denominators once and the exact candidate and stability kernels use FLINT integers rather than rational arithmetic.
- One reduced symmetric candidate system of dimension $k-1$ replaces the bordered system of dimension $k+1$.
- The fraction-free $LDL^T$ factorization returns the exact solution and inertia. Stability reuses that factorization and constructs
  only the smaller scaled reduced $B$ matrix through a Schur complement when outside best replies leave the result unresolved.

The implementation details belong in `../PROJECT.md` and `../reference/EXACT_STABILITY_SCHUR_COMPLEMENT.md`; they are not repeated
as open proposals here.

## Open Heavy Hitters

The rank reflects potential against the retained historical snapshot, not implementation
order or a promise of measured speedup. “Very high” means that the mechanism
can remove an exponential search or thousands of repeated factorizations.
“High” means that it can remove a whole factorization or replace rational
arithmetic throughout a major stage.

| Rank | Opportunity | Plausible reach | Coverage | Risk |
|---:|:---|:---|:---|:---|
| 1 | Recognize graph-payoff families and enumerate maximal cliques | Very high | Structured graph games | Medium |
| 2 | Detect strictly concave games and solve their unique equilibrium directly | Very high | Mathematically qualifying games | Medium |
| 3 | Quotient support search by exact matrix automorphisms | Very high | Symmetric/repeated-block games | Medium to high |
| 4 | Reuse linear algebra across neighboring supports | High | Broad exponential scans | High |
| 5 | Certify positive definiteness numerically before exact stability | High | Candidate-heavy games | Medium |
| 6 | Probe promising full support automatically | Very high when it succeeds | Full-support ESS games | Low |
| 7 | Eliminate never-best-reply strategies before enumeration | Exponential in strategies removed | Reducible games | Low to medium |
| 8 | Replace support enumeration by complementarity reverse search | Potentially very high | General nondegenerate cases | Very high |

### 1. Exact graph-family dispatch to maximal-clique enumeration

Several of the slowest retained matrices are not arbitrary dense matrices.
The complete-multipartite families and published `G(...)` families encode
graph quadratic programs, and their thousands of ESS supports are combinatorial
objects. FracESSA currently rediscovers each object by solving many unrelated
linear systems.

For a theorem-approved graph payoff family, recognize exact affine equivalence
to a regularized adjacency game, then enumerate the corresponding maximal
cliques directly with a bitset branch-and-bound algorithm such as Bron--Kerbosch with pivoting. The existing support representation
keeps the one-word fast path through dimension 64 and uses multiple words above it. Affine recognition must be exact because adding a
common payoff and multiplying by a positive scalar preserve candidates and
ESS, while approximate pattern recognition would not.

The payoff is potentially enormous: a bitset clique generator visits
combinatorially feasible supports, not every subset followed by a KKT solve.
For a complete multipartite game, it can work directly with the parts and emit
one choice from every part. That replaces millions of numerical proofs and
thousands of exact stability tests with a short combinatorial recursion.

Correctness boundaries:

- Enable the path only for a precisely derived payoff family. “The entries look
  graph-like” is not a proof.
- If only ESS output is required, the established strict-local-maximum/maximal-
  clique equivalence can suffice. If all Nash candidates are requested, either
  prove the stronger candidate equivalence for the detected subclass or retain
  the general candidate path.
- During development, compare the complete emitted set with the ordinary exact
  engine, not merely each emitted support. At runtime, fall back if recognition
  or a detector invariant fails.

This is the first workload-specific experiment because graph-like instances, not random matrices, dominated the historical snapshot.
Bomze's regularized maximum-clique formulation is the relevant mathematical
starting point; maximal-clique enumeration itself can use the classic
[Bron--Kerbosch algorithm](https://doi.org/10.1145/362342.362367).

### 2. Strict-concavity fast path: one equilibrium, no support enumeration

Choose any full-rank tangent basis, for example columns $e_i-e_m$, and test

$$
Z^TAZ\prec0.
$$

If this holds, $q(x)$ is strictly concave on the complete simplex. Its KKT point
is the unique symmetric Nash equilibrium and the unique ESS. FracESSA does not
need to enumerate any support at all.

The algorithm for this class is:

1. certify global tangent negative definiteness once;
2. solve the strictly concave simplex quadratic program with an active-set
   method;
3. use the returned support only as a proposal;
4. verify that support with the exact KKT and outside inequalities;
5. return the unique exact result.

The numerical optimizer does not need to be trusted. Strict concavity proves
uniqueness, and the final exact verification proves the proposed answer. If the
proposal fails, the normal exhaustive engine remains available.

This turns an exponential search into a polynomial-size active-set solve for
every strictly concave input, including boundary equilibria that an early
full-support test alone cannot find. In game-theory terminology these are
strictly stable games; negative definiteness on the simplex tangent space is
the key condition.

### 3. Generate support orbits under the actual automorphism group

The circular generator uses a known dihedral symmetry group, but a general
symmetric matrix may have a much larger automorphism group

$$
\mathrm{Aut}(A)=
\{\pi:A_{\pi(i),\pi(j)}=A_{ij}\}.
$$

Candidates and ESS are carried to candidates and ESS by every such
permutation. Therefore only one support per group orbit needs expensive
analysis; the remaining outputs are permutations of the representative.

The simplest useful version is exact twin/block detection. Strategies with the
same interaction pattern form interchangeable classes, and supports can be
generated as count vectors across classes. This directly attacks complete-
multipartite and repeated-block matrices without a general group library.

Only if that leaves substantial measured time should FracESSA encode the
rational matrix as an edge-coloured graph, compute its full automorphism group,
and use canonical augmentation to emit supports isomorph-free. McKay's
[isomorph-free exhaustive generation](https://users.cecs.anu.edu.au/~bdm/papers/orderly.pdf)
is the appropriate framework.

Canonicalizing every support after ordinary generation is not enough: it saves
solves but still visits $2^n$ masks. The generator itself must avoid noncanonical
branches. Candidate output can either be expanded from the orbit or represented
with a generalized multiplier; the latter would be an API decision, not merely
an optimization.

### 4. Reuse factorization state across the subset lattice

Candidate supports overlap heavily, but the current code factors each reduced system from scratch in $O(k^3)$ work. Reuse is not a
drop-in update: the reduced matrix depends on the chosen reference strategy, while the larger bordered KKT systems have the simpler
principal-submatrix relationship.

Two related redesigns can lower the per-support algebra:

- Traverse fixed-cardinality supports in revolving-door order. Neighboring
  supports exchange one strategy, so their KKT matrices differ by one row and
  one column. Update and downdate a factorization in $O(k^2)$, with a complete
  refactorization when the update is singular or fails its verification.
- Traverse a subset tree while carrying a Schur-complement state from parent to
  child. Fast principal-minor algorithms exploit exactly this relation and can
  compute all principal minors in order $O(2^n)$ arithmetic operations rather
  than independently paying a cubic factor. See Griffin and Tsatsomeros,
  [*Principal minors, Part I*](https://doi.org/10.1016/j.laa.2006.04.008).

FracESSA needs more than determinants: it needs the last inverse column for
$(x,u)$ and exact outside inequalities. Nevertheless, the same block inverse or
Schur state can update that solution in quadratic rather than cubic work.

This is a credible general improvement, but it is a high-risk implementation. Singular intermediate principal submatrices, exact
integer growth, stable verified updates, reference changes, and support-pruning order all need explicit fallbacks. Attempt it only if
profiling of the current fraction-free kernel still shows support-by-support factorization as a dominant cost.

### 5. One-sided ball Cholesky before exact reduced-B classification

After exact candidate construction has established the exact extended support,
try to prove the scaled reduced $B$ matrix positive definite with rigorous ball arithmetic. A successful ball Cholesky is
a complete proof and can skip its exact integer positive-definiteness factorization. An inconclusive result falls through to the
unchanged exact reduced-$B$ classification.

FLINT already ships Arb, and
[`arb_mat_cho`](https://flintlib.org/doc/arb_mat.html) returns success only when
the represented symmetric matrix is certainly positive definite. This makes
Arb suitable as a narrow proof filter, not as a generic replacement for every
small double solve.

This path is attractive because exact positive definiteness is overwhelmingly
the successful stability result in retained data. It can also certify the one
global matrix required by the strict-concavity fast path. It does not justify
optimizing full copositivity first: the retained workload reaches successful
full copositivity for only a handful of candidate representatives.

### 6. Automatic, cheap full-support probe

The existing `--fullsupport` path can turn a full-support ESS game from an
exponential search into one candidate and stability test. A completely mixed
ESS is unique, so returning immediately is valid; the current implementation
already relies on that fact when the flag is selected.

The default engine can obtain most of this benefit without always paying for an
expensive exact full-size solve:

1. run an ordinary full-support KKT and curvature probe;
2. only if both look clearly promising, run the rigorous/exact confirmation;
3. return on a proven ESS;
4. otherwise start normal cardinality enumeration.

The preliminary result controls scheduling only, never correctness. This is an
orders-of-magnitude win for the retained full-support examples and a small
bounded overhead elsewhere. It is also far simpler than a new enumeration
engine, so it should be tested early even though its coverage is narrow.

### 7. Eliminate strategies that can never be best replies

Before support generation, prove that a strategy is never a best response to
any point of the simplex. A purely strictly dominated row is the cheapest case;
strict domination by a mixture of other rows is the stronger LP test.

Such a strategy can never belong to a candidate support or extended support.
Remove it from the support universe while retaining its exact outside-payoff
check for final verification. Removing $r$ strategies divides the raw support
space by $2^r$, so even one or two successful eliminations can exceed the
requested 30% threshold.

Use only strict, exactly certified elimination. Weak dominance is unsafe for
equilibrium enumeration, and a floating LP status alone is not a proof. Start
with the short $O(n^3)$ pure-dominance scan; add mixed-dominance certificates
only if real matrices actually benefit.

### 8. Output-sensitive complementarity or reverse-search enumeration

The candidate equations can be written as the complementarity system

$$
x\geq0,\qquad u\mathbf1-Ax\geq0,\qquad
x_i\bigl(u-(Ax)_i\bigr)=0,\qquad \mathbf1^Tx=1.
$$

The current algorithm guesses every active set and asks whether it solves this
system. A reverse-search or complementary-pivot algorithm instead walks from
one feasible basis to adjacent feasible bases and can enumerate candidate KKT
points rather than all $2^n$ supports. Exact basis arithmetic preserves the
present correctness goal.

This is the most general route to an output-sensitive algorithm and could be
orders of magnitude faster when candidates are sparse. It is also the largest
research project. Degeneracy, singular candidate systems, duplicate bases,
proof that every required symmetric candidate is reached, deterministic
output, and interaction with current superset pruning all need treatment.
The lrsNash work of Avis, Rosenberg, Savani, and von Stengel on
[reverse-search equilibrium enumeration](https://doi.org/10.1007/s00199-009-0449-x)
is the closest established model, but it is not a drop-in implementation of
FracESSA's symmetric candidate contract.

## What Is Unlikely To Reach The Target

The following may still be reasonable maintenance work, but they are not
30--50% strategies for the current program:

- saving a few matrix or vector allocations;
- another container reserve policy;
- replacing the tiny matrices with a generic BLAS/Eigen call;
- further bit-scanning or support-generator micro-tuning;
- generic MPFR/Arb arithmetic for every support;
- FFT diagonalization of a circular matrix, because arbitrary principal
  support submatrices are not circulant;
- optimizing exact full copositivity before a workload shows that it dominates;
- compiler flag, inlining, loop-unrolling, or logging cleanup presented as the
  main performance plan.

The allocation audit already shows the general lesson: dramatically fewer allocator calls did not produce a stable application-level
speedup. The next round should change the number of supports examined or remove a measured mathematical kernel, not polish memory
ownership again.

## Recommended Research Sequence

1. Profile the current build on the balanced matrix panel before selecting another project; the historical 87-matrix totals are not a
   current bottleneck measurement.
2. If graph-like families still dominate, add the narrowest theorem-backed exact recognizer and compare direct maximal-clique output
   with the complete exact engine.
3. Test the global strict-concavity and automatic full-support paths as separate small experiments.
4. Add one-sided Arb positive-definiteness only if exact scaled-reduced-$B$ classification remains material.
5. Try exact twin-class orbit generation before considering a full automorphism package.
6. Attempt subset-lattice updates or complementarity reverse search only if simpler changes leave support-by-support factorization as
   the measured bottleneck.

Each experiment should be retained only after the established balanced panel
covers small and large dimensions, circular and non-circular inputs, and
candidate-sparse and candidate-heavy games. Correctness comparison remains
exact and order-independent. A large kernel win without a large end-to-end win
is not sufficient.
