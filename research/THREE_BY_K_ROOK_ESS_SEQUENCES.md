# Three-by-$k$ Rook ESS Sequences

## Purpose and status

This note records two explicitly verified finite ESS sequences, their proposed infinite extensions, and the unresolved third
congruence class of symmetric games on a $3\times k$ rook grid. They arose from the published $3\times8$ game with 15,120
evolutionarily stable strategies (ESS), and from the later integer search over three rook payoffs.

The main conclusions are:

1. A selected strategy is an edge of the complete bipartite graph $K_{3,k}$.
2. In the parameter region used here, a support containing a cycle cannot be an ESS support.
3. Every verified support constructed here is a spanning tree with $k+2$ selected strategies.
4. Compact circular encoding of all three rook relations exists precisely when $3\nmid k$, giving two interleaved sequences:
   $k=3s+1$ and $k=3s+2$.
5. The three verified support orbits in those two sequences have exact combinatorial counts.
6. Exact rational candidate and stability checks establish the Sequence I and II counts as ESS lower bounds through dimension 60.
7. Complete FracESSA searches establish equality with the total ESS count only for the four smallest cases checked exhaustively:
   $3\times4$, $3\times5$, $3\times7$, and $3\times8$.
8. The remaining class $k=3s$ has three natural balanced spanning-tree orbits with an exact support count, but no tested payoff
   triple makes all three orbits ESS. Its displayed values are therefore support counts, not ESS counts.

For Sequences I and II, the formulas below count explicitly constructed ESS orbits. For $k=3s$, they count only support orbits.
None of the formulas proves that a game contains no additional ESS, nor that a construction maximizes the number of ESS among
rook games or arbitrary symmetric games of the same dimension.

The published matrix and 15,120 count come from the extracted
[`bomze_schachinger_ullrich_2018_complexity_of_simple_models.md`](papers/bomze_schachinger_ullrich_2018_complexity_of_simple_models.md).
The earlier symmetry reconstruction is in [`AUTOMORPHISM_GROUPS_AND_GAMMA.md`](AUTOMORPHISM_GROUPS_AND_GAMMA.md); the broader
finite rook-group catalogue is in [`ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md`](ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md). The sequence
formulas and finite exact checks are the new results recorded here.

## From rook automorphisms to the $A,B,C$ family

### The rook board

Index the $N=3k$ pure strategies by cells $(r,c)$ of a $3\times k$ grid, with

$$
r\in\{0,1,2\},\qquad c\in\{0,\ldots,k-1\}.
$$

Equivalently, identify a strategy with an edge of the complete bipartite graph $K_{3,k}$: the cell $(r,c)$ is the edge joining
row vertex $r$ to column vertex $c$. Two cells share a row or column exactly when the corresponding graph edges meet.

Every permutation of the three rows and every independent permutation of the $k$ columns relabels the strategies while
preserving these incidences. Hence the rook action is

$$
G=S_3\times S_k,
\qquad
|G|=3!\,k!=6k!.
$$

For $k\ne3$, the row cliques have size $k$ and the column cliques have size three, so the two kinds are distinguishable. In
particular, the ordinary two-value rook graph has full automorphism group $S_3\times S_k$ when its two values differ. For the
square $3\times3$ board on the two-value slice, exchanging rows and columns supplies an additional $C_2$ symmetry.

The rook action has four orbits on ordered pairs of strategies:

1. the same cell;
2. two different cells in the same row;
3. two different cells in the same column;
4. two cells sharing neither coordinate.

A symmetric zero-diagonal matrix invariant under $S_3\times S_k$ must therefore be constant on the last three orbits. This gives
the three-parameter family $M(A,B,C)$ with diagonal zero and three off-diagonal values:

$$
M_{(r,c),(r',c')}=
\begin{cases}
0, & (r,c)=(r',c'),\\
A, & r=r'\text{ and }c\ne c',\\
B, & c=c'\text{ and }r\ne r',\\
C, & r\ne r'\text{ and }c\ne c'.
\end{cases}
$$

Thus $A$ is the same-row payoff, $B$ is the same-column payoff, and $C$ is the payoff for cells sharing neither coordinate.

Conversely, every matrix of this form has the rook group $S_3\times S_k$ as automorphisms. If the three values are distinct, this
is the full automorphism group, including at $k=3$ because the distinct values distinguish rows from columns. On the two-value
slice $A=B\ne C$, it is the full group for $k\ne3$; at $k=3$, transposition extends it to
$(S_3\times S_3)\rtimes C_2$. Other equalities can merge pair orbits and enlarge the group, but never destroy the rook
symmetries. If $A=B=C$, for example, every permutation of the $3k$ strategies is an automorphism.

### The old two-value family is a slice of the $A,B,C$ family

The original rook construction used only two off-diagonal values:

$$
P=\text{payoff for cells sharing a row or column},
\qquad
Q=\text{payoff for cells sharing neither coordinate}.
$$

In the present notation this is exactly

$$
\boxed{A=B=P,\qquad C=Q.}
$$

The published $3\times8$ game used $(P,Q)=(7,15)$, hence $(A,B,C)=(7,7,15)$. Allowing $A$ and $B$ to differ does not sacrifice
the $S_3\times S_k$ automorphisms: row permutations and column permutations still preserve all three pair relations. It merely
allows the numerical game to distinguish same-row from same-column interaction. The later search over the full three-variable
family found, for example, $(A,B,C)=(6,7,14)$ for the same $3\times8$ ESS orbit and $(8,10,19)$ for the new $3\times10$
construction.

This explains why the three-variable search is not a different symmetry construction. It is the complete zero-diagonal matrix
family supported by the same rook action; the two-variable search had examined only its plane $A=B$. For $k\ne3$, separating
$A$ and $B$ neither enlarges the automorphism group nor makes circular encoding possible in a new case. Its only benefit is the
additional numerical freedom to tune row and column interactions independently.

### Why automorphisms produce many ESS

If a support $I$ is a candidate or an ESS, every transformed support $gI$, $g\in G$, has the same property. Its number of distinct
copies is determined by orbit-stabilizer:

$$
|G\cdot I|=\frac{|G|}{|\operatorname{Stab}_G(I)|}.
$$

A large group therefore produces a large ESS orbit only when the support has a small stabilizer. Group order alone does not imply
many ESS. The published $3\times8$ support has stabilizer order 16, so its orbit has

$$
\frac{3!\,8!}{16}=15{,}120
$$

members. The constructions below count their support orbits directly; orbit-stabilizer gives the same totals.

## The $A,B,C$ quadratic form

For an arbitrary vector $y=(y_{rc})$, define its total, row sums, and column sums by

$$
T=\sum_{r,c}y_{rc},
\qquad
R_r=\sum_c y_{rc},
\qquad
L_c=\sum_r y_{rc}.
$$

The three off-diagonal pair contributions can be written without enumerating pairs:

$$
\begin{aligned}
\text{same row}
&=\sum_rR_r^2-\sum_{r,c}y_{rc}^2,\\
\text{same column}
&=\sum_cL_c^2-\sum_{r,c}y_{rc}^2,\\
\text{neither}
&=T^2-\sum_rR_r^2-\sum_cL_c^2+\sum_{r,c}y_{rc}^2.
\end{aligned}
$$

Multiplying these expressions by $A$, $B$, and $C$ gives the homogeneous identity

$$
\boxed{
y^TMy
=CT^2+(A-C)\sum_rR_r^2+(B-C)\sum_cL_c^2+(C-A-B)\sum_{r,c}y_{rc}^2
}.
$$

For a mixed strategy $x$, $T=1$, so this reduces to

$$
x^TMx
=C+(A-C)\sum_rR_r^2+(B-C)\sum_cL_c^2+(C-A-B)\sum_{r,c}x_{rc}^2.
$$

The coefficient

$$
\delta=C-A-B
$$

controls perturbations that preserve every row and column mass.

## Exact scale normalization

Multiplying every payoff by a positive scalar preserves candidates, supports, ESS, and every payoff comparison. If $\delta>0$,
divide the matrix by $\delta$. The normalized parameters satisfy

$$
\widehat C-\widehat A-\widehat B=1,
$$

or equivalently

$$
\boxed{\widehat C=\widehat A+\widehat B+1}.
$$

This explains why every record triple found in the present rook search lies on the plane $C=A+B+1$. It is an exact
normalization of the complete $\delta>0$ parameter region, not an empirical equation specific to one dimension.

For $\delta\ne0$, positive scaling by $|\delta|$ gives the first two normalized relations below. The zero region is already
scale-invariant:

$$
\begin{array}{c|c}
\delta>0 & \widehat C=\widehat A+\widehat B+1\\
\delta<0 & \widehat C=\widehat A+\widehat B-1\\
\delta=0 & C=A+B.
\end{array}
$$

The normalized values $\widehat A$ and $\widehat B$ may be rational. Dropping hats after normalization, restricting the plane
$C=A+B+1$ to integer $A,B$ is therefore not an exhaustive parametrization of all positive-scale rays. It is, however, a useful
exact coordinate system for analysis and adaptive search.

The normalization does not imply that the individual payoffs are positive. For example, $(A,B,C)=(1,-1,1)$ also has
$\delta=1$.

## Compact circular representation

Map an index $z\in\mathbb Z_{3k}$ to the grid cell

$$
z\longmapsto (z\bmod3,z\bmod k).
$$

The Chinese remainder theorem makes this a bijection exactly when

$$
\gcd(3,k)=1.
$$

For three distinct values, and also for the two-value slice $A=B\ne C$, this condition is necessary rather than merely a
limitation of the displayed indexing. A circular representation would supply a permutation acting transitively as one cycle on
all $3k$ cells. Every matrix automorphism lies in $S_3\times S_k$ for $k\ne3$. A transitive element must project to a
three-cycle on the rows and a $k$-cycle on the columns, so its cell orbit has length $\operatorname{lcm}(3,k)$. This length equals
$3k$ exactly when $\gcd(3,k)=1$.

The exceptional square case does not evade the obstruction. At $k=3$, elements of $S_3\times S_3$ have order at most six.
Writing transposition as $\tau$, an element $(\sigma,\rho)\tau$ has square
$(\sigma\rho,\rho\sigma)$; its two components are conjugate elements of $S_3$, so this element also has order at most six.
Therefore the enlarged two-value group contains no nine-cycle on the cells.

Under this indexing, the entries at a nonzero circular offset $d$ are determined directly:

$$
a_d=
\begin{cases}
A, & 3\mid d,\\
B, & k\mid d,\\
C, & 3\nmid d\text{ and }k\nmid d.
\end{cases}
$$

The first two cases cannot occur simultaneously for a nonzero offset modulo $3k$. Thus the compact symbolic strings for the two
new examples are

```text
3 x 10: 30#C,C,A,C,C,A,C,C,A,B,C,A,C,C,A
3 x 11: 33#C,C,A,C,C,A,C,C,A,C,B,A,C,C,A,C
```

Substituting $(A,B,C)=(8,10,19)$ produces the numerical inputs shown later. The circular representation exposes translations and
reflections, but those form only a small subgroup of the complete $S_3\times S_k$ rook symmetry.

Consequently, the symbolic three-relation $3\times k$ rook pattern has the compact circular representation used by FracESSA
precisely when $k$ is not a multiple of three. Degenerate value equalities can produce simpler matrices with additional
representations, but they are not representations of all three rook relations. The nontrivial widths studied here split into two
sequences:

$$
k=3s+1: 4,7,10,13,16,19,\ldots
$$

and

$$
k=3s+2: 5,8,11,14,17,20,\ldots.
$$

Widths $k=3s$ still define full symmetric rook matrices, but their three rook relations do not have this compact circular
encoding.

## Why the constructed supports are forests

Interpret each selected cell $(r,c)$ as an edge between row vertex $r$ and column vertex $c$ in $K_{3,k}$.

Suppose an ESS support contains a cycle. Assign alternating values $+1,-1,+1,-1,\ldots$ to the cycle edges and zero to the other
support edges. This direction $h$ has

$$
\sum_c h_{rc}=0,\qquad \sum_r h_{rc}=0,
$$

for every row and column, and also $\sum_{r,c}h_{rc}=0$. Hence the row and column terms in the quadratic form vanish, leaving

$$
h^TMh=\delta\sum_{r,c}h_{rc}^2.
$$

If $x$ is the ESS, then $M_Ix=u\mathbf1$ on its support, and therefore $h^TMx=u\mathbf1^Th=0$. For sufficiently small $|t|$,
$x+th$ remains feasible and

$$
(x+th)^TM(x+th)-x^TMx=t^2h^TMh.
$$

For $\delta>0$ the objective rises; for $\delta=0$ it remains unchanged. Either case contradicts strict local maximality.
Therefore:

$$
\boxed{\delta\ge0\quad\Longrightarrow\quad\text{every ESS support is cycle-free}.}
$$

In each construction below, the bridge columns connect all three row vertices and every other column is attached as a leaf. The
support is therefore connected, includes all $k+3$ vertices, and has $k+2$ edges: it is a spanning tree of $K_{3,k}$. In general,
such a spanning tree has

$$
3+k-1=k+2
$$

edges, explaining the observed ESS support size $k+2$ in both sequences.

Because a spanning tree has two more edges than column vertices, its column degrees have only two possibilities:

1. one column has degree three and every other column has degree one; or
2. two columns have degree two and every other column has degree one.

These are the three-way-bridge and two-bridge structures used below.

## Sequence I: $k=3s+1$

Here the game dimension and support size are

$$
N=3k=9s+3,
$$

$$
|I|=k+2=3s+3.
$$

The constructed trees have equal row degrees

$$
\boxed{(s+1,s+1,s+1)}.
$$

Define the equal-part multinomial coefficient

$$
H_s=\binom{3s}{s,s,s}=\frac{(3s)!}{(s!)^3}.
$$

### Type I-a: one three-way bridge

Choose the column incident to all three rows in $k$ ways. The remaining $3s$ columns are leaves, divided into $s$ leaves for each
row. Therefore

$$
E^{(3)}_s=kH_s.
$$

### Type I-b: two two-way bridges

The two bridge columns connect the three rows as a path. Choose the central row in three ways and assign the two distinct bridge
columns to the two endpoint rows in $k(k-1)$ ways. Among the remaining $3s-1$ leaf columns, each endpoint row receives $s$ and
the central row receives $s-1$. Thus

$$
E^{(2,2)}_s
=3k(k-1)\frac{(3s-1)!}{s!\,(s-1)!\,s!}
=k(k-1)H_s.
$$

Adding both support orbits gives

$$
\boxed{
E^{[1]}_s
=k^2H_s
=(3s+1)^2\frac{(3s)!}{(s!)^3}
}.
$$

In terms of the width $k\equiv1\pmod3$,

$$
\boxed{
E^{[1]}(k)
=k^2\frac{(k-1)!}{\left(\left(\frac{k-1}{3}\right)!\right)^3}
}.
$$

### Verified finite sequence

| $s$ | Grid | $N$ | Verified $(A,B,C)$ | Support size | Constructed ESS | $\gamma=E^{1/N}$ |
|---:|---|---:|---|---:|---:|---:|
| 1 | $3\times4$ | 12 | $(4,5,10)$ | 6 | 96 | 1.462814543224 |
| 2 | $3\times7$ | 21 | $(6,7,14)$ | 9 | 4,410 | 1.491230215179 |
| 3 | $3\times10$ | 30 | $(8,10,19)$ | 12 | 168,000 | **1.493402850893** |
| 4 | $3\times13$ | 39 | $(11,11,23)$ | 15 | 5,855,850 | 1.491172710774 |
| 5 | $3\times16$ | 48 | $(13,13,27)$ | 18 | 193,729,536 | 1.488160980927 |
| 6 | $3\times19$ | 57 | $(15,15,31)$ | 21 | 6,192,282,096 | 1.485206546908 |

The $3\times10$ construction is the first displayed sequence member whose $\gamma$ value exceeds the $3\times8$ value reported in
the 2018 paper. The compact circular input for its verified triple is

```text
30#19,19,8,19,19,8,19,19,8,10,19,8,19,19,8
```

The simpler extrapolation $(8,9,18)$ stabilizes the three-way-bridge orbit but not the much larger two-bridge orbit. Using the
tangent basis $Z=(e_1-e_q,\ldots,e_{q-1}-e_q)$ for the chosen $q$-strategy representative, the largest eigenvalue of
$Z^TM_IZ$ is approximately $+0.00809744$. Moving to $(8,10,19)$ keeps $C=A+B+1$ and makes the reduced form exactly negative
definite.

The two orbit counts at $3\times10$ are

$$
10\frac{9!}{3!^3}=16{,}800
$$

and

$$
3\cdot10\cdot9\frac{8!}{3!\,2!\,3!}=151{,}200,
$$

which sum to 168,000.

## Sequence II: $k=3s+2$

Here

$$
N=3k=9s+6,
$$

$$
|I|=k+2=3s+4.
$$

The constructed tree has one distinguished central row. Two degree-two bridge columns connect it to the other rows. Every row has
$s$ additional leaf columns. The row degrees are

$$
\boxed{(s+2,s+1,s+1)}.
$$

Choose the central row in three ways, assign two distinct bridge columns to the endpoint rows in $k(k-1)$ ways, and partition the
remaining $3s$ leaf columns equally among the rows. This gives

$$
\boxed{
E^{[2]}_s
=3k(k-1)H_s
=3(3s+2)(3s+1)\frac{(3s)!}{(s!)^3}
}.
$$

In terms of $k\equiv2\pmod3$,

$$
\boxed{
E^{[2]}(k)
=3k(k-1)\frac{(k-2)!}{\left(\left(\frac{k-2}{3}\right)!\right)^3}
}.
$$

### Verified finite sequence

| $s$ | Grid | $N$ | Verified $(A,B,C)$ | Support size | Constructed ESS | $\gamma=E^{1/N}$ |
|---:|---|---:|---|---:|---:|---:|
| 1 | $3\times5$ | 15 | $(4,5,10)$ | 7 | 360 | 1.480540073044 |
| 2 | $3\times8$ | 24 | $(6,7,14)$ | 10 | 15,120 | **1.493303184701** |
| 3 | $3\times11$ | 33 | $(8,10,19)$ | 13 | 554,400 | 1.492984328761 |
| 4 | $3\times14$ | 42 | $(11,11,23)$ | 16 | 18,918,900 | 1.490251013137 |
| 5 | $3\times17$ | 51 | $(13,13,27)$ | 19 | 617,512,896 | 1.487187095878 |
| 6 | $3\times20$ | 60 | $(15,15,31)$ | 22 | 19,554,575,040 | 1.484296988716 |

The published $3\times8$ matrix uses the different two-value triple $(7,7,15)$. Both $(6,7,14)$ and $(7,7,15)$ produce the same
complete count and the same rook-tree orbit of 15,120 supports, but the matrices are not claimed to be strategically equivalent in
general. The paper states that its dimension-24 search was restricted to a two-parameter family and reports, among the vectors it
found with the maximal count, one minimizing the infinity norm. The later three-variable integer search found the same count at
$(6,7,14)$.

For $3\times11$, the verified compact input is

```text
33#19,19,8,19,19,8,19,19,8,19,10,8,19,19,8,19
```

Its constructed 554,400 ESS are one orbit of trees with row degrees $(5,4,4)$ and column degrees
$(2,2,1,1,1,1,1,1,1,1,1)$.

## Congruence class zero: $k=3s$

The remaining widths have

$$
k=3s,
\qquad
N=3k=9s,
\qquad
|I|=k+2=3s+2.
$$

The full symmetric rook matrix $M(A,B,C)$ still exists and retains the action of $S_3\times S_{3s}$. For distinct payoff
values this is the full automorphism group for every $s$. At $s=1$, transposing the square $3\times3$ grid supplies an
additional symmetry only on the slice $A=B$. What fails throughout this congruence class is the compact circular
representation of all three rook relations:

$$
\gcd(3,3s)=3.
$$

The most balanced spanning trees have row-degree multiset

$$
\boxed{(s+1,s+1,s)}.
$$

As before, their column degrees contain either one degree-three bridge or two degree-two bridges. These possibilities
now split into three natural rook-tree orbits. Recall

$$
H_s=\frac{(3s)!}{(s!)^3}.
$$

### Type 0-a: one three-way bridge

Choose the bridge column in $k$ ways and the row of degree $s$ in three ways. The remaining $3s-1$ leaf columns are
distributed with counts $(s,s,s-1)$. Hence

$$
E^{(3)}_{0,s}
=3k\frac{(3s-1)!}{s!\,s!\,(s-1)!}
=\boxed{kH_s}.
$$

### Type 0-b: two bridges with the central row of degree $s$

Choose the central row and assign the two distinct bridge columns as in Sequence I. The central row then receives $s-2$
leaf columns, while both endpoint rows receive $s$. This orbit is empty for $s=1$ and otherwise has

$$
E^{(2,2)}_{0,s,\mathrm{low}}
=3k(k-1)\frac{(3s-2)!}{s!\,(s-2)!\,s!}
=\boxed{k(s-1)H_s}.
$$

### Type 0-c: two bridges with the central row of degree $s+1$

The central row receives $s-1$ leaves. One endpoint receives $s$ leaves and the other receives $s-1$, with two choices
for which endpoint has the larger degree. Therefore

$$
E^{(2,2)}_{0,s,\mathrm{high}}
=6k(k-1)\frac{(3s-2)!}{s!\,(s-1)!\,(s-1)!}
=\boxed{2ksH_s}.
$$

Adding the three support orbits gives the exact count

$$
\boxed{
E^{[0]}_s
=kH_s+k(s-1)H_s+2ksH_s
=k^2H_s
}.
$$

If all three orbits were ESS for one common payoff triple, this would also be their ESS count. The first support counts are:

| $s$ | Grid | $N$ | Support size | Balanced support count |
|---:|---|---:|---:|---:|
| 1 | $3\times3$ | 9 | 5 | 54 |
| 2 | $3\times6$ | 18 | 8 | 3,240 |
| 3 | $3\times9$ | 27 | 11 | 136,080 |
| 4 | $3\times12$ | 36 | 14 | 4,989,600 |

These are not ESS lower bounds. For the tested cases $s=1,\ldots,7$, the regular continuation

$$
(A_s,B_s,C_s)=(2s+3,2s+3,4s+7)
$$

does not make all existing $k=3s$ orbits ESS. A binary64 search over every normalized integer pair $A,B\in[-50,50]$, with
$C=A+B+1$, for the same values of $s$ found no common solution. In a separate binary64 screen, 20,000 pairs drawn uniformly
from $[-1000,1000]^2$ for each existing support type and each $s=1,2,3$ did not produce even one ESS representative.
Some tested parameters satisfy the candidate conditions and others satisfy negative definiteness, but none tested
satisfies both requirements for every required orbit.

These finite floating-point searches are not an impossibility proof. They indicate that the missing modulo-three case may have a
genuine feasibility obstruction rather than merely lacking a compact circular encoding. A proof would require symbolic candidate
and stability inequalities showing that their feasible regions are disjoint.

## Shared payoff observations

The same verified triple works for each adjacent pair of widths $3s+1$ and $3s+2$:

| $s$ | Widths | Verified triple |
|---:|---|---|
| 1 | 4 and 5 | $(4,5,10)$ |
| 2 | 7 and 8 | $(6,7,14)$ |
| 3 | 10 and 11 | $(8,10,19)$ |
| 4 | 13 and 14 | $(11,11,23)$ |
| 5 | 16 and 17 | $(13,13,27)$ |
| 6 | 19 and 20 | $(15,15,31)$ |

Every row satisfies $C=A+B+1$. A simple regular family

$$
(A_s,B_s,C_s)=(2s+3,2s+3,4s+7)
$$

also satisfies this normalization and numerically stabilizes all three representatives belonging to Sequences I and II
for $1\le s\le15$. For the tested cases $1\le s\le7$, it does not stabilize all existing $k=3s$ orbits described above.
The successful checks are finite evidence for the two proposed infinite ESS families, not an all-$s$ symbolic stability proof.
They were representative-only numerical checks, not complete FracESSA searches. For $s>6$, the game dimensions exceed
FracESSA's supported maximum of 63. The exact finite claims in the tables use the explicitly listed triples and exact arithmetic.
At $s=3$, the two-value choice $(9,9,19)$ also passes exact candidate and stability checks for both Sequence I representatives and
therefore certifies the same constructed 168,000-ESS lower bound as $(8,10,19)$.

## Exact verification method and evidence

For each support-orbit representative and displayed payoff triple, the verification performed the following exact operations:

1. solve the bordered equilibrium system over rational numbers;
2. require every support probability to be strictly positive;
3. require every strategy outside the support to have strictly smaller payoff;
4. form a basis $Z$ of the support tangent space $\{h:\mathbf1^Th=0\}$;
5. verify that $-Z^TM_IZ$ is positive definite by exact Sylvester principal minors.

Because the matrix automorphism group $S_3\times S_k$ maps every support in an orbit to every other support in that orbit, one exact
representative check proves the ESS property for the complete counted orbit.

### Sequence I exact margins

The table lists the minimum support probability and the largest outside payoff difference. Every difference is negative, and both
the three-way and two-bridge reduced matrices passed exact positive-definiteness checks after negation.

| $k$ | Triple | Three-way: $\min x_i$ | Three-way: max outside gap | Two-bridge: $\min x_i$ | Two-bridge: max outside gap |
|---:|---|---:|---:|---:|---:|
| 4 | $(4,5,10)$ | $2/27$ | $-7/27$ | $7/127$ | $-7/127$ |
| 7 | $(6,7,14)$ | $1/23$ | $-10/69$ | $26/815$ | $-26/815$ |
| 10 | $(8,10,19)$ | $4/129$ | $-13/129$ | $188/8205$ | $-188/8205$ |
| 13 | $(11,11,23)$ | $11/453$ | $-35/453$ | $803/43425$ | $-803/43425$ |
| 16 | $(13,13,27)$ | $13/654$ | $-41/654$ | $1287/85298$ | $-1287/85298$ |
| 19 | $(15,15,31)$ | $5/297$ | $-47/891$ | $129/10117$ | $-129/10117$ |

### Sequence II exact margins

| $k$ | Triple | $\min x_i$ | Largest outside gap | Reduced stability |
|---:|---|---:|---:|---|
| 5 | $(4,5,10)$ | $4/79$ | $-9/79$ | negative definite |
| 8 | $(6,7,14)$ | $13/421$ | $-169/2526$ | negative definite |
| 11 | $(8,10,19)$ | $376/16703$ | $-799/16703$ | negative definite |
| 14 | $(11,11,23)$ | $803/44216$ | $-1679/44216$ | negative definite |
| 17 | $(13,13,27)$ | $1287/86375$ | $-2673/86375$ | negative definite |
| 20 | $(15,15,31)$ | $215/17018$ | $-1333/51054$ | negative definite |

These exact checks establish the stated orbit counts as genuine ESS lower bounds without enumerating the complete $2^N$ support
space.

### Complete enumeration boundary

The current exact `safe` search also exhaustively covered all support orbits for the four smallest games. In each case its total
ESS count equals the constructed orbit count:

| Grid | Triple | Constructed ESS | Complete ESS count |
|---|---|---:|---:|
| $3\times4$ | $(4,5,10)$ | 96 | 96 |
| $3\times5$ | $(4,5,10)$ | 360 | 360 |
| $3\times7$ | $(6,7,14)$ | 4,410 | 4,410 |
| $3\times8$ | $(6,7,14)$ | 15,120 | 15,120 |

For $s\ge3$, no complete enumeration is claimed. The larger values in the sequence tables remain rigorous lower bounds obtained
from the exactly verified support orbits.

## Gamma comparison

The best displayed value in Sequence I is

$$
168000^{1/30}\approx1.493402850893.
$$

The published Sequence II value is

$$
15120^{1/24}\approx1.493303184701.
$$

Thus the constructed $3\times10$ lower bound is slightly larger:

$$
1.493402850893>1.493303184701.
$$

At dimension 30, matching the old exponential base requires more than

$$
15120^{30/24}\approx167663.968
$$

ESS. The constructed 168,000 exceed the integer threshold 167,664 by only 336. The improvement is therefore real but very small.

For $3\times11$, matching the old base requires at least 558,321 ESS. Its constructed orbit contains 554,400, short by 3,921.
Additional ESS or another parameter chamber could change that comparison; the orbit formula alone does not.

## Conditional asymptotics

Stirling's formula gives

$$
H_s=\frac{(3s)!}{(s!)^3}
\sim\frac{\sqrt3}{2\pi s}\,27^s.
$$

Therefore the support-count formulas satisfy

$$
E^{[1]}_s\sim\frac{9\sqrt3}{2\pi}s\,27^s
$$

and

$$
E^{[2]}_s\sim\frac{27\sqrt3}{2\pi}s\,27^s.
$$

The congruence-zero balanced-support total has

$$
E^{[0]}_s=k^2H_s
\sim\frac{9\sqrt3}{2\pi}s\,27^s.
$$

Since all three dimensions are asymptotic to $9s$, each count has

$$
\lim_{s\to\infty}\left(E_s\right)^{1/N}
=27^{1/9}
=3^{1/3}
\approx1.442249570307.
$$

For Sequences I and II, interpreting this as an ESS lower bound is conditional on the corresponding support orbits remaining ESS
for all $s$. The regular payoff family above provides numerical evidence through $s=15$, while the exact verification in this
note covers $s=1,\ldots,6$. For $k=3s$, the same limit describes actual balanced-support counts but not established ESS counts,
because no common stabilizing payoff triple has been found.

## What is established and what remains open

Established here:

1. the quadratic-form decomposition and exact $\delta=1$ normalization;
2. the compact-encoding split into $k=3s+1$ and $k=3s+2$, and its obstruction for $k=3s$;
3. exclusion of cyclic ESS supports when $\delta\ge0$;
4. exact combinatorial counts of the displayed spanning-tree support orbits;
5. exact ESS verification of those orbits through $k=20$, hence game dimension 60;
6. complete counts for the four cases already enumerated by FracESSA;
7. a rigorously constructed 168,000-ESS lower bound at dimension 30;
8. the exact count of the three natural balanced $k=3s$ support orbits, together with finite search evidence that does not
   establish them as ESS.

Not established:

1. that the listed games have no additional ESS for $s\ge3$;
2. that either support formula is maximal within the $3\times k$ rook family;
3. that any displayed game maximizes ESS over all symmetric matrices of its dimension;
4. that the regular payoff family remains stable for every $s$;
5. that the $\delta=0$ or $\delta<0$ regions cannot produce larger records;
6. whether any payoff triple makes all three natural $k=3s$ orbits ESS, or whether their candidate and stability regions are
   provably disjoint.

A next mathematical step would be an all-$s$ symbolic candidate and stability proof for
$(A_s,B_s,C_s)=(2s+3,2s+3,4s+7)$. A separate optimization question is whether other forest degree patterns or other normalized
parameter chambers can exceed the verified Sequence I and II counts. For $k=3s$, the first step is instead to decide whether the
three displayed orbit types can be ESS simultaneously at all.
