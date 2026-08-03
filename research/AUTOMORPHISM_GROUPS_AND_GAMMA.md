# Automorphism Groups and High-Gamma ESS Games

This note records the symmetry analysis of the high-ESS matrices discussed in the Bomze-Schachinger-Ullrich papers. The practical
objective is to maximize

$$
\gamma(A)=|\operatorname{ESS}(A)|^{1/n}
$$

for symmetric games of dimension $n\leq 30$. A large automorphism group is useful because it identifies equivalent supports and
can reduce the search, but group size alone is not the objective: a highly symmetric matrix can still have very few ESS.

The matrix, ESS count, support size, and $\gamma$ value below come from
[`bomze_schachinger_ullrich_2018_complexity_of_simple_models.md`](papers/bomze_schachinger_ullrich_2018_complexity_of_simple_models.md).
The rook interpretation, complete automorphism group, quadratic-form decomposition, and spectrum are algebraic derivations from
that matrix. The orbit statement is a separate computation against the project database. The other dimensions are untested
constructions, not results from the paper.

## Published dimension-24 record

The 2018 paper gives the compact circular-symmetric matrix

~~~text
24#15,15,7,15,15,7,15,7,7,15,15,7
~~~

Its complete last-column vector is

$$
(15,15,7,15,15,7,15,7,7,15,15,7,15,15,7,7,15,7,15,15,7,15,15,0).
$$

The paper reports:

- dimension: $24$,
- ESS count: $15{,}120$,
- support size of every ESS: $10$,
- growth value: $15{,}120^{1/24}\approx 1.4933$.

This is the highest published $\gamma$ value considered in this note.

## Offset description

Index rows and columns by $\mathbb Z_{24}$. The circulant matrix is determined by

$$
a_0=0,
\qquad
a_d=
\begin{cases}
7, & d\neq 0\text{ and }(3\mid d\text{ or }8\mid d),\\
15, & \text{otherwise}.
\end{cases}
$$

Thus the value $7$ occurs at the nonzero offsets

$$
\{3,6,8,9,12,15,16,18,21\},
$$

and the value $15$ occurs at

$$
\{1,2,4,5,7,10,11,13,14,17,19,20,22,23\}.
$$

## The hidden 3 by 8 rook structure

The Chinese remainder theorem gives

$$
\mathbb Z_{24}\cong \mathbb Z_3\times\mathbb Z_8.
$$

After this relabeling, each strategy is a cell $(r,c)$ of a $3\times 8$ grid. Two distinct strategies receive payoff

- $7$ when their cells are in the same row or the same column,
- $15$ when they are in different rows and different columns,
- $0$ on the diagonal.

The matrix is therefore a weighted rook graph. In grid order it can be written as

$$
A=15(J_{24}-I_{24})
-8\left[I_3\otimes(J_8-I_8)+(J_3-I_3)\otimes I_8\right],
$$

or equivalently

$$
A=15J_{24}+I_{24}-8(I_3\otimes J_8+J_3\otimes I_8).
$$

## Automorphism group

The ordinary circular representation exposes only the dihedral group $D_{24}$, of order $48$: 24 rotations and 24
reflections. Multiplication of indices by any unit modulo 24 also preserves the low-offset set:

$$
U(24)=\{1,5,7,11,13,17,19,23\}.
$$

Together with translations, this gives the affine subgroup

$$
\mathbb Z_{24}\rtimes U(24)
$$

of order $24\cdot 8=192$.

The complete automorphism group is larger. Any permutation of the three rows and any independent permutation of the eight columns
preserves the matrix, so

$$
\operatorname{Aut}(A)\cong S_3\times S_8,
\qquad
|\operatorname{Aut}(A)|=3!\,8!=241{,}920.
$$

This is the full group: the maximal payoff-7 cliques have two distinguishable sizes, namely the three rows of size 8 and the eight
columns of size 3. An automorphism must therefore permute rows among rows and columns among columns. It cannot exchange the two
families.

Compared with the dihedral symmetry currently visible to a bracelet generator, the complete group is larger by the factor

$$
\frac{241{,}920}{48}=5{,}040.
$$

This is a potential reduction in equivalent support work, not a guaranteed end-to-end speedup.

## ESS orbit verification

The database stores this game as matrix 34. Its matrix row records 15,120 candidates and 15,120 ESS. The candidate table contains
345 rotation/reflection representatives whose multipliers sum to 15,120; every representative is an ESS. As an independent
computational audit, each support was converted into a binary $3\times8$ incidence matrix and canonicalized under all row and
column permutations.

All ESS reduce to one canonical row/column-permutation class. One canonical description is the multiset of column masks

$$
(1,1,2,2,3,4,4,5).
$$

The stabilizer of such a support has order $16$. Its orbit under $S_3\times S_8$ consequently has size

$$
\frac{3!\,8!}{16}=15{,}120.
$$

This equals the complete ESS count. Hence all 15,120 ESS of the published game form one orbit under the full automorphism group.
The 345 stored representatives are only representatives under the smaller rotation/reflection action.

## Quadratic form and spectrum

Write a mixed strategy as a $3\times8$ array $x_{rc}$, with row sums $R_r$, column sums $C_c$, and
$\sum_{r,c}x_{rc}=1$. Direct substitution gives

$$
x^TAx
=15+\sum_{r,c}x_{rc}^2
-8\left(\sum_rR_r^2+\sum_cC_c^2\right).
$$

The eigenvalues of the full $24\times24$ matrix are

| Eigenvalue | Multiplicity |
|---:|---:|
| $273$ | 1 |
| $1$ | 14 |
| $-23$ | 7 |
| $-63$ | 2 |

These values follow directly from the tensor-product decomposition into the constant, row-contrast, column-contrast, and
row-and-column-contrast subspaces.

## General rectangular rook family

The same construction works on an $m\times k$ grid, where $n=mk$. Choose diagonal payoff zero, payoff $p$ for two distinct
cells in the same row or column, and payoff $q$ otherwise.

For $m\neq k$, the complete automorphism group is

$$
S_m\times S_k.
$$

For $m=k$, transposition can also exchange rows and columns, giving

$$
(S_m\times S_m)\rtimes C_2.
$$

When $m$ and $k$ are coprime, the grid can be encoded as a circular matrix on $\mathbb Z_{mk}$: for a nonzero offset $d$,
use $p$ when $m\mid d$ or $k\mid d$, and use $q$ otherwise. If $m$ and $k$ are not coprime, this compact circular
encoding does not identify every grid cell uniquely.

The published values $p=7$ and $q=15$ are inherited from the dimension-24 record. They are not known to be optimal for other
grid dimensions. A serious search must vary the payoff ratio as well as $(m,k)$.

## Candidate inputs through dimension 30

The following compact matrices instantiate the same $p=7,q=15$ rook construction. Except for the published $3\times8$ game,
they are proposed search inputs only: their candidate counts, ESS counts, and $\gamma$ values have not been established.

| Grid | $n$ | Automorphism group | Group order | Compact input |
|---|---:|---|---:|---|
| $2\times9$ | 18 | $S_2\times S_9$ | 725,760 | `18#15,7,15,7,15,7,15,7,7` |
| $2\times11$ | 22 | $S_2\times S_{11}$ | 79,833,600 | `22#15,7,15,7,15,7,15,7,15,7,7` |
| $3\times8$ | 24 | $S_3\times S_8$ | 241,920 | `24#15,15,7,15,15,7,15,7,7,15,15,7` |
| $4\times7$ | 28 | $S_4\times S_7$ | 120,960 | `28#15,15,15,7,15,15,7,7,15,15,15,7,15,7` |
| $2\times15$ | 30 | $S_2\times S_{15}$ | 2,615,348,736,000 | `30#15,7,15,7,15,7,15,7,15,7,15,7,15,7,7` |
| $3\times10$ | 30 | $S_3\times S_{10}$ | 21,772,800 | `30#15,15,7,15,15,7,15,15,7,7,15,7,15,15,7` |
| $5\times6$ | 30 | $S_5\times S_6$ | 86,400 | `30#15,15,15,15,7,7,15,15,15,7,15,7,15,15,7` |

The $2\times10$ construction at $n=20$ has group $S_2\times S_{10}$, of order $7{,}257{,}600$, but
$\gcd(2,10)\neq1$. It therefore requires a full symmetric grid matrix rather than the compact circular encoding.

The group orders in this table indicate symmetry available to an algorithm. They say nothing by themselves about the ESS count.

## Why the largest group is not the answer

At dimension 30, a matrix of the form

$$
A=\alpha I+\beta J
$$

has the complete symmetric group $S_{30}$, of order $30!$. For example, the compact circular input for $I-J$ is

~~~text
30#-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
~~~

Yet this game has only 30 pure ESS, so

$$
\gamma=30^{1/30}\approx1.120,
$$

far below the published $1.4933$. Maximizing $|\operatorname{Aut}(A)|$ is therefore not equivalent to maximizing $\gamma$.

## Search direction

For every grid and payoff choice:

1. Compute the candidates and ESS with exact arithmetic.
2. Canonicalize ESS supports under the full grid automorphism group, not only rotations and reflections.
3. Record orbit representatives and stabilizer sizes so the total count can be reconstructed exactly.
4. Compute $\gamma=|\operatorname{ESS}|^{1/n}$.
5. Retain a construction as a new record only when its exact $\gamma$ exceeds $1.4933$.

The main unresolved algorithmic opportunity is to use the full matrix automorphism group during support generation and pruning.
The mathematics above shows that substantial symmetry is currently hidden by the circular representation. It does not yet provide
an implementation or a measured speedup.
