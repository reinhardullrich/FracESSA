# Largest Rectangular-Rook Automorphism Groups for Dimensions 2 to 30

## Scope

This table continues the rectangular rook construction from
[`AUTOMORPHISM_GROUPS_AND_GAMMA.md`](AUTOMORPHISM_GROUPS_AND_GAMMA.md). It considers grids

$$
n=mk,
\qquad 2\leq m\leq k,
$$

whose strategies are the cells of an $m\times k$ grid. Two distinct cells receive one payoff when they share a row or column and
another payoff otherwise.

For $m\neq k$, the full automorphism group is

$$
S_m\times S_k
$$

and has order $m!k!$. For a square grid, transposition supplies one additional symmetry:

$$
(S_m\times S_m)\rtimes C_2,
$$

of order $2(m!)^2$.

The degenerate factorization $1\times n$ is excluded. It makes every pair of distinct cells share the only row, so the matrix has
the full group $S_n$ and reduces to the already-discussed constant-off-diagonal family. Its large group does not imply a large ESS
count or a large $\gamma$.

This restriction means that three groups do not exist for every dimension. Prime dimensions have no nondegenerate rectangular
factorization, and most composite dimensions have only one or two. The table lists every available group, ranked by order, and
uses `none` rather than mixing unrelated matrix constructions into the ranking.

## Ranking

| $n$ | Largest rook group | Second largest | Third largest |
|---:|---|---|---|
| 2 | none | none | none |
| 3 | none | none | none |
| 4 | $2\times2$: $(S_2\times S_2)\rtimes C_2$ (8) | none | none |
| 5 | none | none | none |
| 6 | $2\times3$: $S_2\times S_3$ (12) | none | none |
| 7 | none | none | none |
| 8 | $2\times4$: $S_2\times S_4$ (48) | none | none |
| 9 | $3\times3$: $(S_3\times S_3)\rtimes C_2$ (72) | none | none |
| 10 | $2\times5$: $S_2\times S_5$ (240) | none | none |
| 11 | none | none | none |
| 12 | $2\times6$: $S_2\times S_6$ (1,440) | $3\times4$: $S_3\times S_4$ (144) | none |
| 13 | none | none | none |
| 14 | $2\times7$: $S_2\times S_7$ (10,080) | none | none |
| 15 | $3\times5$: $S_3\times S_5$ (720) | none | none |
| 16 | $2\times8$: $S_2\times S_8$ (80,640) | $4\times4$: $(S_4\times S_4)\rtimes C_2$ (1,152) | none |
| 17 | none | none | none |
| 18 | $2\times9$: $S_2\times S_9$ (725,760) | $3\times6$: $S_3\times S_6$ (4,320) | none |
| 19 | none | none | none |
| 20 | $2\times10$: $S_2\times S_{10}$ (7,257,600) | $4\times5$: $S_4\times S_5$ (2,880) | none |
| 21 | $3\times7$: $S_3\times S_7$ (30,240) | none | none |
| 22 | $2\times11$: $S_2\times S_{11}$ (79,833,600) | none | none |
| 23 | none | none | none |
| 24 | $2\times12$: $S_2\times S_{12}$ (958,003,200) | $3\times8$: $S_3\times S_8$ (241,920) | $4\times6$: $S_4\times S_6$ (17,280) |
| 25 | $5\times5$: $(S_5\times S_5)\rtimes C_2$ (28,800) | none | none |
| 26 | $2\times13$: $S_2\times S_{13}$ (12,454,041,600) | none | none |
| 27 | $3\times9$: $S_3\times S_9$ (2,177,280) | none | none |
| 28 | $2\times14$: $S_2\times S_{14}$ (174,356,582,400) | $4\times7$: $S_4\times S_7$ (120,960) | none |
| 29 | none | none | none |
| 30 | $2\times15$: $S_2\times S_{15}$ (2,615,348,736,000) | $3\times10$: $S_3\times S_{10}$ (21,772,800) | $5\times6$: $S_5\times S_6$ (86,400) |

## Exact Symbolic Matrix Strings

Every group in the ranking needs only two distinct off-diagonal values:

- `P`: two different grid cells are in the same row or the same column;
- `Q`: the cells are in different rows and different columns.

`P` and `Q` must be different. If they are equal, all off-diagonal entries are equal and the automorphism group grows to $S_n$.
The published dimension-24 matrix uses `P = 7` and `Q = 15`.

A non-square grid also permits the more general three-variable form:

- `A`: same row;
- `B`: same column;
- `C`: different row and different column.

Distinct `A`, `B`, and `C` preserve the exact group $S_m\times S_k$ when $m\neq k$. For a square grid, however, `A` and `B` must
be equal to retain the row-column transposition in $(S_m\times S_m)\rtimes C_2$; a genuinely three-variable form has only the
smaller group $S_m\times S_m$. The common diagonal is zero and absent from the compact strings. A full symmetric matrix could give
the diagonal a fourth variable without changing its automorphism group.

For coprime $m$ and $k$, Chinese remainder indexing makes the rook matrix circular. At compact offset $d$, the three-variable value
is `A` when $m\mid d$, `B` when $k\mid d$, and `C` otherwise. Non-coprime grids require the full symmetric-matrix format. The
$2\times2$ rook graph is the sole square exception for the two-variable compact representation because it is a four-cycle.

These are symbolic templates, not literal parser input: replace their letters with rational values first.

For a support $S$, let

$$
M_G=\max_S|G\cdot S|=\frac{|G|}{\min_S|\operatorname{Stab}_G(S)|}.
$$

The penultimate column gives $M_G$ and $\gamma_G=M_G^{1/n}$. It is the symmetry-only ceiling for one ESS orbit: if a support attaining
the largest orbit is an ESS, all $M_G$ transformed supports are ESS. It is not an absolute maximum when several ESS orbits coexist,
and the listed payoff templates do not guarantee that the ceiling is attainable. Using $|G|$ directly would be wrong whenever every
support has a nontrivial stabilizer. For example, the $3\times8$ action has minimum stabilizer 2, so its ceiling is 120,960 rather
than 241,920; the published support has stabilizer 16 and produces 15,120 ESS.

For every game of dimension $n$, the ESS supports form an inclusion antichain. Sperner's theorem therefore gives the absolute
upper bound

$$
\gamma_{\max}(n)=\binom{n}{\lfloor n/2\rfloor}^{1/n}.
$$

The `Largest ESS count found` column is the maximum stored `ess_count` over all matrices of that dimension in the canonical test
database, regardless of matrix family or circular symmetry. It is an empirical lower bound, not a theorem: it can decrease when
the database has weaker coverage at a larger dimension. The last column records the independent Sperner upper bound, which
approaches 2 from below as $n$ grows.

| $n$ | Rank | Exact rook group | Minimum variables | Two-variable compact string | Three-variable compact string | Maximum one-orbit ESS; $\gamma_G$ | Largest ESS count found | $\gamma_{\max}(n)$ |
|---:|---:|---|---:|---|---|---:|---:|---:|
| 2 | - | none | - | - | - | - | 2 | 1.414214 |
| 3 | - | none | - | - | - | - | 3 | 1.442250 |
| 4 | 1 | $2\times2$: $(S_2\times S_2)\rtimes C_2$ (8) | 2 | `4#P,Q` | not available for the listed full group | 4; 1.414214 | 4 | 1.565085 |
| 5 | - | none | - | - | - | - | 6 | 1.584893 |
| 6 | 1 | $2\times3$: $S_2\times S_3$ (12) | 2 | `6#Q,P,P` | `6#C,A,B` | 12; 1.513086 | 9 | 1.647549 |
| 7 | - | none | - | - | - | - | 14 | 1.661809 |
| 8 | 1 | $2\times4$: $S_2\times S_4$ (48) | 2 | not possible in compact circular format | not possible in compact circular format | 24; 1.487738 | 20 | 1.700737 |
| 9 | 1 | $3\times3$: $(S_3\times S_3)\rtimes C_2$ (72) | 2 | not possible in compact circular format | not available for the listed full group | 36; 1.489095 | 30 | 1.711491 |
| 10 | 1 | $2\times5$: $S_2\times S_5$ (240) | 2 | `10#Q,P,Q,P,P` | `10#C,A,C,A,B` | 120; 1.614054 | 50 | 1.738361 |
| 11 | - | none | - | - | - | - | 70 | 1.746788 |
| 12 | 1 | $2\times6$: $S_2\times S_6$ (1,440) | 2 | not possible in compact circular format | not possible in compact circular format | 360; 1.633147 | 105 | 1.766604 |
| 12 | 2 | $3\times4$: $S_3\times S_4$ (144) | 2 | `12#Q,Q,P,P,Q,P` | `12#C,C,A,B,C,A` | 144; 1.513086 | 105 | 1.766604 |
| 13 | - | none | - | - | - | - | 157 | 1.773409 |
| 14 | 1 | $2\times7$: $S_2\times S_7$ (10,080) | 2 | `14#Q,P,Q,P,Q,P,P` | `14#C,A,C,A,C,A,B` | 1,260; 1.665156 | 233 | 1.788707 |
| 15 | 1 | $3\times5$: $S_3\times S_5$ (720) | 2 | `15#Q,Q,P,Q,P,P,Q` | `15#C,C,A,C,B,A,C` | 720; 1.550561 | 360 | 1.794334 |
| 16 | 1 | $2\times8$: $S_2\times S_8$ (80,640) | 2 | not possible in compact circular format | not possible in compact circular format | 3,360; 1.661102 | 536 | 1.806544 |
| 16 | 2 | $4\times4$: $(S_4\times S_4)\rtimes C_2$ (1,152) | 2 | not possible in compact circular format | not available for the listed full group | 1,152; 1.553606 | 536 | 1.806544 |
| 17 | - | none | - | - | - | - | 784 | 1.811287 |
| 18 | 1 | $2\times9$: $S_2\times S_9$ (725,760) | 2 | `18#Q,P,Q,P,Q,P,Q,P,P` | `18#C,A,C,A,C,A,C,A,B` | 15,120; 1.706858 | 1,164 | 1.821288 |
| 18 | 2 | $3\times6$: $S_3\times S_6$ (4,320) | 2 | not possible in compact circular format | not possible in compact circular format | 4,320; 1.592104 | 1,164 | 1.821288 |
| 19 | - | none | - | - | - | - | 1,694 | 1.825348 |
| 20 | 1 | $2\times10$: $S_2\times S_{10}$ (7,257,600) | 2 | not possible in compact circular format | not possible in compact circular format | 50,400; 1.718389 | 2,560 | 1.833707 |
| 20 | 2 | $4\times5$: $S_4\times S_5$ (2,880) | 2 | `20#Q,Q,Q,P,P,Q,Q,P,Q,P` | `20#C,C,C,A,B,C,C,A,C,B` | 2,880; 1.489257 | 2,560 | 1.833707 |
| 21 | 1 | $3\times7$: $S_3\times S_7$ (30,240) | 2 | `21#Q,Q,P,Q,Q,P,P,Q,P,Q` | `21#C,C,A,C,C,A,B,C,A,C` | 15,120; 1.581344 | 4,410 | 1.837228 |
| 22 | 1 | $2\times11$: $S_2\times S_{11}$ (79,833,600) | 2 | `22#Q,P,Q,P,Q,P,Q,P,Q,P,P` | `22#C,A,C,A,C,A,C,A,C,A,B` | 184,800; 1.735384 | 6,300 | 1.844331 |
| 23 | - | none | - | - | - | - | 9,060 | 1.847419 |
| 24 | 1 | $2\times12$: $S_2\times S_{12}$ (958,003,200) | 2 | not possible in compact circular format | not possible in compact circular format | 554,400; 1.735106 | 15,120 | 1.853537 |
| 24 | 2 | $3\times8$: $S_3\times S_8$ (241,920) | 2 | `24#Q,Q,P,Q,Q,P,Q,P,P,Q,Q,P` | `24#C,C,A,C,C,A,C,B,A,C,C,A` | 120,960; 1.628459 | 15,120 | 1.853537 |
| 24 | 3 | $4\times6$: $S_4\times S_6$ (17,280) | 2 | not possible in compact circular format | not possible in compact circular format | 17,280; 1.501635 | 15,120 | 1.853537 |
| 25 | 1 | $5\times5$: $(S_5\times S_5)\rtimes C_2$ (28,800) | 2 | not possible in compact circular format | not available for the listed full group | 28,800; 1.507911 | 8,748 | 1.856270 |
| 26 | 1 | $2\times13$: $S_2\times S_{13}$ (12,454,041,600) | 2 | `26#Q,P,Q,P,Q,P,Q,P,Q,P,Q,P,P` | `26#C,A,C,A,C,A,C,A,C,A,C,A,B` | 2,402,400; 1.759582 | 13 | 1.861602 |
| 27 | 1 | $3\times9$: $S_3\times S_9$ (2,177,280) | 2 | not possible in compact circular format | not possible in compact circular format | 544,320; 1.630944 | 27 | 1.864040 |
| 28 | 1 | $2\times14$: $S_2\times S_{14}$ (174,356,582,400) | 2 | not possible in compact circular format | not possible in compact circular format | 8,408,400; 1.767304 | 105 | 1.868733 |
| 28 | 2 | $4\times7$: $S_4\times S_7$ (120,960) | 2 | `28#Q,Q,Q,P,Q,Q,P,P,Q,Q,Q,P,Q,P` | `28#C,C,C,A,C,C,B,A,C,C,C,A,C,B` | 120,960; 1.518878 | 105 | 1.868733 |
| 29 | - | none | - | - | - | - | 29 | 1.870924 |
| 30 | 1 | $2\times15$: $S_2\times S_{15}$ (2,615,348,736,000) | 2 | `30#Q,P,Q,P,Q,P,Q,P,Q,P,Q,P,Q,P,P` | `30#C,A,C,A,C,A,C,A,C,A,C,A,C,A,B` | 31,531,500; 1.778108 | 30 | 1.875090 |
| 30 | 2 | $3\times10$: $S_3\times S_{10}$ (21,772,800) | 2 | `30#Q,Q,P,Q,Q,P,Q,Q,P,P,Q,P,Q,Q,P` | `30#C,C,A,C,C,A,C,C,A,B,C,A,C,C,A` | 5,443,200; 1.676982 | 30 | 1.875090 |
| 30 | 3 | $5\times6$: $S_5\times S_6$ (86,400) | 2 | `30#Q,Q,Q,Q,P,P,Q,Q,Q,P,Q,P,Q,Q,P` | `30#C,C,C,C,A,B,C,C,C,A,C,B,C,C,A` | 86,400; 1.460664 | 30 | 1.875090 |

## Important Interpretation

The order of the matrix automorphism group measures how many strategy relabelings preserve the matrix. It does not determine how
many ESS the matrix has. In particular:

- the largest rook group at a given dimension need not produce the largest $\gamma$;
- the published dimension-24 record uses the second-ranked $3\times8$ grid, not the larger $2\times12$ group;
- each support has its own stabilizer, so its orbit multiplier is the group order divided by a support-specific denominator;
- a matrix's numerical payoff values decide which support orbits are candidates and which are ESS.

For search purposes, this table supplies symmetry families to test. It is not a ranking of expected ESS counts.

## Boundary of This Result

This is not a classification of the three largest automorphism groups among all possible symmetric $n\times n$ matrices. Without
the rook-family restriction, the full group $S_n$ is always available, and identifying the next largest exact automorphism groups
requires specifying which matrix structures are admissible. Such a classification should not be mixed with this factorization
table.
