# Detecting All Evolutionarily Stable Strategies

**Author:** I. M. Bomze[^1]<br>
**Journal:** *Journal of Optimization Theory and Applications*, Vol. 75, No. 2, November 1992, pp. 313-329<br>
**Communicated by:** G. Leitmann<br>
**Source PDF:** [bomze_1992.pdf](bomze_1992.pdf)

> Complete text transcription of the paper. Mathematical notation has been reconstructed as LaTeX from the printed pages. Journal running headers, page numbers, and repeated publication footers are omitted.

## Abstract

In evolutionary game theory, the central solution concept is the evolutionarily stable state, which also can be interpreted as an evolutionarily stable population strategy (ESS). As such, this notion is a refinement of the Nash equilibrium concept in that it requires an additional stability property. In the present paper, an algorithm for detecting all ESSs of a given evolutionary game consisting of pairwise conflicts is presented which both is efficient and complete, since it involves a procedure avoiding the search for unstable equilibria to a considerable extent, and also has a finite, exact routine to check evolutionary stability of a given equilibrium. The article also contains the generalization of these results to the playing-the-field setting, where the payoff is nonlinear.

**Key Words:** Copositivity, efficient search, evolutionary games, extremal equilibria, stable equilibria, playing-the-field models.

## 1. Introduction

Consider an animal population consisting of individuals belonging to the same species. We shall assume that, with respect to a certain type of conflict, the behavior of such an individual does not change within its lifetime and is inherited by its offspring. Let us postulate that there are $n$ distinct behavior patterns which we view as pure strategies represented by the standard basis vectors $\{e_1,\ldots,e_n\}$, i.e., the vertices of the standard simplex

$$
S^N
=
\left\{
x\in\mathbb{R}^N:
x_i\geq 0,\ \text{for all }i,\ \text{and }\sum_{i=1}^{n}x_i=1
\right\}
\tag{1}
$$

in $n$-dimensional Euclidean space $\mathbb{R}^N$. Here, $N=\{1,\ldots,n\}$; in order to be compatible with the notation used later in this paper, where not only the number of coordinates, but also their ordering is important, we use the notation $\mathbb{R}^N$ instead of $\mathbb{R}^n$.

The state of a population concerning this conflict is described by $x\in S^N$, where $x_i$ is the frequency of $i$-strategists (i.e., individuals with behavior pattern $e_i$), $1\leq i\leq n$. If an $i$-strategist contends with a $j$-strategist in a pairwise conflict, the payoff to the former will be denoted by $a_{ij}$. This payoff is in terms of incremental Darwinian individual fitness, since we also assume that there are no genetic relations between the individuals. Thus, the payoff structure is fully determined by the payoff or fitness matrix $A=[a_{ij}]_{1\leq i,j\leq n}$. Moreover, we assume that the population is so large that it can ideally be regarded as infinite, and we postulate random mixing: in state $x$, an $i$-strategist will encounter a $j$-strategist with probability $x_j$ (and then receive $a_{ij}$), so that the mean payoff to an $i$-strategist will be

$$
\sum_{j=1}^{n}a_{ij}x_j=(Ax)_i.
\tag{2}
$$

The average mean payoff within a population in state $x$ then will be

$$
\sum_{i=1}^{n}x_i(Ax)_i=x^T Ax,
\tag{3}
$$

$x^T$ always denoting the transpose of $x$.

The prototypical examples in the seminal papers of Maynard Smith and Price (Refs. 1 and 2) dealt with stability of population states $x$ with the aim to explain polymorphism of behavior (a state where several individuals behave differently; cf. also Ref. 3, a forerunner of the above cited articles). However, in recent literature on evolutionary game theory, the monomorphistic approach plays a prominent role, according to which $x$ is interpreted as the average strategy within the population, which, after some time of adaptation, will be played by everyone in the population (as a mixed strategy, say). So a monomorphistic interpretation focuses on evolutionarily stable strategies instead of evolutionarily stable population states. Although the concept of evolutionary stability seems, at least from the point of view of frequency-dependent selection, to be more stringent in the context of polymorphic states than in models under a monomorphistic interpretation, recent developments seem to support the latter to some extent (see Refs. 4 and 5).

Hence, we in the sequel shall stick to the notion of an evolutionarily stable strategy (ESS) as defined in Ref. 1, rendering the evolutionary game as described above as a special bimatrix game $\Gamma(A,A^T)$ with symmetry in payoff and strategies, where the state $x$ is viewed as a mixed strategy over the pure strategies given by the behavior patterns $e_1,\ldots,e_n$; cf. Ref. 6. According to Ref. 1, a strategy $p\in S^N$ is called ESS, if (i) and (ii) below hold:

1. $p$ is a best reply on itself, i.e.,

$$
   x^T Ap\leq p^T Ap,\qquad\text{for all }x\in S^N,
   \tag{4}
$$

   or equivalently $(p,p)$ is a (symmetric) Nash equilibrium point in $\Gamma(A,A^T)$;

2. $p$ fares better against any other reply $x$ on $p$, than $x$ does against itself:

$$
   \text{if }x^T Ap=p^T Ap\text{ and }x\neq p,\qquad\text{then }p^T Ax>x^T Ax.
   \tag{5}
$$

For obvious reasons we refer to conditions (4) and (5) as equilibrium and stability conditions, respectively. The relationship of evolutionary stability to other refinements of the Nash equilibrium concept is discussed in Ref. 6 and in Chapter 9 of Ref. 7.

## 2. Algorithm

There were several attempts to devise a procedure to find ESSs of a given game (Refs. 8 to 10). As noted already in Ref. 11, the algorithms in Refs. 8 and 9 are not able to find all ESSs; also, these articles do not deal with the question of efficient search for ESS candidates among equilibrium strategies. This problem has been attacked, but not been solved explicitly, in Ref. 10 with the help of Proposition 2.1 below. Williams reduces in Ref. 10 the check for stability of a given equilibrium strategy to a linear-quadratic optimization problem, without specifying an exact, finite algorithm to solve it for any given matrix $A$. Recall that the most common procedures do not guarantee convergence even to a local optimum for any matrix $A$, let alone the problem of getting a global one, which is stressed also in Ref. 10.

Here, a procedure is presented which both searches in an explicitly prescribed, efficient way for serious ESS candidates among the equilibria, and is exact in the sense that it specifies a finite routine for checking the stability condition for these candidates, applicable to all possible cases, i.e., all payoff matrices $A$ and all of the candidates generated by the above-mentioned search. For the algorithm described below, we need at first the following support criterion. To this end, we define the support of an arbitrary mixed strategy $x\in S^N$ to be

$$
I(x)=\{i\in N:x_i>0\};
\tag{6}
$$

and we define the extended support of a strategy $p\in S^N$ satisfying condition (4) to be

$$
J(p)=\{j\in N:(Ap)_j=p^T Ap\}.
\tag{7}
$$

Since any equilibrium strategy $p$ satisfies

$$
(Ap)_j\leq p^T Ap,\qquad\text{for all }j\in N,
$$

it always fulfills $I(p)\subseteq J(p)$.

### Proposition 2.1

If $x\in S^N$ is an equilibrium strategy and if $p\in S^N$ is an ESS with $I(x)\subseteq J(p)$, then $x=p$.

**Proof.** See Ref. 6, Theorem 13(iii). $\square$

As a consequence, the set

$$
\mathcal{J}^*=\{I(p):p\text{ is an ESS for }A\}
\tag{8}
$$

forms an antichain in the complete lattice $2^N$ of all subsets of $N$ with respect to set inclusion $\subseteq$. Hence, the number

$$
\binom{n}{\lfloor n/2\rfloor}
$$

is an upper bound for the number of ESSs in a game (see Section 6.3 in Ref. 12). Moreover, the antichain property enables us to devise an efficient routine for traversing the lattice $2^N$ in search of supports $I\in\mathcal{J}^*$ of ESSs. Note that there is at most one ESS $p$ with $I(p)=I$ due to Proposition 2.1, so that the system $\mathcal{J}^*$ completely determines the ESS structure in a given game. Apart from this search technique, the algorithm described below relies on two basic routines, FINDEQ and CHECKSTAB, which will be treated in detail in the subsequent section.

**FINDEQ.** Given a nonvoid subset $I\subseteq N$, this procedure generates an equilibrium strategy $x$ with $I(x)=I$, which is a serious candidate for an ESS (note that there could be an infinity of equilibrium strategies, but only finitely many ESSs), and delivers a negative message if no such candidate exists.

**CHECKSTAB.** This procedure determines whether or not a given equilibrium strategy $p$ also satisfies the stability condition.

### Initialization Step

Check for pure equilibrium strategies $e_i$ (i.e., check the relation $a_{ii}=\max_{j\in N}a_{ji}$), and let

$$
N'=N\setminus\{i\in N:e_i\text{ is an equilibrium strategy}\}.
\tag{9}
$$

Then, any nonpure ESS $p$ satisfies $I(p)\subseteq J(p)\subseteq N'$ due to Proposition 2.1. Hence, for all $i\notin N'$, call CHECKSTAB$(e_i)$, record $e_i$ if the answer is positive, put

$$
\mathcal{J}=\{I\subseteq N':I\text{ has more than one element}\},
\tag{10}
$$

and proceed to the main step.

### Main Step

Denote by $\mathcal{J}_{\max}$ the system of all $I\in\mathcal{J}$ that are maximal with respect to set inclusion in $\mathcal{J}$. For any $I\in\mathcal{J}_{\max}$, call FINDEQ$(I)$; if an equilibrium strategy $p$ is generated, call CHECKSTAB$(p)$; if the answer is positive, record $p$. If there are not any ESSs with support belonging to $\mathcal{J}_{\max}$, then put

$$
\mathcal{J}'=\mathcal{J}\setminus\mathcal{J}_{\max};
\tag{11}
$$

else, denote the ESSs by $p_1,\ldots,p_s$, determine $J(p_i)$, $1\leq i\leq s$, and put

$$
\mathcal{J}'
=
\{I\in\mathcal{J}:I\nsubseteq J(p_i),\ 1\leq i\leq s\}
\setminus\mathcal{J}_{\max}.
\tag{12}
$$

If $\mathcal{J}'$ is void, stop: the list of ESSs recorded up to now is complete. Otherwise, denote by $\mathcal{J}_{\min}$ the system of all $I\in\mathcal{J}'$ which are minimal with respect to set inclusion in $\mathcal{J}'$. For any $I\in\mathcal{J}_{\min}$, call FINDEQ$(I)$. If there are not any equilibrium strategies with support belonging to $\mathcal{J}_{\min}$, then put

$$
\mathcal{J}''=\mathcal{J}'\setminus\mathcal{J}_{\min};
\tag{13}
$$

else, let $x_i$, $1\leq i\leq r$, be all of them; record the ESSs among them by calling CHECKSTAB$(x_i)$; and put

$$
\mathcal{J}''
=
\{I\in\mathcal{J}':I(x_i)\nsubseteq I,\ 1\leq i\leq r\}
\setminus\mathcal{J}_{\min}.
\tag{14}
$$

If $\mathcal{J}''$ is void, stop: the list of ESSs recorded up to now is complete. Otherwise, repeat the main step, replacing $\mathcal{J}$ by $\mathcal{J}''$.

### Example 2.1

Let $n=5$ and

$$
A=
\begin{bmatrix}
1&0&2&2&2\\
0&1&2&2&2\\
2&2&1&0&0\\
2&2&0&1&0\\
2&2&0&0&1
\end{bmatrix}.
$$

Then, according to (9), $N'=N$, because

$$
a_{ii}<\max_{j\in N}a_{ji},\qquad\text{for all }i\in N.
$$

Thus, (10) yields

$$
\mathcal{J}=\{I\subseteq N:I\text{ has at least two elements}\},
$$

with $\mathcal{J}_{\max}=\{N\}$. The only candidate for a completely mixed equilibrium strategy is

$$
p=[5/19,5/19,3/19,3/19,3/19]^T;
$$

indeed, due to the symmetries of $A$, such a candidate has to satisfy

$$
p_1=p_2=a
\qquad\text{and}\qquad
p_3=p_4=p_5=b,
$$

with

$$
2a+3b=1
\qquad\text{and}\qquad
a+6b=4a+b,
$$

due to the condition

$$
N=I(p)\subseteq J(p)\subseteq N.
$$

But since

$$
p^T Ax=\frac{23}{19}<\frac32=x^T Ax,
\qquad
\text{for }x=[1/2,0,0,0,1/2]^T,
$$

$p$ violates the stability condition; hence, due to (11),

$$
\mathcal{J}'
=
\mathcal{J}\setminus\mathcal{J}_{\max}
=
\{I\subseteq N:I\text{ has at least two and at most four elements}\}.
$$

Now,

$$
\mathcal{J}_{\min}
=
\{I\subseteq N:I\text{ has exactly two elements}\}.
$$

Any candidate

$$
q=[t,1-t,0,0,0]^T
$$

with $I(q)=\{1,2\}$, i.e., $0<t<1$, has to fulfill $t=1/2$, again due to the symmetries in $A$. But since

$$
e_3^T A q=2>q^T A q=\frac12,
$$

this candidate $q$ violates the equilibrium condition; similarly, there are no equilibrium strategies $x$ with

$$
I(x)\in\{\{3,4\},\{3,5\},\{4,5\}\}.
$$

Now, the remaining six supports $I_1,\ldots,I_6$ in $\mathcal{J}_{\min}$ are of the form $\{i,k\}$ with $1\leq i\leq2$ and $3\leq k\leq5$. If

$$
r=[t,0,1-t,0,0]^T\in S^N
$$

satisfies

$$
I(r)=I_1=\{1,3\},
$$

then

$$
Ar=[2-t,2(1-t),t+1,2t,2t]^T,
$$

so that $r$ is an equilibrium strategy if and only if $t=1/2$, i.e.,

$$
r=r_1=(1/2)e_1+(1/2)e_3;
$$

furthermore, any other strategy $x\in S^N$ with $I(x)\nsubseteq I_1$ cannot be an alternative best reply to $r_1$, since $J(r_1)=I_1$. Finally, any

$$
x=[u,0,1-u,0,0]^T,
\qquad
\text{with }0\leq u\leq1\text{ and }u\neq1/2,
$$

satisfies

$$
r_1^T Ax=(1/2)(2-u)+(1/2)(u+1)=3/2,
$$

$$
x^T Ax=u(2-u)+(1-u)(u+1)=2u(1-u)+1<3/2,
$$

so that $r_1$ also fulfills the stability condition. By symmetry, we conclude that all $r_j$, $1\leq j\leq6$, of the form $(1/2)e_i+(1/2)e_k$, where $1\leq i\leq2$ and $3\leq k\leq5$, are ESSs with

$$
I(r_j)=J(r_j)=I_j,\qquad 1\leq j\leq6.
$$

Hence, $\{I_j:1\leq j\leq6\}\subseteq\mathcal{J}^*$, and from (14) we obtain

$$
\begin{aligned}
\mathcal{J}''
&=
\{I\subseteq N\text{ with more than two elements}:I_j\nsubseteq I
   \text{ for all }j,\ 1\leq j\leq6\}\\
&=\{\{3,4,5\}\}.
\end{aligned}
$$

Since $\mathcal{J}''\neq\varnothing$, we have to iterate the main step with $\mathcal{J}''$ instead of $\mathcal{J}$, yielding

$$
\mathcal{J}_{\min}=\{\{3,4,5\}\}.
$$

Now, the only candidate for an equilibrium strategy $s$ with

$$
I(s)=\{3,4,5\}
$$

is, again due to the symmetries in $A$, of the form

$$
s=[0,0,1/3,1/3,1/3]^T,
$$

with

$$
As=[2,2,1/3,1/3,1/3]^T,
$$

apparently violating the equilibrium condition. Thus, in the second iteration, $\mathcal{J}'=\varnothing$, and the algorithm stops with $\mathcal{J}^*=\{I_1,\ldots,I_6\}$ corresponding to the six ESSs $r_j$, $1\leq j\leq6$.

Note that the system of all supports has $2^5-1=31$ elements, but that our systematic search procedure involves only five checks for the initialization step (which in this case yielded unfortunately no further reduction: any pure equilibrium strategy found reduces the number of necessary checks by a factor of two). The main procedure calls FINDEQ twelve times and CHECKSTAB only seven times, although there is a total of 21 equilibrium strategies.

[^1]: Associate Professor, Department for Statistics and Computer Science, University of Vienna, Vienna, Austria.

## 3. Selecting Extremal Equilibria and Checking for Evolutionary Stability

To make precise the somewhat fuzzy notion of a serious candidate for an ESS, we use the notion of extremal equilibria in the sense of Vorob'ev and Kuhn, which is helpful in efficiently enumerating those strategies and serves similar purposes here. But instead of following the original approach as described in Ref. 13, pp. 175ff, we use a characterization with the help of a polyhedron introduced in Ref. 10 (but without investigating the relation to extremal equilibria): for any nonvoid $I\subseteq N$, let

$$
\mathcal{P}_I
=
\left\{
\begin{bmatrix}x\\v\end{bmatrix}
\in S^N\times\mathbb{R}:
(Ax)_i\leq v,\ \text{for all }i\in N;\
x_i=0,\ \text{for all }i\notin I
\right\}.
\tag{15}
$$

### Theorem 3.1

Let $p\in S^N$ with $I(p)=I$.

1. If $p$ is an ESS, then $p$ is an extremal equilibrium strategy, i.e., $p$ is an extremal point to the set of all strategies $q$ yielding a Nash equilibrium point $(q,p)$ in the game $\Gamma(A,A^T)$:

$$
   p\in\operatorname{ex}
   \left\{
   q\in S^N:
   x^T Ap\leq q^T Ap
   \text{ and }
   x^T Aq\leq p^T Aq,\ \text{for all }x\in S^N
   \right\}.
   \tag{16}
$$

2. The strategy $p$ is an extremal equilibrium strategy if and only if

$$
   \begin{bmatrix}p\\p^T Ap\end{bmatrix}
$$

   is a vertex of $\mathcal{P}_I$.

**Proof.** For (a), see Ref. 6, Theorems 5 and 11. The proof of (b) is straightforward and can be found in Ref. 14. $\square$

Hence, the procedure FINDEQ$(I)$ can be reduced to the problem to find a vertex $\begin{bmatrix}p\\v\end{bmatrix}$ of $\mathcal{P}_I$ with

$$
v=p^T Ap
\qquad\text{and}\qquad
I(p)=I,
$$

which can be easily accomplished, e.g., by linear programming techniques. Now, it remains to specify the routine CHECKSTAB$(p)$. To this end, we need some further notation: for any index sets $J\subseteq N$ and $K\subseteq J$, let

$$
\mathbb{R}_+^J(K)
=
\{x\in\mathbb{R}^J:x_k\geq0,\ \text{for all }k\in K\}
\tag{17}
$$

be the cone in $\mathbb{R}^J$ which results by imposing sign restrictions on some coordinates.

### Theorem 3.2

Let $p\in S^N$ be an equilibrium strategy; pick an arbitrary $m\in I(p)$; and form $J'(p)=J(p)\setminus\{m\}$.

1. If $J'(p)=\varnothing$, then $p=e_m$ is a pure ESS.
2. Else, form the matrix $B=[b_{ij}]_{i,j\in J'(p)}$, where

$$
   b_{ij}
   =
   a_{mj}+a_{jm}+a_{im}+a_{mi}-a_{ij}-a_{ji}-2a_{mm},
   \qquad i,j\in J'(p),
   \tag{18}
$$

   and put $K=J(p)\setminus I(p)$. Then, $p$ is an ESS if and only if $B$ is strictly $\mathbb{R}_+^{J'(p)}(K)$-copositive; i.e.,

$$
   y^T By>0,
   \qquad
   \text{for all }y\in\mathbb{R}_+^{J'(p)}(K)\text{ with }y\neq o.
   \tag{19}
$$

**Proof.**

1. This follows from the fact that, for any equilibrium strategy $p$, the relation $J(p)=\{m\}$ means that $(p,p)=(e_m,e_m)$ is a strict (or strong, in Harsanyi's terminology) equilibrium point. Now, Theorem 14(i) in Ref. 6 yields the result.
2. A proof for a more general situation can be found in Ref. 12, pp. 95ff. A simpler version adapted to the present setting is specified in Ref. 14. $\square$

Since the only inputs for CHECKSTAB$(p)$ are in fact an element $m\in I(p)$ and the index set $J(p)$, this procedure might serve as well to determine whether or not an ESS $p$ for $A$ exists with $m\in I(p)$, $J(p)=J$, for given $m$ and $J$, without calling FINDEQ$(I)$ before.

The procedure CHECKSTAB, and thus the algorithm proposed, will be complete if it incorporates a routine for deciding whether or not a given symmetric matrix $B$ is strictly $\mathbb{R}_+^J(K)$-copositive or not. There are several algorithms for solving this problem for strict copositivity on general polyhedral cones (Refs. 15-19). However, the cone $\mathbb{R}_+^J(K)$ has a special structure which enables us to reduce calculational effort considerably by means of the following recursive criterion.

### Theorem 3.3

Let $B=[b_{ij}]_{i,j\in J}$ be a symmetric matrix. For $K\subset J$, let

$$
J\setminus K=\{i_1,\ldots,i_r\}.
$$

Put

$$
K(0)=J
\qquad\text{and}\qquad
b_{jk}^{(0)}=b_{jk},\quad j,k\in K(0).
$$

For $1\leq\nu\leq r$, define

$$
K(\nu)=K(\nu-1)\setminus\{i_\nu\}
$$

and

$$
b_{jk}^{(\nu)}
=
b_{i_\nu i_\nu}^{(\nu-1)}b_{jk}^{(\nu-1)}
-
b_{i_\nu j}^{(\nu-1)}b_{i_\nu k}^{(\nu-1)},
\qquad j,k\in K(\nu).
\tag{20}
$$

Then, $B$ is strictly $\mathbb{R}_+^J(K)$-copositive if and only if:

1. $b_{i_\nu i_\nu}^{(\nu-1)}>0$, for all $\nu\in\{1,\ldots,r\}$;
2. $B^{(r)}$ is strictly $\mathbb{R}_+^K(K)$-copositive (if $K=\varnothing$, this condition should be ignored).

**Proof.** For $i\in J\setminus K$, set

$$
b'_{jk}=b_{ii}b_{jk}-b_{ij}b_{ik},
$$

where $j,k\in J$ with $j\neq i$ and $k\neq i$. All that we have to show is that $B$ is strictly $\mathbb{R}_+^J(K)$-copositive if and only if $b_{ii}>0$ and $B'$ is strictly $\mathbb{R}_+^{J\setminus\{i\}}(K)$-copositive. Now, denote by

$$
x(i)=\sum_{j\neq i}x_j e_j\in\mathbb{R}^J
$$

the vector resulting from $x$ by setting its $i$th coordinate zero. Then,

$$
x^T Bx
=
b_{ii}x_i^2+2x_i e_i^T Bx(i)+[x(i)]^T Bx(i),
\tag{21}
$$

and $x\in\mathbb{R}_+^J(K)$ if and only if $x(i)\in\mathbb{R}_+^J(K)$. Therefore, $B$ is strictly $\mathbb{R}_+^J(K)$-copositive if and only if the function

$$
\varphi(t\mid x)
=
b_{ii}t^2+2t e_i^T Bx(i)+[x(i)]^T Bx(i)>0,
\qquad\text{for all }t\in\mathbb{R},
\tag{22}
$$

for all $x\in\mathbb{R}_+^J(K)\setminus\{o\}$. Since $e_i\in\mathbb{R}_+^J(K)\setminus\{o\}$, strict $\mathbb{R}_+^J(K)$-copositivity of $B$ yields $b_{ii}>0$. Furthermore, the function $\varphi(t\mid x)$ is positive for all $t\in\mathbb{R}$ if and only if $\varphi(\bar t\mid x)>0$, where

$$
\bar t=-e_i^T Bx(i)/b_{ii}
$$

denotes the minimizer of $\varphi(\,\cdot\mid x)$. Now, for $y\in\mathbb{R}_+^{J\setminus\{i\}}(K)\setminus\{o\}$ and

$$
x_j=
\begin{cases}
y_j,&j\neq i,\\
\bar t,&j=i,
\end{cases}
\tag{23}
$$

we get

$$
x\in\mathbb{R}_+^J(K)\setminus\{o\}
$$

and

$$
y^T B'y
=
b_{ii}[x(i)]^T Bx(i)-[e_i^T Bx(i)]^2
=
b_{ii}\varphi(\bar t\mid x)
=
b_{ii}x^T Bx.
$$

Thus, the equivalence stated above is true. $\square$

Typically, $K=J(p)\setminus I(p)$ will contain only a few elements; hence, the remaining matrix $B^{(r)}$ will be of low dimension, so that the question of strict $\mathbb{R}_+^K(K)$-copositivity of $B^{(r)}$ can be settled by use of recursive determinant criteria, e.g., in Ref. 17. Note that, if $K$ contains no element or only one element, then strict $\mathbb{R}_+^K(K)$-copositivity is positive definiteness. In this case, Theorem 3.3 is just a recursive version of Sylvester's minorant criterion for positive definiteness. If, however, $K$ contains many elements, then the recursive procedure in Ref. 19 may have some numerical advantages compared to the determinantal criteria as proposed, e.g., in Ref. 17. In any case, the worst-case computational complexity for this problem is the same as for a linear complementarity problem or a linear-quadratic optimization problem (as used in Ref. 10), which all are NP-complete (Ref. 20). These results stress the importance of an efficient search procedure avoiding unnecessary calls of CHECKSTAB. For the sake of completeness, we close this section with some explicit copositivity criteria in low dimensions, where $K$ contains more than one element.

### Theorem 3.4

Let $B=[b_{ij}]_{i,j\in J}$ be a symmetric matrix. For $K\subseteq J$, put $\mathcal{C}=\mathbb{R}_+^J(K)$.

1. If $J=K=\{1,2\}$, then $B$ is strictly $\mathcal{C}$-copositive if and only if:
   - $b_{ii}>0$, $1\leq i\leq2$; and
   - $\det B>0$ or $b_{12}\geq0$.

2. If $J=K=\{1,2,3\}$, then $B$ is strictly $\mathcal{C}$-copositive if and only if:
   - $b_{ii}>0$, $1\leq i\leq3$;
   - $b_{ij}>-\sqrt{b_{ii}b_{jj}}$, $1\leq i,j\leq3$; and
   - $\det B>0$ or

$$
     b_{12}\sqrt{b_{33}}
     +b_{23}\sqrt{b_{11}}
     +b_{13}\sqrt{b_{22}}
     +\sqrt{b_{11}b_{22}b_{33}}
     \geq0.
$$

3. If $K=\{1,2\}$, $J=\{1,2,3\}$, then $B$ is strictly $\mathcal{C}$-copositive if and only if:
   - $b_{33}>0$;
   - $b_{ii}b_{33}>(b_{i3})^2$, $1\leq i\leq2$; and
   - $\det B>0$ or $b_{12}b_{33}\geq b_{13}b_{23}$.

**Proof.** For (a) and (b), see Theorem 9 in Ref. 17 and the remarks preceding this theorem. It remains to show (c), where we apply Theorem 3.3 to get the conditions $b_{33}>0$ and strict $\mathbb{R}_+^K(K)$-copositivity of

$$
B'
=
\begin{bmatrix}
b_{11}b_{33}-(b_{13})^2 & b_{12}b_{33}-b_{13}b_{23}\\
b_{12}b_{33}-b_{23}b_{13} & b_{22}b_{33}-(b_{23})^2
\end{bmatrix}.
\tag{24}
$$

Now, apply (a1) to obtain

$$
b_{ii}b_{33}-(b_{i3})^2>0,
$$

while the first condition in (a2) for $B'$ reads

$$
\begin{aligned}
0<\det B'
&=
b_{33}\left[
b_{11}b_{22}b_{33}
-b_{11}(b_{23})^2
-b_{22}(b_{13})^2
-b_{12}^2b_{33}
+2b_{12}b_{13}b_{23}
\right]\\
&=b_{33}\det B.
\end{aligned}
\tag{25}
$$

Finally, the last condition of (a2) for $B'$ is equivalent to

$$
b_{12}b_{33}\geq b_{23}b_{13}.
$$

$\square$

## 4. Playing-the-Field Models

Up to now, we concentrated on evolutionary games consisting of pairwise conflicts, the payoff structure being determined by the mean payoff resulting from random encounters within the population. However, there are some situations (called "playing the field," e.g., in Ref. 21) where incremental fitness depends in a more complicated way on the average strategy $x$ within the population; cf. Ref. 12, Chapters 5 and 7. So we now assume the $i$-strategist's payoff to be given by $f_i(x)$, where $f:S^N\to\mathbb{R}^N$ is the mean payoff function. If $f(x)$ is linear in $x$, i.e., equal to $Ax$ for some payoff matrix $A$, then this corresponds to the pairwise contest case we considered in the preceding sections, while a nonlinear payoff function $f$ signifies the playing-the-field context, in which a strategy $p\in S^N$ is said to be evolutionarily stable whenever:

1. $p$ is a best reply on itself, i.e.,

$$
   x^T f(p)\leq p^T f(p),\qquad\text{for all }x\in S^N;
   \tag{26}
$$

2. $p$ fares better against any other reply $x$ on $p$, than $x$ does against itself:

$$
   \text{if }x^T f(p)=p^T f(p),\qquad\text{then }p^T f(x)>x^T f(x),
   \tag{27}
$$

   provided $x$ is of the form $x=(1-\epsilon)p+\epsilon q$, where $q\in S^N$ and with $0<\epsilon<1$ represents the relative amount of individuals with different behavior $q\neq p$ (perhaps due to mutation) within the overall population; the inequality above shall be valid for all $\epsilon$ with $0<\epsilon\leq\rho$, where $\rho>0$ may depend (but need not depend) on $q$.

Linearity of $f$ allows us to extend the locality condition in property (27) to the whole of $S^N\setminus\{p\}$ and to obtain the stability condition (4).[^2] The equilibrium condition (26) means that the payoff vector $f(p)$ belongs to the (polyhedral) normal cone of $S^N$ at $p$ (Ref. 12, p. 74),

$$
\begin{aligned}
\mathcal{N}_p
&=
\{y\in\mathbb{R}^N:(x-p)^T y\leq0,\ \text{for all }x\in S^N\}\\
&=
\{y\in\mathbb{R}^N:
y_i\leq p^T y,\ \text{for all }i\in N,\ \text{and }
y_i=p^T y,\ \text{for all }i\in I(p)\},
\end{aligned}
\tag{28}
$$

so that an analogue to the procedure FINDEQ$(I)$ has to generate serious candidates $p\in S^N$ with $I(p)=I$ and $f(p)\in\mathcal{N}_p$. Here, let us concentrate on the counterpart of CHECKSTAB$(p)$: due to nonlinearity, this analogue now relies on a first-order approximation; thus, there is a gap between necessary and sufficient conditions for stability.

### Theorem 4.1

Let $p\in S^N$ satisfy $f(p)\in\mathcal{N}_p$. Pick $m\in I(p)$; define

$$
J(p)=\{i\in N:f_i(p)=p^T f(p)\}\supseteq I(p);
\tag{29}
$$

denote $J'(p)=J(p)\setminus\{m\}$; and put $K=J(p)\setminus I(p)$.

1. If $J'(p)=\varnothing$, then $p=e_m$ is a pure ESS.
2. Otherwise, let

$$
   A_p=\left[\left(\frac{\partial f_i}{\partial x_j}\right)(p)\right]_{i,j\in N}
$$

   and form $B_p=[b_{ij}]_{i,j\in J'(p)}$ from $A_p$ by analogy to (18); then:
   - $p$ is evolutionarily stable if $B_p$ is strictly $\mathbb{R}_+^{J'(p)}(K)$-copositive;
   - if $p$ is evolutionarily stable, then $B_p$ is $\mathbb{R}_+^{J'(p)}(K)$-copositive, i.e.,

$$
     y^T B_p y\geq0,\qquad\text{for all }y\in\mathbb{R}_+^{J'(p)}(K).
$$

**Proof.** This follows from Theorem 35 and Proposition 40 in Ref. 12. $\square$

In order to reduce the gap between necessary and sufficient stability conditions, a procedure seems to be desirable which decides whether or not a symmetric matrix $B$ is $\mathbb{R}_+^J(K)$-copositive. As in the preceding section, also in this case a recursive reduction procedure may be of some value prior to application of the well-known copositivity algorithms in Refs. 15 to 19.

### Theorem 4.2

Let $B=[b_{ij}]_{i,j\in J}$ be a symmetric matrix, and suppose that $K\subset J$. For $0\leq\nu\leq r$, define $i_\nu$, $K(\nu)$, and $B^{(\nu)}$ as in Theorem 3.3.

1. If

$$
   b_{i_\nu i_\nu}^{(\nu-1)}\neq0,
   \qquad\text{for all }\nu\in\{1,\ldots,r\},
$$

   then $B$ is $\mathbb{R}_+^J(K)$-copositive if and only if:
   - $b_{i_\nu i_\nu}^{(\nu-1)}>0$, for all $\nu\in\{1,\ldots,r\}$; and
   - $B^{(r)}$ is $\mathbb{R}_+^K(K)$-copositive.

2. If, however, there is an $s\in\{1,\ldots,r\}$ such that

$$
   b_{i_\nu i_\nu}^{(\nu-1)}\neq0,
   \quad\text{for all }\nu\in\{1,\ldots,s-1\},
   \qquad\text{but}\qquad
   b_{i_s i_s}^{(s-1)}=0,
$$

   then $B$ is $\mathbb{R}_+^J(K)$-copositive if and only if:
   - $b_{i_\nu i_\nu}^{(\nu-1)}>0$, for all $\nu\in\{1,\ldots,s-1\}$;
   - $b_{i_s j}^{(s-1)}=0$, for all $j\in K(s-1)$; and
   - $C=[b_{ij}^{(s-1)}]_{i,j\in K(s)}$ is $\mathbb{R}_+^{K(s)}(K)$-copositive.

Again, condition (b3) shall be considered as void if $K(s)=K(r)=K=\varnothing$.

**Proof.**

1. This follows as in the proof of Theorem 3.3.
2. To prove (b), we proceed similarly by induction, so that all that we have to show is the following assertion: whenever $b_{ii}=0$, for some $i\in J\setminus K$, then $B$ is $\mathbb{R}_+^J(K)$-copositive if and only if $b_{ij}=0$, for all $j\in J\setminus\{i\}$, and $C=[b_{jk}]_{j,k\in J\setminus\{i\}}$ is $\mathbb{R}_+^{J\setminus\{i\}}(K)$-copositive. To this end, consider again the expression $\varphi(t\mid x)$ as in (22). Now, copositivity of $B$ yields via $\varphi(t\mid x)\geq0$, for all $t\in\mathbb{R}$, and $b_{ii}=0$, the relation

$$
   e_i^T Bx(i)=0,\qquad\text{for all }x\in\mathbb{R}_+^J(K),
$$

   in particular

$$
   b_{ij}=e_i^T Be_j=0,\qquad\text{for all }j\in J\setminus\{i\}.
$$

   Hence, from (21) we obtain

$$
   \varphi(t\mid x)=[x(i)]^T Bx(i)=y^T Cy,
$$

   where for any $x\in\mathbb{R}_+^J(K)$ we put

$$
   y=[x_j]_{j\in J\setminus\{i\}}\in\mathbb{R}_+^{J\setminus\{i\}}(K).
   \tag{30}
$$

   Since every vector in $\mathbb{R}_+^{J\setminus\{i\}}(K)$ can be obtained in this way, we conclude that $C$ has to be $\mathbb{R}_+^{J\setminus\{i\}}(K)$-copositive. Conversely, copositivity of $C$ and the relation $b_{ij}=0$, for all $j\in J\setminus\{i\}$, together yield via $\varphi(t\mid x)=y^T Cy\geq0$ copositivity of $B$. $\square$

Also if case (b) occurs, Theorem 4.2 exhibits a recursive method to detect copositivity, since the order of $C$ is smaller than the order of $B$, because $K(s)$ does not contain $i_s$. To keep the order of $C$ small, it is reasonable to impose maximality of $s$. If the search for an arrangement of the indices in $J\setminus K=\{i_1,\ldots,i_r\}$ yielding a maximal $s$ turns out to be too long to be efficient, one might proceed inductively as follows: choose an index $i_1\in K(0)=J$ such that $b_{i_1i_1}\neq0$. If, for some $\nu$ with $1\leq\nu\leq r$, the quantities $K(\nu-1)$ and $B^{(\nu-1)}$ are already determined, then choose $i_\nu\in K(\nu)$ such that $b_{i_\nu i_\nu}^{(\nu-1)}\neq0$.[^3] In the case where

$$
b_{ii}^{(\nu-1)}=0,\qquad\text{for all }i\in K(\nu-1),
$$

we have $s=\nu-1$ and proceed as in Theorem 4.2(b).[^4] According to its proof, the order of $C$ can be reduced further by omitting the $i$th row and column if

$$
b_{ii}^{(s-1)}=0,\qquad\text{for some }i\in K(s)\setminus K.
$$

Finally, we again specify some explicit criteria in low dimensions.

### Theorem 4.3

Let $B=[b_{ij}]_{i,j\in J}$ be a symmetric matrix. For $K\subseteq J$, put $\mathcal{C}=\mathbb{R}_+^J(K)$.

1. If $J=K=\{1,2\}$, then $B$ is $\mathcal{C}$-copositive if and only if:
   - $b_{ii}\geq0$, $1\leq i\leq2$; and
   - $\det B\geq0$ or $b_{12}\geq0$.

2. If $J=K=\{1,2,3\}$, then $B$ is $\mathcal{C}$-copositive if and only if:
   - $b_{ii}\geq0$, $1\leq i\leq3$;
   - $b_{ij}\geq-\sqrt{b_{ii}b_{jj}}$, $1\leq i,j\leq3$; and
   - $\det B\geq0$ or

$$
     b_{12}\sqrt{b_{33}}
     +b_{23}\sqrt{b_{11}}
     +b_{13}\sqrt{b_{22}}
     +\sqrt{b_{11}b_{22}b_{33}}
     \geq0.
$$

3. If $K=\{1,2\}$, $J=\{1,2,3\}$, then $B$ is $\mathcal{C}$-copositive if and only if either
   - $b_{33}>0$;
   - $b_{ii}b_{33}\geq(b_{i3})^2$, $1\leq i\leq2$; and
   - $\det B\geq0$ or $b_{12}b_{33}\geq b_{13}b_{23}$;

   or
   - $b_{13}=b_{23}=b_{33}=0$;
   - $b_{ii}\geq0$, $1\leq i\leq2$; and
   - $b_{11}b_{22}\geq(b_{12})^2$ or $b_{12}\geq0$.

**Proof.** The assertions (a) and (b) are again taken from Ref. 17. To prove (c), apply Theorem 4.2 and distinguish the cases $b_{33}\neq0$ or $b_{33}=0$. If $b_{33}\neq0$, then $b_{33}>0$ and, as in the proof of Theorem 3.4, conditions (a1) and (a2) applied to $B^{(1)}=B'$ instead of $B$ yield (c2) and (c3). If, however, $b_{33}=0$, conditions (c2') and (c3') are (a1) and (a2) applied to $C=[b_{ij}]_{1\leq i,j\leq2}$. $\square$

## 5. Conclusions

Utilizing the support criterion (Proposition 2.1), an algorithm for finding all ESSs in a given game is proposed which is both efficient and exact in Section 2. The method incorporates two main subroutines which are described in Section 3: FINDEQ generates serious candidates for ESSs using extremality properties (Theorem 3.1); and CHECKSTAB checks the stability condition by rephrasing it into a strict copositivity condition in Theorem 3.2. Theorem 3.3 helps to reduce problem dimension recursively, while Theorem 3.4 exhibits explicit copositivity criteria in low dimensions. Section 4 deals with the generalization of these results for playing-the-field models where payoff depends nonlinearly on the state. This nonlinearity results in a bifurcation of necessary and sufficient conditions for evolutionary stability specified in Theorem 4.1, where the necessary part is related to copositivity. Procedures for checking this weaker property are described in Theorems 4.2 and 4.3, which are considerably more complex than their counterparts in Theorems 3.3 and 3.4.

## References

1. Maynard Smith, J., and Price, G., "The Logic of Animal Conflict," *Nature*, Vol. 246, pp. 15-18, 1973.
2. Maynard Smith, J., "The Theory of Games and the Evolution of Animal Conflict," *Journal of Theoretical Biology*, Vol. 47, pp. 209-221, 1974.
3. Stewart, F. M., "Evolution of Dimorphism in a Predator-Prey Model," *Theoretical Population Biology*, Vol. 2, pp. 493-506, 1971.
4. Bomze, I. M., and van Damme, E. E. C., "A Dynamical Characterization of Evolutionary Stability," *Annals of Operations Research*, Vol. 37, pp. 229-244, 1992.
5. Bomze, I. M., "Dynamical Aspects of Evolutionary Stability," *Monatshefte für Mathematik*, Vol. 110, pp. 189-206, 1990.
6. Bomze, I. M., "Noncooperative Two-Person Games in Biology: A Classification," *International Journal of Game Theory*, Vol. 15, pp. 31-57, 1986.
7. van Damme, E. E. C., *Stability and Perfection of Nash Equilibria*, Springer, Berlin, Germany, 1987.
8. Haigh, J., "Game Theory and Evolution," *Advances of Applied Probability*, Vol. 7, pp. 8-11, 1975.
9. Bishop, D. T., and Cannings, C., "Models of Animal Conflict," *Advances of Applied Probability*, Vol. 8, pp. 616-621, 1976.
10. Williams, H. P., "Evolution, Game Theory, and Polyhedra," *Journal of Mathematical Biology*, Vol. 25, pp. 393-409, 1987.
11. Abakuks, A., "Conditions for Evolutionarily Stable Strategies," *Journal of Applied Probability*, Vol. 17, pp. 559-562, 1980.
12. Bomze, I. M., and Pötscher, B. M., *Game Theoretic Foundations of Evolutionary Stability*, Springer, Berlin, Germany, 1989.
13. Parthasarathy, T., and Raghavan, T. E. S., *Some Topics in Two-Person Games*, Elsevier, New York, New York, 1971.
14. Bomze, I. M., *Detecting All Evolutionarily Stable Strategies*, Technical Report No. 88, Institut für Statistik und Informatik, Universität Wien, 1990.
15. Diananda, P. H., "On Nonnegative Forms in Real Variables Some or All of Which Are Nonnegative," *Proceedings of the Cambridge Philosophical Society*, Vol. 58, pp. 17-25, 1962.
16. Cottle, R. W., Habetler, G. J., and Lemke, C. E., "Quadratic Forms Semi-Definite over Convex Cones," *Proceedings of the Princeton Symposium on Mathematical Programming*, edited by H. W. Kuhn, Princeton University Press, Princeton, New Jersey, pp. 551-565, 1970.
17. Hadeler, K. P., "On Copositive Matrices," *Linear Algebra and Applications*, Vol. 49, pp. 79-89, 1983.
18. Bomze, I. M., "Remarks on the Recursive Structure of Copositivity," *Journal of Informational and Optimization Sciences*, Vol. 8, pp. 243-260, 1987.
19. Danninger, G., "A Recursive Algorithm for Determining (Strict) Copositivity of a Symmetric Matrix," *Methods of Operations Research*, Hain, Meisenheim, Germany, Vol. 62, pp. 45-52, 1990.
20. Murty, K. G., and Kabadi, S. N., "Some NP-Complete Problems in Quadratic and Nonlinear Programming," *Mathematical Programming*, Vol. 39, pp. 117-129, 1987.
21. Maynard Smith, J., *Evolution and the Theory of Games*, Cambridge University Press, Cambridge, England, 1982.

[^2]: Equation (4) is the equilibrium condition in Section 1; equation (5) is the stability condition. The printed reference to the "stability condition (4)" therefore appears to mean condition (5).
[^3]: The printed paper says $i_\nu\in K(\nu)$. Because $K(\nu)=K(\nu-1)\setminus\{i_\nu\}$, this appears to be a typographical error; the recursive definition requires choosing $i_\nu\in K(\nu-1)$.
[^4]: If every prospective pivot in $B^{(\nu-1)}$ is zero, Theorem 4.2(b) applies with its first zero pivot at $s=\nu$, not $s=\nu-1$. The printed index appears to be off by one.
