# On Copositive Matrices

**Author:** K. P. Hadeler<br>
**Affiliation:** Lehrstuhl für Biomathematik der Universität Tübingen, Institut für Biologie II, Auf der Morgenstelle 28,
7400 Tübingen 1, West Germany<br>
**Journal:** *Linear Algebra and Its Applications*, Vol. 49, 1983, pp. 79-89<br>
**Submitted by:** Hans Schneider<br>
**Received:** 8 September 1981; revised 12 April 1982<br>
**Source PDF:** [Hadeler_1983.pdf](Hadeler_1983.pdf)

> Complete text transcription of the paper. Mathematical notation has been reconstructed as LaTeX from the printed pages.
> Journal running headers, page numbers, and repeated publication footers are omitted.

## Abstract

A criterion for copositive matrices is given and for \(n=3\) the set of all copositive matrices is determined in terms of
matrix elements. Copositive matrices are applied to the problem of excluding periodic solutions of certain algebraic
differential equations.

## Introduction

Since Motzkin [31] introduced the notion of a copositive matrix, the original results have been extended in several
directions: characterization of copositive matrices [6, 12, 35], extension to quadratic programming [10, 11, 21, 22, 23,
29] and related combinatorial problems [1, 2, 3, 14, 15, 20], and extension to spaces of infinite dimension [26]. However,
it is still difficult to determine whether a given matrix is copositive; and there are only few useful applications. In the
following I shall derive some additional criteria for copositivity, and find all copositive matrices for dimension 3.
Furthermore the concept of copositivity will be applied to algebraic differential equations.

First I give an account of definitions and earlier results. Let \(\mathbb{R}^n\) be the usual coordinate space, and let
\(K\) be the closed convex cone of nonnegative vectors. Consider the space \(H\) of real symmetric matrices
\(A=(a_{jk})\) of order \(n\). Its dimension is \(n(n+1)/2\). Among the many cones in this space the following four are of
particular interest (\(T\) denotes transpose):

\[
P=\{A:a_{jk}\geq 0\text{ for }j,k=1,\ldots,n\},
\]

\[
S=\{A:x^T Ax\geq 0\text{ for }x\in\mathbb{R}^n\},
\]

\[
C=\{A:x^T Ax\geq 0\text{ for }x\in K\},
\]

\[
B=
\left\{
A:A=\sum_{j=1}^{m}x_jx_j^T,\ x_j\in K,\ m\text{ finite}
\right\}.
\]

\(P\) is the cone of positive matrices in the sense of Perron, \(S\) are the positive definite matrices, the matrices in
\(C\) are copositive, and those in \(B\) are the completely positive matrices. Each of these cones is convex, closed,
pointed, and has nonempty interior. A matrix \(A\) with \(x^T Ax>0\) for \(x\in K\), \(x\neq 0\), is called strictly
copositive.

A natural inner product in \(H\) is given by \(\operatorname{tr}XY\). Denote the dual cone by an asterisk. Obviously
\(P=P^*\). From the identity \(\operatorname{tr}Axx^T=x^T Ax\) one sees \(S=S^*\), \(B=C^*\), and thus \(B^*=C\).
Always \(P+S\subset C\); thus by duality, \(P\cap S\supset B\). For \(n=2\) one has \(P\cup S=C\); for \(n=3,4\) one
has \(P+S=C\) (Diananda [9]); for \(n=5\) there are elements in \(P\cap S\) which are not in \(B\) (examples due to Hall
and Horn are given in [14, p. 265 ff.]), and hence \(P+S\neq C\). A detailed survey on copositive matrices has been given
by M. Hall [14], in particular on the unsolved problem of determining the extremal rays of the cone \(C\) (those of
\(P,S,B\) are easily found). Many results on extremal rays are contained in [1, 2, 3, 4, 14, 15, 20]. An interesting
relation between copositivity with respect to a general cone and the Perron property has been shown in [18]. In a series of
recent papers (see [29]) quadratic forms with general linear restrictions are considered. In a fundamental paper Cottle,
Habetler, and Lemke [6] have given an inductive criterion for copositivity. In the following Theorem 1 I show a more
general criterion (actually a whole family of criteria), which also provides a coordinate-free proof of their result
(Theorem 2).

The sign \(\geq\) indicates that a vector or matrix has nonnegative elements; the sign \(>\) says that all elements are
positive.

### Definition 1

Let \(A\) be a symmetric matrix, and let \(B\) be a strictly copositive matrix of the same order. The pair \(A,B\) is called
(strictly) codefinite iff

\[
Ax=\lambda Bx,\quad x>0
\]

implies \(\lambda\geq 0\) (\(\lambda>0\), respectively).

### Definition 2

The symmetric matrix \(A\) of order \(n\) is called (strictly) copositive of order \(m\), \(1\leq m\leq n\), iff every
principal submatrix of order \(m\) is (strictly) copositive.

## Results

### Theorem 1

Let \(A\) of order \(n\) be (strictly) copositive of order \(n-1\). The following statements are equivalent:

1. \(A\) is (strictly) copositive.
2. There is a strictly copositive matrix \(B\) of order \(n\) such that the pair \(A,B\) is (strictly) codefinite.
3. For every strictly copositive matrix \(B\) of order \(n\) the pair \(A,B\) is (strictly) codefinite.

**Proof.** \((1)\Rightarrow(3)\): Let \(B\) be strictly copositive of order \(n\). Suppose \(Ax=\lambda Bx\), \(x>0\).
Then \(x^T Ax=\lambda x^T Bx\), and thus \(\lambda\geq 0\) (or \(\lambda>0\), respectively). Since
\((3)\Rightarrow(2)\) is obvious, it remains to show \((2)\Rightarrow(1)\). Suppose the form \(x^T Ax\) assumes negative
(nonpositive, respectively) values for certain \(x\geq 0\), \(x\neq 0\). Since \(A\) is (strictly) copositive of order
\(n-1\), such \(x\) are necessarily positive. Consider the form \(x^T Ax\) on the set

\[
D=\{x\geq 0:x^T Bx=1\}.
\]

It assumes its minimum at a point \(\bar{x}\) of the relative interior of \(D\). The minimum is negative (nonpositive).
From the necessary condition for a minimum of \(x^T Ax\) on \(D\) it follows that
\(A\bar{x}=\lambda B\bar{x}\) with \(\lambda=\bar{x}^T A\bar{x}\). In view of statement (2), \(\lambda\geq 0\) (or
\(\lambda>0\)), which leads to a contradiction. \(\square\)

From Theorem 1 one can easily proceed to Theorem 2, which is essentially the result of Cottle, Habetler, and Lemke [6].
For convenience the theorem is stated in a negative form.

### Theorem 2

Let \(A\) of order \(n\) be copositive of order \(n-1\). Then the following statements are equivalent:

1. \(A\) is not copositive.
2. For every \(b>0\) there is an \(x>0\) such that \(Ax=\lambda b\), \(\lambda<0\).
3. The matrix \(-A^{-1}\) exists and is nonnegative.
4. \(\det A<0\) and \(\operatorname{adj}A\geq 0\).

Here \(\det A\) is the determinant and \(\operatorname{adj}A\) is the adjugate, i.e. the matrix of signed cofactors, such
that

\[
A\operatorname{adj}A=(\det A)I.
\]

**Proof.** For the matrix \(B\) in Theorem 1 one can choose the identity or matrices \(B=bb^T\) of rank 1, where \(b>0\)
is a positive vector. Then the equivalence of (1) and (2) follows from Theorem 1 by contraposition. Assume (2), and suppose
there is \(y\) with \(Ay=0\). Since \(Ax=b\) has a solution for \(b>0\), it follows \(y^T b=0\) for \(b>0\); thus
\(y=0\). Hence \(A^{-1}\) exists. For every \(b>0\) the vector \(x=A^{-1}b\) is elementwise negative. Thus \(-A^{-1}\)
is a nonnegative matrix, which shows (3). Of course from (3) follows (2).

From Theorem 1, choosing \(B=I\), it follows that there is \(x>0\) such that \(Ax=\lambda x\), \(\lambda<0\). Suppose
\(A\) has a second negative eigenvalue, \(Ay=\mu y\), \(y^T x=0\), \(y\neq 0\), \(\mu<0\). The vector \(y\) has
positive and negative components in view of \(x>0\), \(y^T x=0\). Choose \(\alpha>0\) such that \(x+\alpha y\geq 0\),
and \(x+\alpha y\) has at least one zero component. Since \(A\) is copositive of order \(n-1\),

\[
0\leq (x+\alpha y)^T A(x+\alpha y)
=\lambda x^T x+\mu\alpha^2y^T y<0.
\]

Thus \(A\) has exactly one negative eigenvalue, and \(\det A<0\), which shows the equivalence with (4). \(\square\)

For completeness I include the corresponding result on strictly copositive matrices. The proof basically follows [6], but
is somewhat shorter.

### Theorem 3

Let \(A\) of order \(n\) be strictly copositive of order \(n-1\). Then the following statements are equivalent:

1. \(A\) is not strictly copositive.
2. For every \(b>0\) there is an \(x>0\) such that \(Ax=\lambda b\), \(\lambda\leq 0\).
3. \(\det A\leq 0\) and \(\operatorname{adj}A>0\).

**Proof.** The equivalence of (1) and (2) follows from Theorem 1. Now assume (1) and (2).

*Case 1: \(A\) is not copositive.* Then \(\det A<0\) and \(-A^{-1}\geq 0\) in view of Theorem 2. For any vector \(u\)
and \(v=(\operatorname{adj}A)u\),

\[
v^T Av=(\det A)u^T(\operatorname{adj}A)u.
\tag{*}
\]

If \(\operatorname{adj}A\) is not elementwise positive, then there is a vector \(u\geq 0\), \(u\neq 0\), such that
\(v\geq 0\), \(v\neq 0\) is not positive. By strict copositivity of order \(n-1\) it follows that \(v^T Av>0\), whereas
\((\det A)u^T(\operatorname{adj}A)u\leq 0\).

*Case 2: \(A\) is copositive.* By Theorem 1, \(Ax=\lambda b\), \(b>0\), \(x>0\) implies \(\lambda\geq 0\). At least
for one \(x>0\) it must occur that \(\lambda=0\) (thus \(\det A=0\)); otherwise \(A\) would be strictly copositive.
Suppose there is a second nonpositive eigenvalue, \(Ay=\mu y\), \(\mu\leq 0\), \(y\neq 0\), \(y^T x=0\). Choose
\(\alpha\) such that \(x+\alpha y\geq 0\) with at least one zero component. Then

\[
0<(x+\alpha y)^T A(x+\alpha y)
=\lambda x^T x+\mu\alpha^2y^T y\leq 0.
\]

Thus \(A\) is positive semidefinite of rank \(n-1\). For any semidefinite matrix \(A\) of rank \(n-1\) the adjugate is
\(\operatorname{adj}A=xx^T\), where \(x\) is a (normalized) eigenvector for \(\lambda=0\). Thus
\(\operatorname{adj}A>0\).

These arguments yield statement (3). Now suppose (3) and that \(A\) is strictly copositive. Then from \((*)\) a
contradiction immediately follows. \(\square\)

In [6] it has been observed that the determinant criteria of Motzkin [32, 33] and Garsia [13] can be derived from Theorems 2
and 3.

From Theorem 2 it would appear that, among the matrices which are copositive of order \(n-1\), the copositive matrices are
characterized by a large set of inequalities, corresponding to the cofactors. In the cases \(n=2\) and \(n=3\) there is in
fact only one inequality in addition to \(\det A\geq 0\). The result for \(n=2\) is well known. The matrix
\(A=(a_{jk})\) is copositive iff \(a_{11}\geq 0\), \(a_{22}\geq 0\), and

\[
a_{11}a_{22}-a_{12}^2\geq 0
\qquad\text{or}\qquad
a_{12}\geq 0.
\]

It is strictly copositive iff \(a_{11}>0\), \(a_{22}>0\), and

\[
a_{11}a_{22}-a_{12}^2>0
\qquad\text{or}\qquad
a_{12}\geq 0.
\]

The case \(n=3\) is covered in the following theorem.

### Theorem 4

Let \(n=3\). The matrix \(A\) is copositive if and only if the conditions

\[
\begin{aligned}
\text{(1)}\quad &a_{11}\geq 0,\quad a_{22}\geq 0,\quad a_{33}\geq 0,\\
\text{(2)}\quad &a_{12}\geq-\sqrt{a_{11}a_{22}},\quad
a_{23}\geq-\sqrt{a_{22}a_{33}},\quad
a_{31}\geq-\sqrt{a_{33}a_{11}}
\end{aligned}
\]

are satisfied, as well as at least one of the following conditions:

\[
\text{(3a)}\quad
a_{12}\sqrt{a_{33}}+a_{23}\sqrt{a_{11}}+a_{31}\sqrt{a_{22}}+\sqrt{a_{11}a_{22}a_{33}}\geq 0.
\]

\[
\text{(3b)}\quad \det A\geq 0.
\]

The matrix is strictly copositive if and only if these conditions are satisfied with strict inequality in (1), (2), (3b).

**Proof.** First assume \(a_{11}>0\), \(a_{22}>0\), \(a_{33}>0\). By diagonal scaling we can put the matrix into the form

\[
A=
\begin{pmatrix}
1 & \alpha & \beta\\
\alpha & 1 & \gamma\\
\beta & \gamma & 1
\end{pmatrix}.
\tag{1}
\]

Now the necessary and sufficient conditions for copositivity of order two are

\[
\alpha\geq-1,\qquad \beta\geq-1,\qquad \gamma\geq-1.
\tag{2}
\]

Under this hypothesis the necessary and sufficient condition for \(A\) not to be copositive is

\[
\det A<0
\qquad\text{and}\qquad
\operatorname{adj}A\geq 0,
\]

i.e.,

\[
1+2\alpha\beta\gamma<\alpha^2+\beta^2+\gamma^2
\tag{3}
\]

and

\[
\alpha^2\leq 1,\qquad \beta^2\leq 1,\qquad \gamma^2\leq 1,
\tag{4a}
\]

\[
\alpha\beta\geq\gamma,\qquad \beta\gamma\geq\alpha,\qquad \alpha\gamma\geq\beta.
\tag{4b}
\]

The assertion of Theorem 4 claims that \(A\) is not copositive iff (3) is satisfied and

\[
\alpha+\beta+\gamma+1<0.
\tag{5}
\]

One has to show that, under the hypothesis of (2), (3), the conditions (4) and (5) are equivalent.

*Case 1: Let \(\alpha,\beta,\gamma\geq 0\).* From (4a) it follows that \(\alpha,\beta,\gamma\in[0,1]\), and from (4b) it
follows that \(\alpha\beta\gamma\geq\alpha^2,\beta^2,\gamma^2\); thus
\(3\alpha\beta\gamma\geq\alpha^2+\beta^2+\gamma^2\). From (3) it follows that
\(1+2\alpha\beta\gamma<3\alpha\beta\gamma\); thus \(\alpha\beta\gamma>1\), in contradiction to
\(\alpha,\beta,\gamma\in[0,1]\). Hence (4) is empty. Of course (5) is empty.

*Case 2: \(\alpha,\beta\geq 0\), \(\gamma<0\).* From (4b) follows
\(0\geq\beta\gamma\geq\alpha\geq 0\), thus \(\alpha=\beta=0\). From (3) it follows that \(\gamma^2>1\), which
contradicts (4a). Hence (4) is empty. But (5) is also empty in view of
\(1\leq 1+\alpha+\beta<-\gamma\), in contradiction to (2).

*Case 3: \(\alpha>0\), \(\beta,\gamma\leq 0\).* Then \(\beta,\gamma\in[-1,0]\). First we show
\((4)\Rightarrow(5)\). We add the inequality (3) and the inequality
\(2(1+\alpha)(\alpha-\beta\gamma)\leq 0\) from (4b) to obtain
\((1+\alpha)^2<(\beta+\gamma)^2\), from which follows \(1+\alpha<-(\beta+\gamma)\). Now assume (2), (3), (5). From
\(\alpha<-\beta-\gamma-1\leq-\beta\) and \(\beta,\gamma\in[-1,0]\) it follows that \(\alpha\leq 1\) and
\(\alpha\gamma\geq\beta\). Similarly \(\alpha\beta\geq\gamma\). It remains to show \(\alpha\leq\beta\gamma\). Suppose
\(\alpha>\beta\gamma\). Then from (3) it follows that

\[
\alpha>\beta\gamma+\sqrt{(1-\beta^2)(1-\gamma^2)}=: \alpha_0.
\]

From (5), \(\alpha<-\beta-\gamma-1=: \alpha_1\). But \(\alpha_1\leq\alpha_0\) in view of

\[
-(1+\beta)(1+\gamma)\leq\sqrt{(1-\beta^2)(1-\gamma^2)}.
\]

Thus (4) is proved.

*Case 4: \(\alpha,\beta,\gamma\leq 0\).* From (2) it follows that \(\alpha,\beta,\gamma\in[-1,0]\), and (4) is trivial.
It remains to show that (5) is a consequence of (3). Suppose \(\alpha+\beta+\gamma\geq-1\). Then
\((\alpha+\beta+\gamma)^2\leq 1\); thus

\[
1-(\alpha^2+\beta^2+\gamma^2)
\geq 2(\alpha\beta+\beta\gamma+\gamma\alpha)
\geq-2\alpha\beta\gamma,
\]

which contradicts (3). \(\square\)

If some of the diagonal elements are zero, then the matrix can be diagonally scaled as in (1) with some of the 1's replaced
by 0's. Then the proof of equivalence is very simple.

The proof of the assertion on strict copositivity is verbally the same; one has to pay careful attention to equality signs.

## Application: Quadratic Differential Equations

A quadratic differential equation is an equation of the form

\[
\dot y=f(y),\qquad f:\mathbb{R}^n\longrightarrow\mathbb{R}^n,
\tag{6}
\]

where

\[
f_i(y)=\sum_{j,k=1}^{n}b_{ijk}y_jy_k.
\tag{7}
\]

Equations of this type occur in many applications, e.g. in population genetics. I assume that the \(n^3\) coefficients are
nonnegative. Then the cone \(K\) is positively invariant with respect to the equations (6). I introduce relative frequencies
by

\[
x=\frac{y}{e^T y},\qquad e^T=(1,\ldots,1),
\tag{8}
\]

and obtain, after rescaling the time variable, a differential equation

\[
\dot x=f(x)-e^T f(x)\cdot x.
\tag{9}
\]

The cone \(K\) as well as the simplex

\[
T=\{x\in K:e^T x=1\}
\]

is positively invariant with respect to the system (9). The Jacobian of the system (9) on \(K\) is

\[
J(x)=f'(x)-xe^T f'(x)-e^T f(x)I.
\tag{10}
\]

For \(x\in T\) the Jacobian has the left eigenvector \(e^T\),

\[
e^T J(x)=-e^T f(x)\cdot e^T.
\tag{11}
\]

Thus the divergence of the function \(f(x)-e^T f(x)\cdot x\), considered as a vector field on the simplex \(T\), is

\[
\begin{aligned}
D_0(x)
&=\operatorname{tr}J(x)+e^T f(x)\\
&=\operatorname{tr}f'(x)-e^T f'(x)x-(n-1)e^T f(x).
\end{aligned}
\]

On the set \(T\) the function \(D_0\) coincides with the function \(D:K\to\mathbb{R}\),

\[
D(x)=e^T x\operatorname{tr}f'(x)-e^T f'(x)x-(n-1)e^T f(x).
\tag{12}
\]

The function \(D\) is quadratic in \(x\). It is not the divergence of the vector field (9) on \(K\). In coordinate notation

\[
D(x)=
\sum_{j,k=1}^{n}
\left(
\sum_{i=1}^{n}b_{iji}
+\sum_{i=1}^{n}b_{iik}
-(n+1)\sum_{i=1}^{n}b_{ijk}
\right)x_jx_k.
\tag{13}
\]

Now I consider the case \(n=3\). Then the simplex \(T\) is planar. By the criterion of Dulac (or the negative criterion of
Bendixson) an autonomous differential system in the plane does not have periodic solutions (except constants) in a given
simply connected domain if the divergence does not change sign. Thus \(D(x)>0\) on \(T\) [or \(D(x)<0\) on \(T\)]
excludes the existence of periodic orbits of the system (9).

The function \(D\) is a quadratic form on the cone \(K\). One can apply the results of the previous section in order to
decide whether \(D\) changes sign. In particular, for \(n=3\) Theorem 4 provides the exact conditions on the \(b_{ijk}\).

There are vector fields for which this criterion applies, e.g. \(b_{ijk}=1\), for which \(D(x)=-6\). The question whether a
system of the form (9) can have periodic solutions at all can be easily answered. Consider a system of high symmetry
depending on only 6 parameters (instead of 18):

\[
\dot y_1
=ay_1^2+by_2^2+cy_3^2+2dy_1y_2+2ey_2y_3+2fy_3y_1,
\]

\[
\dot y_2
=cy_1^2+ay_2^2+by_3^2+2fy_1y_2+2dy_2y_3+2ey_3y_1,
\]

\[
\dot y_3
=by_1^2+cy_2^2+ay_3^2+2ey_1y_2+2fy_2y_3+2dy_3y_1.
\]

Straightforward computation shows that the Jacobian at the stationary point \(e/3\) is

\[
J=
\frac{2}{3}
\begin{pmatrix}
\alpha & \beta & \gamma\\
\gamma & \alpha & \beta\\
\beta & \gamma & \alpha
\end{pmatrix}
-\frac{2\chi}{9}ee^T
-\frac{\chi}{3}I,
\]

where

\[
\alpha=a+d+f,\qquad
\beta=b+d+c,\qquad
\gamma=c+e+f,\qquad
\chi=\alpha+\beta+\gamma.
\][^1]

The eigenvalues are \(\lambda_0=-\chi/3\) with eigenvector \(e^T\) and

\[
\lambda_{1,2}
=\frac{\alpha-2\beta-2\gamma}{3}
\pm\frac{\beta-\gamma}{\sqrt{3}}i.
\]

Thus for \(\beta\neq\gamma\), i.e. \(b+d\neq c+f\), and varying \(\alpha\) (e.g. \(a\)) a Hopf bifurcation occurs. A
simple choice is \(b=1\), \(c=d=e=f=0\),

\[
\lambda_{1,2}=\frac{a-2}{3}\pm\frac{i}{\sqrt{3}}.
\]

*The author wishes to thank the referee for helpful comments.*

## References

1. V. J. D. Baston, "Extreme copositive quadratic forms," *Acta Arith.* 15:319-327 (1969).
2. L. D. Baumert, "Extreme copositive quadratic forms," *Pacific J. Math.* 19:197-204 (1966).
3. L. D. Baumert, "Extreme copositive quadratic forms-II," *Pacific J. Math.* 20:1-20 (1967).
4. A. Berman, *Cones, Matrices, and Mathematical Programming*, Lecture Notes in Economics and Mathematical Systems 79,
   Springer, 1973.
5. A. Berman and R. J. Plemmons, *Non-negative Matrices in the Mathematical Sciences*, Academic, 1979.
6. R. W. Cottle, G. J. Habetler, and C. E. Lemke, "Quadratic forms semi-definite over convex cones," in *Proceedings of
   the Princeton Symposium on Mathematical Programming* (H. W. Kuhn, Ed.), Princeton U. P., Princeton, N.J., 1970,
   pp. 551-565.
7. R. W. Cottle, "Manifestations of the Schur complement," *Linear Algebra Appl.* 8:189-211 (1974).
8. R. W. Cottle, G. J. Habetler, and C. E. Lemke, "On classes of copositive matrices," *Linear Algebra and Appl.*
   3:295-310 (1970).
9. P. H. Diananda, "On nonnegative forms in real variables some or all of which are non-negative," *Proc. Cambridge
   Philos. Soc.* 58:17-25 (1962).
10. R. W. Farebrother, "Necessary and sufficient conditions for a quadratic form to be positive whenever a set of linear
    constraints is satisfied," *Linear Algebra Appl.* 16:39-42 (1977).
11. J. W. Gaddum, "Linear inequalities and quadratic forms," *Pacific J. Math.* 8:411-414 (1958).
12. A. Garsia, "Remarks about copositive forms," working paper, Dept. of Math., Univ. of California, San Diego, 1964.
13. A. Garsia, in [14], p. 274; [6], Theorem 4.1.
14. M. Hall, Jr., *Combinatorial Theory*, Blaisdell, Waltham, Mass., 1967.
15. M. Hall, Jr., in *Discrete Problems, A Survey of Numerical Analysis* (J. Todd, Ed.), New York, 1962, pp. 518-542.
16. M. Hall, Jr., and M. Newman, "Copositive and completely positive quadratic forms," *Proc. Cambridge Philos. Soc.*
    59:329-339 (1963).
17. H. Hancock, *Theory of Maxima and Minima*, Ginn, 1917; Dover, New York, 1960.
18. E. Haynsworth and A. J. Hoffman, "Two remarks on copositive matrices," *Linear Algebra Appl.* 2:387-392 (1969).
19. M. R. Hestenes and E. J. MacShane, "A theorem on quadratic forms and application in the calculus of variations,"
    *Trans. Amer. Math. Soc.* 47:501-512 (1940).
20. A. J. Hoffman and F. Pereira, "On copositive matrices with \(-1,0,1\) entries," *J. Combin. Theory Ser. A*
    14:302-309 (1973).
21. D. H. Jacobson, *Extensions of Linear Quadratic Control, Optimization, and Matrix Theory*, Academic, New York, 1977.
22. D. H. Jacobson, "A generalization of Finsler's theorem for quadratic inequalities and equalities," *Quaestiones Math.*
    1:19-28 (1976).
23. E. L. Keller, "Quadratic optimization and linear complementarity," Ph.D. Thesis, Univ. of Michigan, Ann Arbor, 1969.
24. R. Loewy and H. Schneider, "Positive operators on the \(n\)-dimensional ice-cream cone," *J. Math. Anal. Appl.*
    49:375-392 (1975).
25. T. L. Markham, "Factorizations of completely positive matrices," *Proc. Cambridge Philos. Soc.* 69:53-58 (1971).
26. D. H. Martin, "Conditional positivity of quadratic forms in Hilbert space," *SIAM J. Math. Anal.* 11:1047-1057
    (1980).
27. D. H. Martin, "Finite criteria for conditional definiteness of quadratic forms," *Linear Algebra Appl.* 39:9-21
    (1981).
28. D. H. Martin, M. J. D. Powell and D. H. Jacobson, "On a decomposition of conditionally positive-semidefinite
    matrices," *Linear Algebra Appl.* 39:51-59 (1981).
29. D. H. Martin and D. H. Jacobson, "Copositive matrices and definiteness of quadratic forms subject to homogeneous
    linear inequality constraints," *Linear Algebra Appl.* 35:227-258 (1981).
30. J. E. Maxfield and H. Minc, "On the matrix equation \(X'X=A\)," *Proc. Edinburgh Math. Soc.* 13:125-129 (1962/3).
31. T. S. Motzkin, "Copositive quadratic forms," Natl. Bur. of Standards Report 1818, 1952, pp. 11-22.
32. T. S. Motzkin, "Quadratic forms positive for nonnegative variables not all zero," *Notices Amer. Math. Soc.* 12:224
    (1965).
33. T. S. Motzkin and E. G. Strauss, "Maxima for graphs and a new proof of a theorem of Turan," *Canad. J. Math.*
    17:533-540 (1965).
34. T. S. Motzkin, "Signs of minors," in *Inequalities* (O. Shisha, Ed.), Academic, New York, 1967, pp. 225-240.
35. F. J. Pereira R., "On characterizations of copositive matrices," Ph.D. Thesis, Dept. of Operations Research, Stanford
    Univ., 1972.
36. A. C. Williams, "Boundedness relations for linear constraint sets," *Linear Algebra Appl.* 3:129-141 (1970).

[^1]: The printed paper reads \(\beta=b+d+c\). Direct differentiation of the displayed system gives
      \(\beta=b+d+e\), which is also the value required by the paper's following equivalence
      \(\beta\neq\gamma\iff b+d\neq c+f\); the printed \(c\) therefore appears to be a typographical error.
