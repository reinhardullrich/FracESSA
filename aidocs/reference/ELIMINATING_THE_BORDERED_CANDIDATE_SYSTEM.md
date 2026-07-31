# Eliminating the Border from the Candidate System

## Purpose

This note explains one mathematical simplification of the candidate equations.
It assumes only elementary matrix multiplication, transposes, systems of linear
equations, and the meaning of a singular matrix.

The subject is deliberately narrow. We only derive a smaller linear system for
the candidate probabilities. We do not discuss numerical algorithms, speed,
outside-support inequalities, or evolutionary stability.

## 1. The candidate equations on one support

Suppose a support contains $k$ strategies. Restrict the payoff matrix to these
strategies and call the resulting $k\times k$ matrix $A$.

We assume that $A$ is symmetric:

$$
A^T=A.
$$

We want to find a probability vector

$$
x=
\begin{bmatrix}
x_1\\
x_2\\
\vdots\\
x_k
\end{bmatrix}
$$

and a common payoff $u$.

Every strategy in the support must receive the same payoff. In matrix form,
this means

$$
Ax=u\mathbf 1,
$$

where

$$
\mathbf 1=
\begin{bmatrix}
1\\
1\\
\vdots\\
1
\end{bmatrix}.
$$

Because $x$ is a probability vector, its entries must add up to one:

$$
\mathbf 1^T x=1.
$$

For the moment, we are only studying these equalities. The additional
requirements $x_i>0$ and the inequalities for strategies outside the support
are separate questions.

## 2. The bordered system

The equal-payoff and normalization equations can be written as one system:

$$
\begin{bmatrix}
A&-\mathbf 1\\
\mathbf 1^T&0
\end{bmatrix}
\begin{bmatrix}
x\\
u
\end{bmatrix}
=
\begin{bmatrix}
0\\
1
\end{bmatrix}.
$$

The last row and last column are called the **border**. They were added to the
matrix to express the normalization equation and the unknown common payoff.

Multiplying the last equation by $-1$ gives the equivalent symmetric system

$$
K
\begin{bmatrix}
x\\
u
\end{bmatrix}
=
\begin{bmatrix}
0\\
-1
\end{bmatrix},
\qquad
K=
\begin{bmatrix}
A&-\mathbf 1\\
-\mathbf 1^T&0
\end{bmatrix}.
$$

This system has $k+1$ unknowns: the $k$ entries of $x$ and the payoff $u$.
We will now eliminate both $u$ and one entry of $x$. The remaining system will
have only $k-1$ unknowns.

## 3. Use normalization to eliminate one probability

Choose one strategy as a reference strategy. We call it strategy $m$. Relabeling
the strategies does not change the mathematics, so we may choose

$$
m=k.
$$

The normalization equation is

$$
x_1+x_2+\cdots+x_{k-1}+x_k=1.
$$

Therefore,

$$
x_k=1-\sum_{j=1}^{k-1}x_j.
$$

Introduce a shorter vector containing the first $k-1$ probabilities:

$$
y=
\begin{bmatrix}
y_1\\
y_2\\
\vdots\\
y_{k-1}
\end{bmatrix}
=
\begin{bmatrix}
x_1\\
x_2\\
\vdots\\
x_{k-1}
\end{bmatrix}.
$$

Then the complete probability vector is

$$
x=
\begin{bmatrix}
y_1\\
y_2\\
\vdots\\
y_{k-1}\\
1-y_1-y_2-\cdots-y_{k-1}
\end{bmatrix}.
$$

This substitution builds normalization directly into $x$. Every choice of $y$
produces a vector whose entries add up to one.

## 4. Write the substitution with basis vectors

Let $e_i$ denote the vector whose $i$-th entry is one and whose other entries
are zero. For example, in dimension three,

$$
e_1=
\begin{bmatrix}1\\0\\0\end{bmatrix},
\qquad
e_2=
\begin{bmatrix}0\\1\\0\end{bmatrix},
\qquad
e_3=
\begin{bmatrix}0\\0\\1\end{bmatrix}.
$$

For every $i<k$, define

$$
z_i=e_i-e_k.
$$

Each $z_i$ has one entry equal to $1$, one entry equal to $-1$, and all other
entries equal to zero. Consequently, its entries add up to zero:

$$
\mathbf 1^Tz_i=0.
$$

Put these vectors next to one another as the columns of a matrix:

$$
Z=
\begin{bmatrix}
z_1&z_2&\cdots&z_{k-1}
\end{bmatrix}.
$$

With reference strategy $k$, this matrix has the simple form

$$
Z=
\begin{bmatrix}
1&0&\cdots&0\\
0&1&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
0&0&\cdots&1\\
-1&-1&\cdots&-1
\end{bmatrix}.
$$

The first $k-1$ rows form an identity matrix. Therefore, the columns of $Z$
are linearly independent.

The probability vector can now be written compactly as

$$
x=e_k+Zy.
$$

This formula is exactly the same substitution as before. It also makes the
normalization especially easy to verify:

$$
\mathbf 1^Tx
=\mathbf 1^Te_k+\mathbf 1^TZy
=1+0
=1.
$$

The columns of $Z$ describe all changes to a vector that leave the sum of its
entries unchanged.

## 5. Eliminate the common payoff

The equal-payoff equation is

$$
Ax=u\mathbf 1.
$$

Written one row at a time, it says

$$
(Ax)_1=u,
\quad
(Ax)_2=u,
\quad\ldots\quad,
(Ax)_k=u.
$$

Subtract the equation for the reference strategy $k$ from each of the first
$k-1$ equations. The unknown $u$ cancels:

$$
(Ax)_i-(Ax)_k=0,
\qquad i=1,\ldots,k-1.
$$

Because $z_i=e_i-e_k$, this difference can be written as

$$
z_i^TAx=0.
$$

Putting all $k-1$ equations together gives

$$
Z^TAx=0.
$$

Now substitute $x=e_k+Zy$:

$$
Z^TA(e_k+Zy)=0.
$$

Distributing the matrix multiplication gives

$$
Z^TAe_k+Z^TAZy=0.
$$

Move the first term to the right-hand side:

$$
Z^TAZy=-Z^TAe_k.
$$

Define

$$
H=Z^TAZ
$$

and

$$
r=-Z^TAe_k.
$$

The complete bordered system has therefore been reduced to

$$
Hy=r.
$$

This is a $(k-1)\times(k-1)$ system with $k-1$ unknowns.

## 6. The entries of the reduced system

The abstract expression $H=Z^TAZ$ has a direct entry-by-entry formula. For
$i,j<k$,

$$
H_{ij}
=(e_i-e_k)^TA(e_j-e_k).
$$

Multiplying out the four terms gives

$$
H_{ij}
=A_{ij}-A_{ik}-A_{kj}+A_{kk}.
$$

The right-hand side has entries

$$
r_i
=-(e_i-e_k)^TAe_k
=A_{kk}-A_{ik}.
$$

Thus the reduced equations can be written without $Z$ as

$$
\sum_{j=1}^{k-1}
\left(
A_{ij}-A_{ik}-A_{kj}+A_{kk}
\right)y_j
=A_{kk}-A_{ik},
\qquad i=1,\ldots,k-1.
$$

## 7. Why the reduced matrix is symmetric

We assumed that $A$ is symmetric. Therefore,

$$
A^T=A.
$$

Take the transpose of $H$:

$$
H^T
=(Z^TAZ)^T
=Z^TA^TZ
=Z^TAZ
=H.
$$

Hence $H$ is symmetric.

The same fact is visible from the entry formula. Since $A_{ij}=A_{ji}$,

$$
H_{ij}=H_{ji}.
$$

## 8. Reconstruct the original unknowns

Suppose the reduced system

$$
Hy=r
$$

has been solved.

First reconstruct the probabilities:

$$
x_i=y_i,
\qquad i=1,\ldots,k-1,
$$

and

$$
x_k=1-\sum_{i=1}^{k-1}y_i.
$$

Then calculate the common payoff from the reference strategy:

$$
u=(Ax)_k.
$$

The reduced equations guarantee that every other strategy in the support has
the same payoff:

$$
(Ax)_i-(Ax)_k=0.
$$

Therefore,

$$
(Ax)_i=u
$$

for every $i=1,\ldots,k$.

## 9. Why the two systems are equivalent

We now prove that the reduction neither loses nor invents solutions.

### From the bordered system to the reduced system

Assume that $x$ and $u$ satisfy

$$
Ax=u\mathbf 1,
\qquad
\mathbf 1^Tx=1.
$$

The normalization equation lets us write exactly one vector $y$ such that

$$
x=e_k+Zy.
$$

Multiplying $Ax=u\mathbf 1$ by $Z^T$ gives

$$
Z^TAx=uZ^T\mathbf 1=0.
$$

Substituting $x=e_k+Zy$ gives

$$
Hy=r.
$$

Thus every solution of the bordered system produces a solution of the reduced
system.

### From the reduced system to the bordered system

Now assume that $y$ satisfies

$$
Hy=r.
$$

Set

$$
x=e_k+Zy.
$$

This automatically gives

$$
\mathbf 1^Tx=1.
$$

The reduced equations give

$$
Z^TAx=0.
$$

Each entry of this vector equation says

$$
(Ax)_i-(Ax)_k=0.
$$

Define

$$
u=(Ax)_k.
$$

Then every entry of $Ax$ equals $u$, so

$$
Ax=u\mathbf 1.
$$

Thus every solution of the reduced system produces a solution of the bordered
system.

The two systems therefore describe exactly the same probability vectors and
common payoffs.

## 10. Singularity is preserved

For the candidate equations, it is also important to know whether the solution
is unique. A square matrix is singular exactly when its homogeneous system has
a nonzero solution.

The homogeneous system belonging to the symmetric bordered matrix $K$ is

$$
K
\begin{bmatrix}
d\\
\alpha
\end{bmatrix}
=0.
$$

Written as separate equations, this is

$$
Ad-\alpha\mathbf 1=0
$$

and

$$
\mathbf 1^Td=0.
$$

The second equation says that $d$ has entries whose sum is zero. Therefore,
there is a vector $z$ such that

$$
d=Zz.
$$

To see this directly, take $z_1=d_1,\ldots,z_{k-1}=d_{k-1}$. Since the entries
of $d$ add up to zero, its last entry must be

$$
d_k=-d_1-d_2-\cdots-d_{k-1}.
$$

This is exactly the vector produced by $Zz$.

### If the bordered matrix is singular, then $H$ is singular

Suppose the bordered homogeneous system has a nonzero solution. Multiplying

$$
Ad=\alpha\mathbf 1
$$

by $Z^T$ gives

$$
Z^TAd=\alpha Z^T\mathbf 1=0.
$$

Using $d=Zz$, we obtain

$$
Z^TAZz=Hz=0.
$$

The vector $d$ cannot be zero, because $d=0$ in
$Ad-\alpha\mathbf 1=0$ would also force $\alpha=0$. The columns of $Z$ are
linearly independent, so $d=Zz\neq0$ also implies $z\neq0$. Therefore, $H$
has a nonzero null vector and is singular.

### If $H$ is singular, then the bordered matrix is singular

Now suppose that $H$ is singular. Then there is a nonzero vector $z$ such that

$$
Hz=0.
$$

Set

$$
d=Zz.
$$

Because $z\neq0$ and the columns of $Z$ are linearly independent, $d\neq0$.

Then

$$
Z^TAd=0.
$$

This equation says that all entries of $Ad$ are equal. Call their common value
$\alpha$. Then

$$
Ad=\alpha\mathbf 1.
$$

Also,

$$
\mathbf 1^Td=\mathbf 1^TZz=0.
$$

Consequently,

$$
K
\begin{bmatrix}
d\\
\alpha
\end{bmatrix}
=0
$$

has a nonzero solution, and $K$ is singular.

We have proved

$$
K\text{ is nonsingular}
\quad\Longleftrightarrow\quad
H\text{ is nonsingular}.
$$

Therefore, the smaller system makes exactly the same unique-solution decision
as the original bordered system.

## 11. A complete three-strategy example

Consider the symmetric matrix

$$
A=
\begin{bmatrix}
3&0&0\\
0&2&0\\
0&0&1
\end{bmatrix}.
$$

Choose strategy $3$ as the reference strategy. Then

$$
Z=
\begin{bmatrix}
1&0\\
0&1\\
-1&-1
\end{bmatrix}.
$$

The reduced matrix is

$$
H=Z^TAZ
=
\begin{bmatrix}
4&1\\
1&3
\end{bmatrix}.
$$

The right-hand side is

$$
r=-Z^TAe_3
=
\begin{bmatrix}
1\\
1
\end{bmatrix}.
$$

We solve

$$
\begin{bmatrix}
4&1\\
1&3
\end{bmatrix}
\begin{bmatrix}
y_1\\
y_2
\end{bmatrix}
=
\begin{bmatrix}
1\\
1
\end{bmatrix}.
$$

The solution is

$$
y_1=\frac{2}{11},
\qquad
y_2=\frac{3}{11}.
$$

Reconstructing the third probability gives

$$
x_3
=1-y_1-y_2
=1-\frac{2}{11}-\frac{3}{11}
=\frac{6}{11}.
$$

Thus

$$
x=
\begin{bmatrix}
\frac{2}{11}\\
\frac{3}{11}\\
\frac{6}{11}
\end{bmatrix}.
$$

Finally,

$$
Ax=
\begin{bmatrix}
3\cdot\frac{2}{11}\\
2\cdot\frac{3}{11}\\
1\cdot\frac{6}{11}
\end{bmatrix}
=
\begin{bmatrix}
\frac{6}{11}\\
\frac{6}{11}\\
\frac{6}{11}
\end{bmatrix}.
$$

Every strategy has the same payoff, so

$$
u=\frac{6}{11}.
$$

The original four unknowns

$$
x_1,x_2,x_3,u
$$

were found by solving a symmetric $2\times2$ system.

## 12. Why we do not simply invert $A$

If $A$ were known to be nonsingular, the equation

$$
Ax=u\mathbf 1
$$

could be written as

$$
x=uA^{-1}\mathbf 1.
$$

Normalization would then determine $u$. However, this argument requires
$A^{-1}$ to exist.

That requirement is too strong. The matrix $A$ can be singular even when the
bordered candidate system is nonsingular. In the one-strategy example

$$
A=\begin{bmatrix}0\end{bmatrix},
$$

the matrix $A$ is singular, but

$$
K=
\begin{bmatrix}
0&-1\\
-1&0
\end{bmatrix}
$$

has determinant $-1$ and is nonsingular.

The reduction through $Z$ and $H=Z^TAZ$ does not assume that $A$ is invertible.
It is therefore equivalent to the bordered system in every case.

## 13. The one-strategy case

When $k=1$, there are no free probability variables. Normalization immediately
gives

$$
x_1=1,
$$

and the common payoff is

$$
u=A_{11}.
$$

The reduced matrix has size $0\times0$, so no linear system needs to be solved.
This is the natural endpoint of the same reduction.

## Summary

Starting from

$$
Ax=u\mathbf 1,
\qquad
\mathbf 1^Tx=1,
$$

choose one reference strategy and write

$$
x=e_k+Zy.
$$

The normalization equation is then automatic, and subtracting the reference
payoff eliminates $u$. The remaining system is

$$
Hy=r,
\qquad
H=Z^TAZ,
\qquad
r=-Z^TAe_k.
$$

It has size $(k-1)\times(k-1)$, its matrix $H$ is symmetric, and

$$
K\text{ is nonsingular}
\quad\Longleftrightarrow\quad
H\text{ is nonsingular}.
$$

After solving for $y$, the original probabilities and common payoff are
reconstructed exactly.
