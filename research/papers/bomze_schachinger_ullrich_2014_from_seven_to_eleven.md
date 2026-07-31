# From Seven to Eleven: Completely Positive Matrices with High Cp-Rank

**Authors:** Immanuel M. Bomze, Werner Schachinger, Reinhard Ullrich  
**Affiliation:** University of Vienna, Department of Statistics and Operations Research, 1090 Vienna, Austria  
**Journal:** *Linear Algebra and its Applications*, Vol. 459, 2014, pp. 208-221  
**Source version:** Preprint submitted July 16, 2014  
**Source PDF:** [bomze_schachinger_ullrich_2014_from_seven_to_eleven.pdf](bomze_schachinger_ullrich_2014_from_seven_to_eleven.pdf)  
**Exact LaTeX source:** [bomze_schachinger_ullrich_2014_from_seven_to_eleven.tex](bomze_schachinger_ullrich_2014_from_seven_to_eleven.tex)

> Complete text transcription based on the exact LaTeX source that generated the stored PDF. The formulas, theorem and
> equation numbers, tables, and references were checked against the PDF. Journal running headers and page numbers are
> omitted.

## Abstract

We study $n\times n$ completely positive matrices ${\mathsf M}$ on the boundary of the completely positive cone, namely
those orthogonal to a copositive matrix ${\mathsf S}$ which generates a quadratic form with finitely many zeroes in the
standard simplex. Constructing particular instances of ${\mathsf S}$, we are able to construct counterexamples to the
famous Drew-Johnson-Loewy conjecture (1994) for matrices of order seven through eleven.

**Keywords:** Copositive optimization; completely positive matrices; cp-rank; nonnegative factorization; circular
symmetry.  
**2010 Mathematics Subject Classification:** 15B48; 90C25; 15A23.

## 1. Introduction

In this article we consider completely positive matrices ${\mathsf M}$ and their cp-rank. An $n\times n$ matrix ${\mathsf M}$ is said to be *completely positive* if there exists a nonnegative (not necessarily square) matrix ${\mathsf V}$ such that ${\mathsf M}={\mathsf V}{\mathsf V}^\top$. Typically, a completely positive matrix ${\mathsf M}$ may have many such factorizations, and the *cp-rank* of ${\mathsf M}$, $\mbox{cpr }{\mathsf M}$, is the minimum number of columns in such a nonnegative factor ${\mathsf V}$ (for completeness, we define $\mbox{cpr }{\mathsf M}= 0$ if ${\mathsf M}$ is a square zero matrix and $\mbox{cpr }{\mathsf M}= \infty$ if ${\mathsf M}$ is not completely positive). Completely positive matrices form a cone dual to the cone of *copositive matrices*. An $n\times n$ matrix ${\mathsf S}$ is said to be copositive if $\mathbf x^\top{\mathsf S} \mathbf x\ge 0$ for every nonnegative vector $\mathbf x\in {\mathbb R}^n_+$. Both cones are central in the rapidly evolving field of *copositive optimization* which links discrete and continuous optimization, and has numerous real-world applications. For recent surveys and structured bibliographies, we refer to [5; 6; 8; 12], and for a fundamental text book to [2].

Determining the maximum possible cp-rank of $n\times n$ completely positive matrices, $$p_n := \max\left\{\mbox{cpr }{\mathsf M}: {\mathsf M}\mbox{ is a completely positive }n\times n\mbox{ matrix}\right \} ,$$ is still an open problem for general $n$. It is known [2, Theorem 3.3] that $p_n=n$ if $n\le 4$, whereas this equality does no longer hold for $n\ge 5$. Let $d_n := \left\lfloor \frac{n^2}4\right\rfloor$ and $s_n := {n+1\choose 2}-4$. For $n\ge 5$, it is known that [16] $$
d_n\le p_n \le s_n \, ,
\tag{1}
$$ and that $d_n=p_n$ in case $n=5$ [17]. It is still unknown whether $d_6=p_6$ although the bracket (1) was reduced in the recent paper [16] where also the upper bound $p_n\le s_n$ was established for the first time.

The famous Drew-Johnson-Loewy (DJL) conjecture [11] is by now twenty years old. It states that $d_n=p_n$ is true for all $n\ge 5$, and some evidence in support of the DJL conjecture is found in [1; 10; 11; 15], see also [2, Section 3.3]. However, we will show in this paper that the DJL conjecture does not hold for $n\in \left\{7,8,9,10,11\right \}$ by constructing examples which establish $p_n > d_n$.

The paper is organized as follows: In Section 2 we look at copositive matrices ${\mathsf S}$ which allow for finitely many (but many) zeroes $\mathbf q_i$ of the quadratic form $\mathbf x^\top{\mathsf S} \mathbf x$ over the standard simplex. Such matrices ${\mathsf S}$ lie on the boundary of the copositive cone, and elementary conic duality therefore tells us that there are nontrivial completely positive matrices ${\mathsf M}$ such that ${\mathsf M}\perp {\mathsf S}$ in the Frobenius inner product sense, and we will study the cp-rank of these ${\mathsf M}$. Section 3 deals with a particular construction of above mentioned copositive matrices ${\mathsf S}$ (they will be cyclically symmetric) in a way that many $\mathbf q_i$ can coexist, and in Section 4 we present the second main result -- counterexamples to the DJL conjecture for $7\le n\le 11$. Let us mention here that such a counterexample for $n=7$ with cp-rank 14 was announced in 2002, according to [2, p.177]. The matrix there (which never got public) should have rank 5; by contrast, our matrix ${\mathsf M}$ in Example 1 will have full rank $7$, but also $\mbox{cpr }{\mathsf M}=14$ by mere coincidence.

Some notation and terminology: we abbreviate $[r\!:\!s]=\left\{r,r+1,\ldots ,s\right \}$ for integers $r\le s$, and by $| S|$ the number of elements of a finite set $S$. For a function $f:T\to {\mathbb R}$ we abbreviate $$\textup{Argmin}\left\{f(t):t\in T\right \} := \left\{\bar t\in T : f(\bar t) \le f(t) \mbox{ for all }t\in T\right \} \, .$$ The nonnegative orthant is denoted by ${\mathbb R}^n_+$. For a vector $\mathbf x\in {\mathbb R}^n_+$, the index set $$I(\mathbf x)= \left\{i\in [1\!:\!n] : x_i > 0\right \}$$ is the *support* of $\mathbf x$. Let $\mathbf e_i$ be the $i$th column vector of the $n\times n$ identity matrix ${\mathsf I}_n$ and $\mathbf e=\sum_{i=1}^n \mathbf e_i$. The zero vector and the zero matrix (of appropriate sizes) are denoted by $\mathbf o$ and ${\mathsf O}$, respectively, and $\Delta = \{\mathbf x\in {\mathbb R}^n_+ : \mathbf e^\top\mathbf x= 1 \}$ stands for the standard simplex. The vector space of real symmetric $n\times n$ matrices is denoted by ${{\cal S}^n}$, and the Frobenius inner product of two matrices $\left\{{\mathsf A},{\mathsf B}\right \} \subset {{\cal S}^n}$ by $\langle {\mathsf A},{\mathsf B}\rangle := \mbox{trace }({\mathsf A}{\mathsf B})$. For an $n\times p$ matrix ${\mathsf V}=[\mathbf v_1, \ldots ,\mathbf v_p]$, the relation ${\mathsf M}={\mathsf V}{\mathsf V}^\top$ is equivalent to ${\mathsf M}=\sum_{i=1}^p \mathbf v_i\mathbf v_i^\top$. We will refer to this sum as a "cp decomposition" of ${\mathsf M}$, if ${\mathsf V}$ has no negative entries. Given a square matrix ${\mathsf S}$, we will, by slight abuse of language, use the phrase "zero(es) of ${\mathsf S}$" as an abbreviation of "zero(es) of the quadratic form $\mathbf x^\top{\mathsf S} \mathbf x$ over $\mathbf x\in \Delta$"; this terminology differs slightly from that in [14].

By ${{\cal C}^{n*}}$ we denote the cone of completely positive matrices, $${{\cal C}^{n*}}= \mbox{conv\,}\{\mathbf x\mathbf x^\top:  \mathbf x\in {\mathbb R}^n_+ \}\, .$$ Both, ${{\cal C}^{n*}}$ and its dual, the cone of copositive matrices $${\cal C}^{n} = \left\{{\mathsf S}\in {{\cal S}^n}: \mathbf x^\top{\mathsf S}\mathbf x\ge0 \mbox{ for all }\mathbf x\in {\mathbb R}^n_+\right \} \, ,$$ are pointed closed convex cones with nonempty interior. The copositive cone ${\cal C}^{n}$ and, in particular, its extremal rays, are important as any matrix on the boundary $\partial {{\cal C}^{n*}}$ of ${{\cal C}^{n*}}$ is orthogonal to an extremal ray of ${\cal C}^{n}$. So, studies of the extremal rays of ${\cal C}^{n}$ like in [9; 13; 14] lead to conclusions on all matrices on $\partial {{\cal C}^{n*}}$, which allow for inference on *upper* bounds on $p_n$. This was an essential ingredient of the arguments in [16; 17]. Here we employ a somewhat reverse approach: we start from (appropriate) matrices ${\mathsf S}\in \partial {\cal C}^{n}$ and construct ${\mathsf M}\in \partial {{\cal C}^{n*}}$ where we can calculate the cp-rank $\mbox{cpr }{\mathsf M}$, improving upon *lower* bounds on $p_n$. Eventually, this will lead to examples refuting the DJL conjecture.

## 2. Iterative Reduction of the Cp-Rank

Consider a copositive matrix ${\mathsf S}\in \partial {\cal C}^{n}$ and assume that $\left\{\mathbf q_1,\ldots ,\mathbf q_m\right \}$ are all the zeroes of ${\mathsf S}$. Since ${\mathsf S}\in \partial {\cal C}^{n}$, there is a matrix ${\mathsf M}\in {{\cal C}^{n*}}\setminus \left\{{\mathsf O}\right \}$ such that $\langle {{\mathsf M}} , {{\mathsf S}} \rangle=0$, e.g., any matrix of the form $${\mathsf M}=\sum_{i=1}^m y_i \mathbf q_i\mathbf q_i^\top$$ for some $\mathbf y\in {\mathbb R}^m_+\setminus \left\{\mathbf o\right \}$. The next result shows the converse of this statement, so that the set of possible cp decompositions of matrices orthogonal to ${\mathsf S}$ is quite restricted:

### Lemma 2.1

*Let $Q=\left\{\mathbf x\in \Delta : \mathbf x^\top{\mathsf S} \mathbf x= 0\right \}$ be all the zeroes of ${\mathsf S}\in \partial {\cal C}^{n}$. Then any matrix ${\mathsf M}\in {{\cal C}^{n*}}$ orthogonal to ${\mathsf S}$ must be of the form $$
{\mathsf M}=\sum_{j=1}^m y_j \mathbf q_j\mathbf q_j^\top\quad \mbox{with }\left\{\mathbf q_1,\ldots , \mathbf q_m\right \} \subseteq Q
\tag{2}
$$ for some $\mathbf y\in {\mathbb R}^m_+$.*

*Proof.*
Let ${\mathsf M}$ have the cp decomposition ${\mathsf M}=\sum\limits_{i=1}^m \mathbf v_i\mathbf v_i^\top$ with $\mathbf v_i\in {\mathbb R}^n_+\setminus\left\{\mathbf o\right \}$ for all $i\in [1\!:\! m]$. Then ${\mathsf M}\perp {\mathsf S}$ implies $$0=\langle {{\mathsf M}} , {{\mathsf S}} \rangle= \sum_{i=1}^m \mathbf v_i^\top{\mathsf S} \mathbf v_i\, ,$$ and as ${\mathsf S}$ is copositive, every term in above sum must be zero. So all $\mathbf q_i := \frac 1{\mathbf e^\top\mathbf v_i }\, \mathbf v_i\in  Q$, and the result (2) follows with $y_i:= (\mathbf e^\top\mathbf v_i)^2$. $\Box$

Although we have restricted the possible cp decompositions by above observation, there still could be infinitely many, but they can be obtained in a linear way. To be more precise, suppose that $Q=\left\{\mathbf q_1,\ldots, \mathbf q_m\right \}$, fix any $\mathbf y\in {\mathbb R}^m_+$ such that (2) holds, and consider $$
X_\mathbf y:= \left\{\mathbf x\in {\mathbb R}^m_+: \sum_{i=1}^m x_i\mathbf q_i\mathbf q_i^\top= \sum_{j=1}^m y_j\mathbf q_j\mathbf q_j^\top\right \} \, .
\tag{3}
$$ A particular case is obtained if $X_\mathbf y= \left\{\mathbf y\right \}$, because then $\mbox{cpr }{\mathsf M}= | I(\mathbf y)|$ is immediate from Lemma 2.1. However, this may not always be the case, but the next theorem will show how to fix some variables $x_k$ of points $\mathbf x\in X_\mathbf y$ to $y_k$, with some consequences for the construction of matrices of high cp-rank. To apply that theorem in more general situations, we need some further notation. First define for $Q\subseteq\Delta$ the set $\overline Q:=\left\{\mathbf q\mathbf q^\top: \mathbf q\in Q\right \} \subset {{\cal S}^n}$ and by $\mbox{cone }\! \overline Q := {\mathbb R}_+\mbox{conv\,}\overline Q$ the *convex conic hull* of $\overline Q$; moreover, for finite $P\subset \Delta$, we denote by $$\operatorname{cone}^{\circ} \!\! \overline P := \left\{\sum_{\mathbf f\in P} y_\mathbf f\mathbf f{}\mathbf f^\top: y_\mathbf f>0 \mbox{ for all }\mathbf f\in P\right \} \, .$$ Finally, we abbreviate the set of completely positive matrices whose cp decompositions can only use multiples of vectors from $Q$ by $$\mathcal{E}(Q) := \left\{{\mathsf M}\in \mbox{cone }\! \overline Q : \mbox{if }\,{\mathsf M}\in \operatorname{cone}^{\circ} \!\! \overline P \textup{ for finite }P\subseteq\Delta\, , \textup{ then } P\subseteq Q\right \} \, .$$ So Lemma 2.1 would read: If $Q$ is the set of all zeroes of ${\mathsf S}\in{\cal C}^{n}$, then any ${\mathsf M}\in{{\cal C}^{n*}}$ with $\langle{\mathsf M},{\mathsf S}\rangle=0$ satisfies ${\mathsf M}\in\mathcal{E}(Q)$.

### Theorem 2.1

*For a finite subset $Q\subset\Delta$ consider ${\mathsf M}\!=\!\sum\limits_{\mathbf f\in Q} y_\mathbf f\mathbf f{}\mathbf f^\top$ with $\mathbf y\!\in\!\mathbb{R}^{|Q|}_+$ and assume ${\mathsf M}\in\mathcal{E}(Q)$. Suppose that there is $\mathbf q\in Q$ such that for two different indices $r,s$, we have $$
\left\{r,s\right \} \subseteq I(\mathbf q) \quad\mbox{but}\quad \left\{r,s\right \} \not\subseteq I(\mathbf q')\mbox{ for all }\mathbf q'\in Q\setminus\{\mathbf q\} \, .
\tag{4}
$$ Then*

1.  *$x_\mathbf q= \frac {\mathbf e_r^\top{\mathsf M}\mathbf e_s}{(\mathbf e_r^\top\mathbf q)(\mathbf e_s^\top\mathbf q)}= y_\mathbf q$ holds for all $\mathbf x\in X_\mathbf y$,*

2.  *$\widehat{\mathsf M}{}: = {\mathsf M}- y_\mathbf q\mathbf q\mathbf q^\top\in\mathcal{E}\left(Q\setminus\{\mathbf q\}\right)$,*

3.  *$\mbox{cpr }{\mathsf M}= \mbox{sgn\,}(y_\mathbf q) + \mbox{cpr }\widehat{\mathsf M}$.*

*Proof.*
Condition (4) implies $(\mathbf e_r^\top\mathbf q)(\mathbf e_s^\top\mathbf q)>0$ and further that $$x_k (\mathbf e_r^\top\mathbf q)(\mathbf e_s^\top\mathbf q)=\mathbf e_r^\top{\mathsf M}\mathbf e_s\quad \mbox{for all }\mathbf x\in X_\mathbf y\,.$$ Hence $x_k= \frac {\mathbf e_r^\top{\mathsf M}\mathbf e_s}{(\mathbf e_r^\top\mathbf q)(\mathbf e_s^\top\mathbf q)}$ is fixed, which proves (a). Now define $\widehat y_\mathbf q=0$ and $\widehat y_{\mathbf q'}=y_{\mathbf q'}$ for $\mathbf q'\in Q\setminus\{\mathbf q\}$, and observe $\widehat {\mathsf M}=\sum_{\mathbf f\in Q}\widehat y_\mathbf f\mathbf f{}\mathbf f^\top\in\mathcal{E}(Q)$. Assertion (a), applied to $\widehat\mathbf y$, tells us $\widehat x_\mathbf q=\widehat y_\mathbf q=0$ for all $\widehat\mathbf x\in X_{\widehat y}$, therefore (b) holds. By (a), any minimal cp decomposition of ${\mathsf M}$ is of the form ${\mathsf M}=y_\mathbf q\mathbf q\mathbf q^\top+\widehat{\mathsf M}$, which implies (c). $\Box$

If the hypotheses of Theorem 2.1 including condition (4) are satisfied for ${\mathsf M},Q,\mathbf q$, then by (b) of that theorem we know that $\widehat{\mathsf M}\in \mathcal{E}(\widehat Q)$ for $\widehat Q:=Q\setminus\{\mathbf q\}$. Now, if we want to apply Theorem 2.1 iteratively, then we may replace $Q$ with $\widehat Q$ so that condition (4) may be satisfied more easily for some $\widehat \mathbf q\in \widehat Q$.

So if we arrange the supports of (many) $\mathbf q$'s such that condition (4), or a similar one, continues to hold during the iterations, we can construct ${\mathsf M}$ with high $\mbox{cpr }{\mathsf M}$, as long as $y_\mathbf q>0$ continues to hold, too. This will be done in the next section.

## 3. Zeroes of Cyclically Symmetric Matrices

We will employ symmetry transformations of the coordinates given by cyclic permutation, denoting by $a\oplus b$ and $a\ominus b$ the result of addition and subtraction modulo $n$. To keep in line with previous and standard notation, we consider the remainders $[1\!:\! n]$ instead of $[0\!:\! n-1]$, e.g. $1\oplus (n-1) = n$. To be more precise, let ${\mathsf P}_i$ be the square $n\times n$ permutation matrix which effects ${\mathsf P}_i\mathbf x= [x_{i\oplus j}]_{j\in [1:n]}$ for all $\mathbf x\in {\mathbb R}^n$ (for example, if $n=3$ then ${\mathsf P}_2\mathbf x= [x_3,x_1,x_2]^\top$). Obviously ${\mathsf P}_i=({\mathsf P}_1)^i$ for all integers $i$ (recall ${\mathsf P}^{-3}$ is the inverse matrix of ${\mathsf P}{\mathsf P}{\mathsf P}$), ${\mathsf P}_i^\top= {\mathsf P}_{n-i}= {\mathsf P}_i^{-1}$ and ${\mathsf P}_n= {\mathsf I}_n$. A *circulant matrix* ${\mathsf S}={\mathsf C}({\mathbf a})$ based on a vector ${\mathbf a}\in {\mathbb R}^n$ (as its last column rather than the first) is given by $${\mathsf S}=[{\mathsf P}_ {n-1}{\mathbf a},{\mathsf P}_ {n-2}{\mathbf a}, \ldots, {\mathsf P}_ 1{\mathbf a}, {\mathbf a}]\, .$$ If ${\mathsf S}={\mathsf C}({\mathbf a})\in {{\cal S}^n}$, i.e., if ${\mathsf C}({\mathbf a})$ is symmetric, it is called *cyclically symmetric*.

### Lemma 3.1

*Any circulant matrix ${\mathsf S}= {\mathsf C}({\mathbf a})$ satisfies ${\mathsf P}_i^\top{\mathsf S}{\mathsf P}_i={\mathsf S}$ for all $i\in [1\!:\! n]$. Furthermore, if $$
a_i=a_{n-i}\quad\mbox{for all }i\in[1\!:\! n-1]\, ,
\tag{5}
$$ then ${\mathsf S}={\mathsf C}({\mathbf a})$ is cyclically symmetric.*

*Proof.*
The first relation is evident. To show the remaining assertion, assume (5) and let $\mathbf e_j^\top{\mathsf S}\mathbf e_i=\mathbf e_j^\top{\mathsf P}_ {n-i}{\mathbf a}= a_k$ with $k\oplus i = j$ while $\mathbf e_i^\top{\mathsf S}\mathbf e_j= a_{ \ell }$ with $\ell \oplus j=i$. Thus $i\oplus j = k\oplus \ell \oplus i\oplus j$ and $\left\{k,\ell\right \} \subseteq[1\!:\! n]$, so we get $k+\ell\in \left\{n,2n\right \}$ and therefore $a_k=a_\ell$. Hence ${\mathsf C}({\mathbf a})\in {{\cal S}^n}$. $\Box$

Copositive cyclically symmetric matrices ${\mathsf S}={\mathsf C}({\mathbf a}) \in \partial{\cal C}^{n}$ can have many zeroes (which then are global minimizers of the quadratic form $\mathbf x^\top{\mathsf S} \mathbf x$ over $\Delta$; for local minimizers this has already been observed earlier, see [7] and references therein). To facilitate the argument, let us denote by ${\mathsf R}\in {{\cal S}^n}$ the *reflection matrix* which transforms every $\mathbf x\in {\mathbb R}^n$ into its *mirror image* ${\mathsf R}\mathbf x:=[x_{n+1-i}]_{i\in [1:n]}$. Note that ${\mathsf R}^\top={\mathsf R}\in {{\cal S}^n}$.

In the sequel, it will be convenient to denote, for any $\mathbf q\in\mathbb{R}^n_+$, the set $D_\mathbf q$ of differences and the set $U_\mathbf q$ of unique differences of the elements of $I(\mathbf q)$:\
$D_\mathbf q\! :=\!\{d\!\in\![1\!:\!n\!-\!1]:\, d\!=\!r\!\ominus\! s \textup{ has at least one solution with \!}\left\{r\!,\!s\right \} \! \subseteq\!I(\mathbf q)\}\,,$\
$U_\mathbf q\! :=\!\{d\!\in\![1\!:\!n\!-\!1]:\, d\!=\!r\!\ominus\! s \textup{ has exactly one solution with \!}\left\{r\!,\!s\right \} \! \subseteq I(\mathbf q)\}\,.$

### Lemma 3.2

*Let ${\mathsf S}={\mathsf C}({\mathbf a})\in {{\cal S}^n}$ be a cyclically symmetric matrix.*

1.  *We have ${\mathsf R}{\mathsf S}{\mathsf R}= {\mathsf S}$. Further, fixing $\mathbf q\in\mathbb{R}^n_+$, for any shift $\mathbf q'={\mathsf P}_i\mathbf q$, and for its mirror image $\mathbf q''={\mathsf R}\mathbf q$, we have $$
(\mathbf q')^\top{\mathsf S}\mathbf q'  = (\mathbf q'')^\top{\mathsf S}\mathbf q''=\mathbf q^\top{\mathsf S} \mathbf q\, .
\tag{6}
$$*

2.  *For any zero $\mathbf q$ of ${\mathsf S}$ there are actually up to $2n$ zeroes: the shifts ${\mathsf P}_i\mathbf q$ for $i\in[1:n]$ and their mirror images, if they are all different.*

3.  *The supports of zeroes are shifted cyclically, $I({\mathsf P}_i\mathbf q) =  \left\{j\ominus i : j\in I(\mathbf q)\right \}$. However, the relative differences within the support of course remain: if $\left\{r,s\right \} \subseteq I(\mathbf q)$, then $r\ominus s = r'\ominus s'$ if $r'=r\oplus i$ and $s'=s\oplus i$.*

4.  *For any $\mathbf q\in\mathbb{R}^n_+$ the sets $D_\mathbf q$ and $U_\mathbf q$ are invariant under shifts and reflection: with $\mathbf q',\,\mathbf q''$ as in (a), we have $D_\mathbf q\!=\!D_{\mathbf q'}\!=\!D_{\mathbf q''}$ and $U_\mathbf q\!=\!U_{\mathbf q'}\!=\!U_{\mathbf q''}$.\
    Moreover $d=\frac n2\not\in U_\mathbf q$ for even $n$, since then $r\!\ominus\! s\!=\!d$ implies $s\!\ominus\! r\!=\!d$.*

*Proof.*
The relation ${\mathsf R}{\mathsf S}{\mathsf R}= {\mathsf S}$ can be checked in a straightforward manner while the equations in (6) follow from $$(\mathbf q')^\top{\mathsf S}\mathbf q'  =\mathbf q^\top{\mathsf P}_i^\top{\mathsf S}{\mathsf P}_i\mathbf q=\mathbf q^\top{\mathsf S} \mathbf q$$ and from $$(\mathbf q'')^\top{\mathsf S}\mathbf q''   =\mathbf q^\top{\mathsf R}^\top{\mathsf S}{\mathsf R}\mathbf q=  \mathbf q^\top{\mathsf R}{\mathsf S}{\mathsf R}\mathbf q=  \mathbf q^\top{\mathsf S} \mathbf q\, .$$ The assertions about the supports are evident. $\Box$

### Theorem 3.1

*For a finite subset $Q\subset\Delta$ consider ${\mathsf M}\!=\!\sum\limits_{\mathbf f\in Q} y_\mathbf f\mathbf f{}\mathbf f^\top$ with $\mathbf y\!\in\!\mathbb{R}^{|Q|}_+$ and assume ${\mathsf M}\in\mathcal{E}(Q)$. Fix $\mathbf q\in Q$, define $Q_1:=\{{\mathsf P}_i\mathbf q:i\in[1\!:\! n]\}$ and define $Q_2:=Q\setminus Q_1$. Assume $Q_1\subseteq Q$, and that there is $d\in[1\!:\! n-1]$ such that $d\in U_{\mathbf q}\setminus\bigcup\limits_{\mathbf q'\in Q_2}D_{\mathbf q'}$. Then*

1.  *$x_{\mathbf f}  = y_{\mathbf f}$ holds for all $\mathbf x\in X_\mathbf y$ and $\mathbf f\in Q_1$,*

2.  *$\widehat{\mathsf M}{}: = {\mathsf M}- \sum\limits_{\mathbf f\in Q_1}y_\mathbf f\mathbf f{}\mathbf f^\top
         \in\mathcal{E} \left(Q_2\right)$,*

3.  *$\mbox{cpr }{\mathsf M}= \sum\limits_{\mathbf f\in Q_1}\mbox{sgn\,}(y_\mathbf f) + \mbox{cpr }\widehat{\mathsf M}$.*

*Proof.*
Let $\{r,s\}\subseteq I(\mathbf q)$ satisfy $d=r\ominus s$. By the assumptions it is clear that $\left\{r,s\right \} \subseteq I(\mathbf q')$ can never hold for any $\mathbf q' \in Q_2$. So consider instead $\mathbf f={\mathsf P}_i\mathbf q$ for $i\in [1\!:\! n-1]$. We argue by contradiction: if $\left\{r,s\right \} \subseteq I( \mathbf f)$, then $\left\{r\oplus i,s\oplus i\right \} \subseteq I(\mathbf q)$ but differs from the pair $\left\{r,s\right \}$ (note that $r=s\oplus i$ and simultaneously $s=r\oplus i$ is impossible since $d\neq\frac n2$). Obviously the difference would be the same, namely $d$, which by assumption is absurd. So condition (4) holds for $Q$ and $\mathbf q$. Since $U_{\mathbf f}=U_{\mathbf q}$ for all $\mathbf f\in Q_1$, by Lemma 3.2 (d), we similarly obtain that condition (4) holds for $Q$ and $\mathbf f\in Q_1\setminus\{\mathbf q\}$. Finally we obtain (a), (b) and (c) by iterating the reduction step of Theorem 2.1 in total $|Q_1|$ times. $\Box$

### Corollary 3.1

*Let all hypotheses of Theorem 3.1 be satisfied with $Q_2=\emptyset$. Let ${\mathsf M}:=\sum\limits_{\mathbf f\in Q}y_\mathbf f\,\mathbf f\,\!\mathbf f^\top$. If $y_\mathbf f> 0$ for all $\mathbf f\in Q$, then the minimal cp decomposition of ${\mathsf M}$ is unique and $\mbox{cpr }{\mathsf M}=| Q|$.*

The next two results deal with instances where there is more than one minimal cp decomposition of a similarly constructed matrix:

### Lemma 3.3

*Consider $\mathbf q\in\mathbb{R}^n_+$ such that $Q:=\{{\mathsf P}_i\mathbf q:i\in[1\!:\! n]\}$ satisfies $|Q|=n$ and ${\mathsf R}\mathbf q\notin Q$. Suppose there are $d_1,d_2\in\ U_\mathbf q$ with $d_1=r\ominus s$ and $d_2=\rho\ominus \sigma$, such that $\rho+\sigma-r-s$ and $n$ are coprime. We consider the following subset of ${{\cal S}^n}$: $$\mathcal{F}:=\left\{\mathbf f\,\!\mathbf f^\top:\mathbf f\in Q\right \} \cup\{{\mathsf R}\mathbf f({\mathsf R}\mathbf f)^\top:\mathbf f\in Q\}\, .$$ Then every $(2n-1)$-element subset of $\mathcal{F}$ is linearly independent, moreover $\mathcal{F}$ itself (as a subset of the vector space ${{\cal S}^n}$) has rank $2n-1$.*

*Proof.*
We first observe that our assumptions on $Q$ imply $|\mathcal{F}|=2n$. Moreover, $U_{{\mathsf R}\mathbf q}=U_\mathbf q$. We now claim that $$
\sum_{\mathbf f\in Q }{\mathsf R}\mathbf f({\mathsf R}\mathbf f)^\top={\mathsf R}\left(\sum_{\mathbf f\in Q }\mathbf f\,\!\mathbf f^\top\right) {\mathsf R}
=\sum_{\mathbf f\in Q }\mathbf f\,\! \mathbf f^\top\, .
\tag{7}
$$ The last equality can be established in the following way. Note that ${\mathsf A}:=\sum_{\mathbf f\in Q }\mathbf f\,\! \mathbf f^\top$ can be rewritten as ${\mathsf C}({\mathbf a})$ with $a_i =\mathbf q^\top{\mathsf P}_i \mathbf q$, because $$\mathbf e_r^\top\left [ \sum_{i=1}^n ({\mathsf P}_i\mathbf q)({\mathsf P}_i\mathbf q)^\top\right]\mathbf e_s = \sum_{i=1}^n q_{i\oplus r} q_{i\oplus s} = \sum_{i=1}^n q_i q_{i \oplus r \ominus s}= \mathbf q^\top{\mathsf P}_{r \ominus s} \mathbf q$$ depends on $(r,s)$ only via $r\ominus s$. Symmetry of ${\mathsf A}={\mathsf C}({\mathbf a})$ follows from Lemma 3.1 because condition (5) is satisfied due to $$a_i = \mathbf q^\top{\mathsf P}_i \mathbf q=\mathbf q^\top{\mathsf P}_{n-i}^{-1} {\mathsf P}_{n-i} {\mathsf P}_i \mathbf q= \mathbf q^\top{\mathsf P}_{n-i}^\top\mathbf q= a_{n-i}\quad \mbox{for all }i\in [1\!:\! n-1]\, .$$ Equality (7) is now established by Lemma 3.2(a). Hence the rank of $\mathcal{F}$ can be at most $2n-1$. Let $\mathbf q_i:={\mathsf P}_i\mathbf q$ for $i\in[1\!:\! n]$. Then $\{r\ominus i,s\ominus i\}\subseteq I(\mathbf q_i)$. Further define $\mathbf q'_i:={\mathsf R}\mathbf q_i={\mathsf R}{\mathsf P}_i\mathbf q={\mathsf P}_{n-i}{\mathsf R}\mathbf q$ and note that $\mathbf e_a^\top\mathbf q_b= q_{a\oplus b}$ as well as $\mathbf e_a^\top\mathbf q'_b= q_{1 \oplus b\ominus a}$ for all $a,b\in [1\!:\! n]$. Next consider the equation $$
\sum_{i=1}^n x_i\mathbf q_i\mathbf q_i^\top+\sum_{i=1}^n x'_i\mathbf q'_i{\mathbf q'_i}^\top={\mathsf O}.
\tag{8}
$$ Multiplying with $\mathbf e_{r\ominus j}^\top$ from the left and with $\mathbf e_{s\ominus j}$ from the right, we obtain $$\sum_{i=1}^n x_i q_{r\ominus j \oplus i} q_{s \ominus j \oplus i} +\sum_{i=1}^n x'_i q_{1\oplus i \ominus \left ( r \ominus j \right )} q_{1 \oplus i \ominus \left ( s \ominus j \right )}=0.$$ By the assumptions on $U_{\mathbf q}$ we see that the only terms contributing to the sum are achieved by choosing $i=j$ in the first term and $i=  r \oplus s\ominus j \ominus 1$ in the second term. This results in

$$
x_j q_r q_s + x'_{r \oplus s\ominus j \ominus 1} q_s q_r =0 \quad \textup{ for all }j\in[1:n]\, .
\tag{9a}
$$ Similarly multiply with $\mathbf e_{\rho\ominus j}^\top$ from the left and with $\mathbf e_{\sigma\ominus j}$ from the right, yielding $$
x_j q_\rho q_\sigma + x'_{\rho \oplus \sigma\ominus j \ominus 1} q_\sigma q_\rho =0 \quad \textup{ for all }j\in[1\!:\! n]\, .
\tag{9b}
$$ From these equations we conclude that $x_j'=x'_{j\oplus \rho\oplus \sigma\ominus r\ominus s} =x'_{j\oplus (\rho+\sigma-r-s)}$ for all $j\in[1\!:\! n]$. Fixing $x_1'=\xi$, and employing coprimality of $\rho+\sigma-r-s$ and $n$, we see that our system $9$ of $2n$ equations has the unique solution $x_i=-x_i'=-\xi$ for $i\in[1\!:\! n]$. So there is a one parameter family of solutions parameterized by $\xi$, showing that if any of the coefficients in (8) is zero, all others also must be zero, so indeed every $(2n-1)$-element subset of $\mathcal{F}$ has to be linearly independent, as asserted. $\Box$

### Theorem 3.2

*Let $\mathbf q$ satisfy the hypotheses of Lemma 3.3 and define the (therefore disjoint) sets $Q:=\{{\mathsf P}_i\mathbf q:i\in[1\!:\! n]\}$ and $Q':=\{{\mathsf R}\mathbf q:\mathbf q\in Q\}$. Consider the matrix ${\mathsf M}= \!\! \!\! \sum\limits_{\mathbf f\in Q\cup Q'}\!\! y_\mathbf f\,\mathbf f\,\!\mathbf f^\top$ with $\mathbf y\in\mathbb{R}_+^{|Q\cup Q'|}$, and assume that ${\mathsf M}\in\mathcal{E}(Q\cup Q')$ holds. Then we have:*

1.  *If all $y_\mathbf f>0$ and if $|\textup{Argmin} \left\{y_\mathbf f:\mathbf f\in Q\right \} |=|\textup{Argmin} \left\{y_\mathbf f:\mathbf f\in Q'\right \} |=1$, then there are exactly two different minimal cp decompositions of ${\mathsf M}$ and $\mbox{cpr }{\mathsf M}=2| Q |-1$.*

2.  *If $y_\mathbf f=0$ for at least one $\mathbf f\in Q$ and at least one $\mathbf f\in Q'$, then the minimal cp decomposition of ${\mathsf M}$ is unique and $\mbox{cpr }{\mathsf M}=| I(\mathbf y)|$.*

*Proof.*
Define $u_\mathbf f:=1$ for all $\mathbf f\in Q$ and $u_\mathbf f:=-1$ for all $\mathbf f\in Q'$. Then, by Lemma 3.3 and equation (7), the solutions $\mathbf x$ of the equation ${\mathsf M}= \sum\limits_{\mathbf f\in Q\cup Q'}x_\mathbf f\,\mathbf f\,\!\mathbf f^\top$ are given by $\mathbf x=\mathbf y+\xi\mathbf u$. In case (a), the solutions $\mathbf x\ge\mathbf o$ additionally require $\xi\in[-\min\{y_\mathbf f:\mathbf f\in Q\},\min\{y_\mathbf f:\mathbf f\in Q'\}]$, with $|I(\mathbf x)|=2|Q|-1$ (resp. $|I(\mathbf x)|=2|Q|$) for $\xi$ on the boundary (resp. in the interior) of that interval. In case (b), the condition $\mathbf x\ge\mathbf o$ is violated for any $\xi\ne0$, so $\mathbf x=\mathbf y$ is unique. $\Box$

## 4. Counterexamples to the Drew-Johnson-Loewy Conjecture

For the examples to follow, we selected matrices ${\mathsf S}$ with integer entries, where we could determine all minimizers of the quadratic form $\mathbf x^\top{\mathsf S} \mathbf x$ by exact arithmetic, solving the first-order conditions and checking the values for nonnegativity with the help of (6), cf. also [3; 4]. To be more precise, we first checked (by exact arithmetic to avoid any numerical errors) for all possible supports $I\subseteq [1:n]$ with $I\neq\emptyset$, whether there could be a *local* solution $\mathbf q$ to the optimization problem $\min\limits_{\mathbf q\in \Delta} \mathbf q^\top{\mathsf S} \mathbf q$ with $I(\mathbf q)=I$. To this end, ignoring the variables $q_i=0$ for $i\in [1\!:\! n]\setminus I$, we see there is only one locally binding constraint, namely $\mathbf e^\top\mathbf q=1$. So, denoting the multiplier of this constraint by $2m$, we arrive at the first-order conditions $$
\left. \begin{array}{rcl}  {\mathsf S}_I\mathbf x&= &m\mathbf e\\[0.2em]
                                \mathbf e^\top\mathbf x&= &1  \end{array}\right \}\, ,
\tag{10}
$$ where ${\mathsf S}_I$ denotes the principal submatrix of ${\mathsf S}$ on $I\times I$. Since all constraints of the optimization problem $\min\limits_{\mathbf q\in \Delta} \mathbf q^\top{\mathsf S} \mathbf q$ are linear, any local minimizer of $\mathbf q^\top{\mathsf S} \mathbf q$ over $\Delta$ must solve (10) for some $I$, putting $\mathbf x=[q_i]_{i\in I}$ and $m = m\mathbf x^\top\mathbf e= \mathbf x^\top{\mathsf S}_I \mathbf x = \mathbf q^\top{\mathsf S} \mathbf q$. Now suppose $\mathbf e^\top\mathbf x=\mathbf e^\top\mathbf y=1$ and ${\mathsf S}_I\mathbf x=m\mathbf e$ and ${\mathsf S}_I\mathbf y=t\mathbf e$. Then $$
t= ( t\mathbf e)^\top\mathbf x= \mathbf y^\top{\mathsf S}_I^\top\mathbf x= \mathbf y^\top{\mathsf S}_I\mathbf x= \mathbf y^\top(m\mathbf e) = m \, ,
\tag{11}
$$ so that the value $m = \mathbf x^\top{\mathsf S}_I \mathbf x =: m_I$ at any solution $(m,\mathbf x)\in {\mathbb R}\times {\mathbb R}^{|I|}$ to (10) is uniquely determined by $I$. We solved (10) by exact arithmetic for all $I$. If there is a unique solution $(m_I,\mathbf x_I)$, we discarded $I$ where $\mathbf x_I \notin {\mathbb R}^{|I|}_+$. For the remaining $I$, we confirmed that $m_I\ge0$ if (10) has a solution at all. This established copositivity of ${\mathsf S}$. The next step is to determine all zeroes of ${\mathsf S}$, i.e., all solutions $(0,\mathbf x)$ to (10) with $\mathbf x\in {\mathbb R}^{|I|}_+$. While there could be multiple solutions to (10) for $m_I>0$, this is ruled out in case $m_I=0$ for the matrices ${\mathsf S}$ considered below. Indeed, consider again two solutions $(0,\mathbf x)$ and $(0,\mathbf y)$ to (10). Then ${\mathsf S}_I(\mathbf x-\mathbf y) = \mathbf o$ and $\mathbf x-\mathbf y\in \mbox{ker }{\mathsf S}_I \cap \mathbf e^\perp$, so that the condition $\mbox{ker }{\mathsf S}_I\cap \mathbf e^\perp=\left\{\mathbf o\right \}$ rules out multiple solutions to (10); this is in fact true for any value of $m_I$, due to (11). Now [4, Lemma 1] shows that $\mbox{ker }{\mathsf S}_I\cap \mathbf e^\perp=\left\{\mathbf o\right \}$ holds if $\mathbf e\mathbf e^\top- {\mathsf S}_I$ is nonsingular, which we confirmed (again by exact arithmetic) for all $I$ which admit a solution $(0,\mathbf x)$ to (10) with $\mathbf x\in {\mathbb R}^{|I|}_+$. Note that if $(0,\mathbf x)$ solves (10), then $(\mathbf e\mathbf e^\top-{\mathsf S}_I)\mathbf x= \mathbf e$ and $\mathbf e^\top\mathbf x=1$ implies $\mathbf e^\top(\mathbf e\mathbf e^\top-{\mathsf S}_I)^{-1} \mathbf e=1$. So we considered the unique solution $\mathbf x_I := (\mathbf e\mathbf e^\top-{\mathsf S}_I)^{-1} \mathbf e\in {\mathbb R}^{|I|}_+$. Finally, we filled the entries with indices in $[1\!:\! n]\setminus I$ by zeros to get a vector which we call $\mathbf q_I\in {\mathbb R}^n$ and collected these as the set of all zeroes of ${\mathsf S}$. In this way the assumptions of the previous sections were ensured. As a final remark, note that [4, Lemma 1] says that $\mbox{ker }{\mathsf S}_I\cap \mathbf e^\perp\neq\left\{\mathbf o\right \}$ holds if and only if *both* ${\mathsf S}_I$ *and* ${\mathsf S}_I-\mathbf e\mathbf e^\top$ are singular; however, as noted by the Senior Editor, in case $m_I=0$ the principal submatrices ${\mathsf S}_I$ necessarily have to be singular as they must be positive-semidefinite (see, e.g. [9, Lemma 2.4]), so the above argument involving ${\mathsf S}_I-\mathbf e\mathbf e^\top$ is essential.

**Example 1 $(p_7\ge14)$:** Let ${\mathsf S}={\mathsf C}([-153,127,-27,-27,127,-153,162]^\top)$. Then the set of zeroes of ${\mathsf S}$ in $\Delta \subset\mathbb{R}^7$ consists of 14 vectors: $\mathbf q_i={\mathsf P}_i\mathbf u, i\in[1\!:\! 7]$, where $\mathbf u=\frac17[3,3,0,0,1,0,0]^\top$, and $\mathbf q_i={\mathsf P}_i\mathbf v, i\in[8\!:\! 14]$, where $\mathbf v=\frac1{35}[9,17,9,0,0,0,0]^\top$. Let $${\mathsf M}:=\sum_{i=1}^{14}\mathbf q_i\mathbf q_i^\top={\textstyle\frac1{1225}}\,{\mathsf C}([531,81,150,150,81,531,926]^\top).$$ The difference sets of $I(\mathbf u)=\{1,2,5\}$ and $I(\mathbf v)=\{1,2,3\}$ are $$D_\mathbf u=\{1,3,4,6\},\quad U_\mathbf u=\{1,6\},\quad D_\mathbf v=\{1,2,5,6\},\quad U_\mathbf v=\{2,5\}.$$ We note that $d=2\in U_\mathbf v\setminus D_\mathbf u$, so we may apply Theorem 3.1 with $Q=\{\mathbf q_1,\ldots,\mathbf q_{14}\}$, $\mathbf q=\mathbf v$ and $Q_1=\{\mathbf q_8,\ldots,\mathbf q_{14}\}$, to conclude that in any cp decomposition ${\mathsf M}=\sum\limits_{i=1}^{14}x_i\mathbf q_i\mathbf q_i^\top$ we must have $x_k=1$ for $k\in[8\!:\! 14]$. Moreover Theorem 3.1 states that $\widehat{\mathsf M}:={\mathsf M}-\sum\limits_{i=8}^{14}\mathbf q_i\mathbf q_i^\top=\sum\limits_{i=1}^{7}x_i\mathbf q_i\mathbf q_i^\top$ satisfies $\widehat{\mathsf M}\in\mathcal{E}(\widehat Q)$, where $\widehat Q=\{\mathbf q_1,\ldots,\mathbf q_{7}\}$. Therefore $d=1\in U_{\mathbf u}$ allows us to invoke Theorem 3.1 with ${\mathsf M}=\widehat{\mathsf M}$, $Q=\widehat Q$ and $\mathbf q=\mathbf u$. We conclude that $x_k=1$ also for $k\in[1\!:\! 7]$, that ${\mathsf M}$ has a unique minimal cp decomposition, and that $\mbox{cpr }{\mathsf M}=14$. Another matrix of this sort, having small integer entries, is $$\widetilde M_7:={\tfrac23} 7^2\sum_{i=1}^{7}\mathbf q_i\mathbf q_i^\top+{\tfrac13} 35^2\sum_{i=8}^{14}\mathbf q_i\mathbf q_i^\top
 =
 {\left [
  \begin{array}{ccccccc}
  163&108&27&4&4&27&108\\108&163&108&27&4&4&27\\27&108&163&108&27&4&4\\4&27&108&163&108&27&4\\4& 4&27&108&163&108&27\\27&4&4&27&108&163&108\\108&27&4&4&27&108&163
  \end{array}
\right]}
 \, .$$ Note that both, above matrix and ${\mathsf M}$, have no zero entries and full rank.

**Example 2 $(p_9\ge26)$:** Let $${\mathsf S}={\mathsf C}([-1056, 959, -484, 231, 231, -484, 959, -1056, 1089]^\top).$$ Then the set of zeroes of ${\mathsf S}$ in $\Delta\in\mathbb{R}^9$ consists of 27 vectors: indeed, let $$\left . \begin{array}{ccc}
\!\!\mathbf u\!\! &=\!\! &\frac1{26}[11,12,0,0,3,0,0,0,0]^\top\\[0.2em]
\!\!\mathbf v\!\! &=\!\! &\frac1{26}[12,11,0,0,0,0,3,0,0]^\top\\[0.2em]
\!\!\mathbf w\!\!&= \!\!&\frac1{130}[33,64,33,0,0,0,0,0,0]^\top\end{array}\!\!\!\!\right \} \mbox{ and define }\; \mathbf q_i := \left\{\begin{array}{ll}
\!\!{\mathsf P}_i\mathbf u\, ,\!\! &\mbox{if }i\in[1\!:\! 9]\, ,\\[0.2em]
\!\!{\mathsf P}_i\mathbf v\, , \!\! &\mbox{if }i\in[10\!:\! 18]\, ,\\[0.2em]
\!\!{\mathsf P}_i\mathbf w\, ,\!\!   &\mbox{if } i\in[19\!:\! 27]\, .
\end{array}\right .$$ The set of zeroes of ${\mathsf S}$ is $\left\{\mathbf q_i : i\in [1\!:\! 27]\right \}$ and ${\mathsf P}_2\mathbf v={\mathsf R}\mathbf u\notin\left\{{\mathsf P}_i\mathbf u:i\in [1\!:\! 9]\right \}$. Put $$\begin{align*}
{\mathsf M}:=&2\sum_{i=1}^{18}\mathbf q_i\mathbf q_i^\top-\mathbf q_{9}\mathbf q_{9}^\top-\mathbf q_{11}\mathbf q_{11}^\top+\sum_{i=19}^{27}\mathbf q_i\mathbf q_i^\top\\
=&{\tfrac1{16900}}
{
\left [
  \begin{array}{ccccccccc}
30649&14124&1089&3600&2475&3300&3600&1089&17424\\14124&30074&17424&1089&2700&3300& 3300&3600&1089\\1089&17424&33674&17424&1089&3600&3300&3300&3600\\3600&1089&17424& 33674&17424&1089&3600&3300&3300\\2475&2700&1089&17424&33224&17424&1089&2700& 2475\\3300&3300&3600&1089&17424&33674&17424&1089&3600\\3600&3300&3300&3600&1089& 17424&33674&17424&1089\\1089&3600&3300&3300&2700&1089&17424&30074&14124\\17424&1089& 3600&3300&2475&3600&1089&14124&30649\\
   \end{array}
\right]}.
\end{align*}$$ The difference sets of $I(\mathbf u)=\{1,2,5\}$, $I(\mathbf v)=\{1,2,7\}$ and $I(\mathbf w)=\{1,2,3\}$ are $$D_\mathbf u=U_\mathbf u=D_\mathbf v=U_\mathbf v=\{1,3,4,5,6,8\},\quad D_\mathbf w=\{1,2,7,8\},\quad U_\mathbf w=\{2,7\}.$$ We note that $d=2\in U_\mathbf w\setminus (D_\mathbf u\cup D_\mathbf v)$, so we may apply Theorem 3.1 to conclude that in any cp decomposition ${\mathsf M}=\sum\limits_{i=1}^{27}x_i\mathbf q_i\mathbf q_i^\top$ we must have $x_i=1$ for $i\in[19\!:\! 27]$. Next, consider the matrix $\widehat{\mathsf M}:={\mathsf M}-\sum\limits_{i=19}^{27}\mathbf q_i\mathbf q_i^\top$, which satisfies $\widehat{\mathsf M}\in\mathcal{E}(\widehat Q)$ by Theorem 3.1, where $\widehat Q=\{\mathbf q_1,\ldots,\mathbf q_{18}\}$. Since the differences $d_1=1=2-1$ and $d_2=3=5-2$ appear only once in $U_{\mathbf u}$, and $5+2-2-1=4$ and $9$ are coprime allows to invoke Lemma 3.3 and Theorem 3.2 with ${\mathsf M}=\widehat {\mathsf M}$, $Q\cup Q'=\widehat Q$ and $\mathbf q=\mathbf u$. We conclude that there are exactly two vectors $\mathbf x\in\mathbb{R}^{18}_+$ of support of size 17, (and no such vectors of smaller support,) that give rise to minimal cp decompositions of $\widehat {\mathsf M}$, and that $\mbox{cpr }{\mathsf M}=26$. Another matrix of this sort, having small integer entries, is $$\begin{align*}
\widetilde {\mathsf M}_9:=& {\tfrac56 \, 26^2}\left(\sum_{i=1}^{18}\mathbf q_i\mathbf q_i^\top-{\tfrac35}\,(\mathbf q_{7}\mathbf q_{7}^\top+\mathbf q_{13}\mathbf q_{13}^\top)\right)
 +{\tfrac13}\, 130^2\sum_{i=19}^{27}\mathbf q_i\mathbf q_i^\top\\
=&{
\left [
  \begin{array}{ccccccccc}
      2548&1628&363&60&55&55&60&363&1628\\1628&2548&1628&363&60&55&55&60&363\\363&1628&2483&1562&363&42&22&55&60\\60&363&1562&2476& 1628&363&42&55&55\\55&60&363&1628&2548&1628&363&60&55\\55&55&42&363&1628&2476&1562&363&60\\60&55&22&42&363&1562&2483&1628& 363\\363&60&55&55&60&363&1628&2548&1628\\1628&363&60&55&55&60&363&1628&2548\\
   \end{array}
\right]}.
\end{align*}$$

Note that neither of these matrices of cp-rank 26 are cyclically symmetric, they have no zero entries and full rank.

**Example 3 $(p_8\ge18)$:** Continuing Example 2, we observe that the upper left $8\times8$-submatrix ${\mathsf S}_8$ of ${\mathsf S}$ has 18 zeroes. These are obtained by taking the first 8 coordinates of those zeroes $\mathbf q$ of ${\mathsf S}$ satisfying $\mathbf e_9^\top\mathbf q=0$. Indeed, if $\mathbf z^\top{\mathsf S}_8\mathbf z= 0$, then also $[\mathbf z^\top| 0]{\mathsf S}[\mathbf z^\top| 0]^\top= 0$, so that the zero $[\mathbf z^\top| 0]^\top$ must appear in the list of Example 2. Define the set $S_8:=\{\mathbf q\in\{\mathbf q_1,\ldots,\mathbf q_{27}\}:\mathbf e_9^\top\mathbf q=0\}$. Then the matrix ${\mathsf M}:=\sum\limits_{\mathbf q\in S_8}\mathbf q\mathbf q^\top$ satisfies $\mbox{cpr }{\mathsf M}=18$, by Theorem 3.1, Lemma 3.3 and Theorem 3.2. Moreover all entries in the last row and the last column of ${\mathsf M}$ are zero, therefore also ${\mathsf M}_8$, the upper left $8\times8$-submatrix of ${\mathsf M}$, has $\mbox{cpr }{\mathsf M}_8=18$. Again, by adjusting weights, we came up with a matrix with small integer entries: $$\widetilde {\mathsf M}_8:=
{\left [
  \begin{array}{ccccccccc}
541&880&363&24&55&11&24&0\\880&2007&1496&363&48&22&22&24\\363&1496&2223&1452&363&24&22&11\\24&363&1452&2325&1584&363&48&55\\55&48& 363&1584&2325&1452&363&24\\11&22&24&363&1452&2223&1496&363\\24&22&22&48&363&1496&2007&880\\0&24&11&55&24&363&880&541\\
   \end{array}
\right]}.$$ Note that $\widetilde {\mathsf M}_8$ is, again, not cyclically symmetric, and that it has full rank.

**Example 4 $(p_{10}\ge27)$:** Continuing Example 2, let ${\mathsf M}\in {{\cal C}^{10*}}$ be the matrix obtained from $\widetilde {\mathsf M}_9$ by appending a zero column $\mathbf o\in {\mathbb R}^9$ and completing this to a symmetric $10\times 10$ matrix by adding one row $\mathbf e_{10}^\top$ as the last one. Then, by [17, Prop.2.2], we get $$\mbox{cpr }{\mathsf M}= \mbox{cpr }\widetilde {\mathsf M}_9 + 1 = 27\, .$$

**Example 5 $(p_{11}\ge32)$:** Consider $${\mathsf S}={\mathsf C}([32, 18, 4, -24, -31, -31, -24, 4, 18, 32, 32]^\top).$$ There are 33 zeroes of ${\mathsf S}$; indeed, let $$\left . \begin{array}{ccc}
\!\!\mathbf u\!\! &=\!\! &\frac1{21}[8, 0, 3, 0, 0, 0, 10, 0, 0, 0, 0]^\top\\[0.2em]
\!\!\mathbf v\!\! &=\!\! &\frac1{21}[10, 0, 0, 0, 3, 0, 8, 0, 0, 0, 0]^\top\\[0.2em]
\!\!\mathbf w\!\!&= \!\!&\frac1{7}[2, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0]^\top\end{array}\!\!\!\!\right \} \mbox{ and define }\; \mathbf q_i := \left\{\begin{array}{ll}
\!\!{\mathsf P}_i\mathbf u\, ,\!\! &\mbox{if }i\in[1\!:\! 11]\, ,\\[0.2em]
\!\!{\mathsf P}_i\mathbf v\, , \!\! &\mbox{if }i\in[12\!:\! 22]\, ,\\[0.2em]
\!\!{\mathsf P}_i\mathbf w\, ,\!\!   &\mbox{if } i\in[23\!:\! 33]\, ,
\end{array}\right .$$ then the set of zeroes can be written as $\left\{\mathbf q_i : i\in [1\!:\! 33]\right \}$. Now put $$\begin{align*}
{\mathsf M}:=&2\sum_{i=1}^{22}\mathbf q_i\mathbf q_i^\top-\mathbf q_{11}\mathbf q_{11}^\top-\mathbf q_{13}\mathbf q_{13}^\top+\sum_{i=23}^{33}\mathbf q_i\mathbf q_i^\top\\
=&{\tfrac1{441}}
{
\left [
\begin{array}{ccccccccccc}
 781 & 0 & 72 & 36 & 228 & 320 & 240 & 228 & 36 & 96 & 0 \\
 0 & 845 & 0 & 96 & 36 & 228 & 320 & 320 & 228 & 36 & 96 \\
 72 & 0 & 827 & 0 & 72 & 36 & 198 & 320 & 320 & 198 & 36 \\
 36 & 96 & 0 & 845 & 0 & 96 & 36 & 228 & 320 & 320 & 228 \\
 228 & 36 & 72 & 0 & 781 & 0 & 96 & 36 & 228 & 240 & 320 \\
 320 & 228 & 36 & 96 & 0 & 845 & 0 & 96 & 36 & 228 & 320 \\
 240 & 320 & 198 & 36 & 96 & 0 & 745 & 0 & 96 & 36 & 228 \\
 228 & 320 & 320 & 228 & 36 & 96 & 0 & 845 & 0 & 96 & 36 \\
 36 & 228 & 320 & 320 & 228 & 36 & 96 & 0 & 845 & 0 & 96 \\
 96 & 36 & 198 & 320 & 240 & 228 & 36 & 96 & 0 & 745 & 0 \\
 0 & 96 & 36 & 228 & 320 & 320 & 228 & 36 & 96 & 0 & 845 \\
\end{array}
\right]},
\end{align*}$$ and again ${\mathsf M}$ has full rank. We get $I(\mathbf u)=\{1,3,7\}$, $I(\mathbf v)=\{1,5,7\}$, $I(\mathbf w)=\{1,4,8\}$, and we calculate $$D_\mathbf u=U_\mathbf u=D_\mathbf v=U_\mathbf v=\{2,4,5,6,7,9\},\quad D_\mathbf w=\{3,4,7,8\},\quad U_\mathbf w=\{3,8\}.$$ Analogously to Example 2 we now show that the cp-rank is $32$. Since $d=3\in U_\mathbf w\setminus (D_\mathbf u\cup D_\mathbf v)$, we must have $x_i=1$ for $i\in[22\!:\! 33]$ by Theorem 3.1. Therefore consider $\widehat{\mathsf M}:={\mathsf M}-\sum\limits_{i=22}^{33}\mathbf q_i\mathbf q_i^\top$. We can see that the differences $d_1=6=7-1$ and $d_2=4=7-3$ appear in $U_{\mathbf u}$, and knowing that $7+3-7-1=2$ and $11$ are coprime allows to invoke Lemma 3.3 and Theorem 3.2. Hence there are exactly two vectors $\mathbf x\in\mathbb{R}^{22}_+$ of support of size $21$ for $\widehat {\mathsf M}$ and this leads to a total of $32$ for $\mbox{cpr }{\mathsf M}$.

    $n$    5       6          7          8          9          10         11
  ------- ---- ---------- ---------- ---------- ---------- ---------- ----------
   $d_n$   6       9          12         16         20         25         30
   $p_n$   6    $\le 15$   $\ge 14$   $\ge 18$   $\ge 26$   $\ge 27$   $\ge 32$
   $s_n$   11      17         24         32         41         51         62

  : (Ranges for) maximal cp-rank $p_n$ of cp matrices of order $n$.

Table 1 summarizes the known bracket and consequences from above examples. A tighter upper bound $p_6\le 15$ was proved in [16, Thm.6.1], but up to now no ${\mathsf M}\in {{\cal C}^{6*}}$ with $\mbox{cpr }{\mathsf M}>9=d_6$ is known.

## Acknowledgment

The authors are indebted to an anonymous referee for valuable suggestions which significantly improved presentation, and to the Senior Editor Raphael Loewy for his diligence in the evaluation process.

## References

1. Abraham Berman and Naomi Shaked-Monderer. Remarks on completely positive matrices. *Linear and Multilinear Algebra*,
44:149--163, 1998.

2. Abraham Berman and Naomi Shaked-Monderer. *Completely positive matrices*. World Scientific Publishing Co. Inc., River
Edge, NJ, 2003.

3. Immanuel M. Bomze. Detecting all evolutionarily stable strategies. *J. Optim. Theory Appl.*, 75(2):313--329, 1992.

4. Immanuel M. Bomze. Regularity versus degeneracy in dynamics, games, and optimization: a unified approach to different
aspects. *SIAM Rev.*, 44(3):394--414, 2002.

5. Immanuel M. Bomze. Copositive optimization -- recent developments and applications. *European J. Oper. Res.*,
216:509--520, 2012.

6. Immanuel M. Bomze, Werner Schachinger, and Gabriele Uchida. Think co(mpletely )positive ! -- matrix properties, examples
and a clustered bibliography on copositive optimization. *J. Global Optim.*, 52:423--445, 2012.

7. Mark Broom, Chris Cannings, and Glenn T. Vickers. On the number of local maxima of a constrained quadratic form. *Proc.
R. Soc. Lond. A*, 443:573--584, 1993.

8. Samuel Burer. Copositive programming. In Miguel F. Anjos and Jean Bernard Lasserre, editors, *Handbook of Semidefinite,
Cone and Polynomial Optimization: Theory, Algorithms, Software and Applications*, International Series in Operations
Research and Management Science, pages 201--218. Springer, New York, 2012.

9. Peter J. C. Dickinson, Mirjam Dür, Luuk Gijben, and Roland Hildebrand. Irreducible elements of the copositive cone.
*Linear Algebra Appl.*, 439(6):1605--1626, 2013.

10. John H. Drew and Charles R. Johnson. The no long odd cycle theorem for completely positive matrices. In *Random discrete
structures*, volume 76 of *IMA Vol. Math. Appl.*, pages 103--115. 1996.

11. John H. Drew, Charles R. Johnson, and Raphael Loewy. Completely positive matrices associated with $M$-matrices. *Linear
Multilinear Algebra*, 37(4):303--310, 1994.

12. Mirjam Dür. Copositive programming --- a survey. In Moritz Diehl, Francois Glineur, Elias Jarlebring, and Wim Michiels,
editors, *Recent Advances in Optimization and its Applications in Engineering*, pages 3--20. Springer, Berlin Heidelberg
New York, 2010.

13. Roland Hildebrand. The extremal rays of the $5\times 5$ copositive cone. *Linear Algebra Appl.*, 437(7):1538--1547,
2012.

14. Roland Hildebrand. Minimal zeros of copositive matrices. Preprint, <http://arxiv.org/abs/1401.0134>, 2014.

15. Raphael Loewy and Bit-Shun Tam. CP rank of completely positive matrices of order 5. *Linear Algebra Appl.*,
363:161--176, 2003.

16. Naomi Shaked-Monderer, Abraham Berman, Immanuel M. Bomze, Florian Jarre, and Werner Schachinger. New results on the cp
rank and related properties of co(mpletely )positive matrices. *Linear Multilinear Algebra*, to appear. Also available
at: arxiv.org/abs/1305.0737, 2013.

17. Naomi Shaked-Monderer, Immanuel M. Bomze, Florian Jarre, and Werner Schachinger. On the cp-rank and minimal cp
factorizations of a completely positive matrix. *SIAM J. Matrix Anal. Appl.*, 34(2):355--368, 2013.
