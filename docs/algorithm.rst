Algorithm Overview
##################

The mathematical problem
************************

FracESSA analyzes the standard quadratic optimization problem

.. math::

   \max_{x\in\Delta^n} x^\mathsf{T} A x,

where ``A`` is a symmetric rational ``n``-by-``n`` matrix and

.. math::

   \Delta^n=\left\{x\in\mathbb{R}^n:\sum_{i=1}^n x_i=1,\;x_i\geq0\right\}.

The same matrix is the payoff matrix of a symmetric partnership game. The evolutionarily stable strategies are exactly the strict
local maximizers of the quadratic form on the simplex.

For a mixed strategy ``x``, define:

* ``I(x)`` as its support, the indices with ``x_i > 0``;
* ``u = x^T A x`` as its payoff;
* ``J(x)`` as its extended support, the pure strategies whose payoff against ``x`` equals ``u``.

A symmetric Nash equilibrium has ``(Ax)_i = u`` for every index in ``I(x)`` and ``(Ax)_j <= u`` outside it. FracESSA calls a
support solution that passes these conditions a *candidate*. Exact stability checks then decide whether that candidate is an ESS.

Candidate solve on one support
******************************

Let ``S`` be a proposed support of size ``s`` and choose one reference index ``m`` in ``S``. Eliminating its probability with
``x_m = 1 - sum(x_i)`` turns the equal-payoff equations into a symmetric ``(s-1)``-by-``(s-1)`` system

.. math::

   \begin{aligned}
   H y &= r,\\
   H_{ij} &= a_{ij}-a_{im}-a_{mj}+a_{mm},\\
   r_i &= a_{mm}-a_{im}.
   \end{aligned}

for ``i,j`` in ``S`` other than ``m``. FracESSA solves this system with fraction-free exact integer arithmetic, recovers ``x_m``,
and rejects the support unless every support probability is strictly positive. It then computes the common support payoff and checks
every strategy outside ``S``. A larger outside payoff rejects the support; an equal outside payoff adds that strategy to ``J(x)``.

A singular reduced system is not emitted as a candidate. This does not lose an ESS: strict local stability requires ``H`` to be
negative definite, which already implies that ``H`` is nonsingular. Candidate output is therefore an intermediate record from the
ESS search, not an enumeration of every possibly degenerate Nash equilibrium.

Search pipeline
***************

For ``n`` pure strategies there are ``2^n - 1`` nonempty supports in the worst case. FracESSA never stores this complete frontier.
It generates one support at a time, in increasing support size, and performs the following work:

1. **Prepare the matrix.** The coposit parser clears common denominators once and returns one exact integer matrix plus its common
   positive denominator. ``fast`` additionally prepares an equilibrated binary64 copy or switches the matrix to ``safe`` when that
   preparation is not trustworthy.
2. **Generate a support.** Non-circular and circular matrices use different generators, described below.
3. **Find an equilibrium candidate.** The selected candidate route solves the equilibrium equations on the support and verifies the
   probability and outside-payoff conditions. A candidate returned publicly is always verified with exact arithmetic.
4. **Prune strict supersets.** Once ``I(x)`` is an exact equilibrium support, no later strict superset needs to be searched for a new
   ESS. The generator therefore skips that subtree. Such a superset can still be an unstable equilibrium; FracESSA deliberately does
   not enumerate it as another candidate. This pruning depends on exact equilibrium status, not on whether ``x`` itself is an ESS.
5. **Decide stability exactly.** The exact stability path described below labels the candidate and updates the ESS totals.

Candidate methods
*****************

``safe`` directly solves and verifies every required support with fraction-free exact integer arithmetic. It gives a complete ESS
result under the search and pruning rule above. Public candidate probabilities and payoffs are exact rational values;
floating-point conversions never participate in those decisions.

``fast`` first applies a binary64 candidate filter. The complete matrix is converted and equilibrated once, and each reduced
symmetric support system is solved with a pivoted symmetric factorization. A matrix-wide preparation failure switches the complete
search to ``safe`` and records the reason in ``safe_fallback``. An inconclusive solve is instead retried exactly for only that support
and does not set ``safe_fallback``. The remaining floating-point probability and outside-payoff rejections are heuristic, so
``fast`` can miss a valid candidate and is not a completeness certificate. Candidates that survive the filter are rechecked exactly,
and all final stability decisions are exact.

Exact stability
***************

Let ``s = |I(x)|`` and ``k = |J(x) \\ I(x)|``. The first exact stability object is the same ``(s-1)``-by-``(s-1)`` reduced Hessian
``H`` used by the candidate solve. Its directions move probability only among strategies in ``I(x)`` while preserving a total of
one. If ``H`` is not negative definite, one of those support-only directions already violates strict local stability, so the
candidate is rejected.

When ``k = 0``, negative definiteness completes the ESS decision. Otherwise the outside best replies ``K = J(x) \\ I(x)`` introduce
``k`` nonnegative coordinates. FracESSA reuses the exact factorization of ``H`` to eliminate the unrestricted support coordinates in
one Schur complement. The result is a ``k``-by-``k`` integer-scaled reduced Bomze matrix. Multiplying this matrix by a positive
integer does not change strict copositivity, and the candidate is an ESS exactly when this smaller matrix is strictly copositive.

The strict-copositivity path delegates to the embedded ``coposit::safe`` API. coposit applies exact shared prechecks, splits the graph
of negative off-diagonal entries into independent connected components, and sends only unresolved components to a finite Dickinson
certificate traversal. Each uncovered principal matrix uses one exact fraction-free solve in the nonsingular case or one exact null
vector in the singular case. Floating-point arithmetic is never accepted as a stability certificate.

Support generation and pruning
******************************

Non-circular matrices
---------------------

The non-circular generator performs fixed-cardinality binary depth-first search. It visits support sizes from one through ``n`` and
builds each mask one bit at a time. A forbidden exact-equilibrium support prunes the recursion as soon as the partial support contains
it, so the complete subtree of strict supersets is skipped before those masks are materialized.

Circular-symmetric matrices
---------------------------

A circular-symmetric matrix is unchanged after a simultaneous cyclic shift of rows and columns. Symmetry also makes reflected
supports equivalent. Binary supports related by rotation or reflection form a *bracelet*. The circular generator produces one
fixed-cardinality bracelet representative rather than every support in its orbit.

When a representative is an exact equilibrium, the generator forbids every distinct rotation and reflection of that support. The
output row stores a ``multiplier`` equal to the number of those distinct supports represented by the row; a support with its own
symmetry can have a multiplier smaller than ``2n``. Candidate and ESS totals, as well as their support-size structures, include this
multiplier. If repeated payoff entries create further exact cyclic-index symmetries, FracESSA also avoids solving equivalent bracelet
representatives while preserving the same public dihedral rows and weighted totals.

Support representation and limits
*********************************

Every support is one ``uint64_t`` handle interpreted by the matrix-wide ``SupportContext``. Through dimension 64 it stores bits
directly. At larger dimensions it addresses a fixed-width context-owned word array. The generators retain direct scalar internals so
ordinary small games do not pay for dynamic storage.

Context-owned large masks remove the former support-storage limit, not the exponential search cost. Running time depends on dimension,
candidate placement, circular symmetry, and how early exact equilibria prune their supersets. One matrix remains single-core; use
Python multiprocessing to analyze independent matrices concurrently.
