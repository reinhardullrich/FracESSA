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

A symmetric Nash equilibrium has ``(Ax)_i = u`` for every index in ``I(x)`` and ``(Ax)_j <= u`` outside it. FracESSA calls every
equilibrium found in this way a *candidate*. Exact stability checks then decide whether that candidate is an ESS.

Search pipeline
***************

For ``n`` pure strategies there are ``2^n - 1`` nonempty supports in the worst case. FracESSA never stores this complete frontier.
It generates one support at a time, in increasing support size, and performs the following work:

1. **Prepare the matrix.** The parser preserves exact rational values. The exact solver clears common denominators once and works
   with integers. ``fast`` additionally prepares an equilibrated binary64 copy or switches the matrix to ``safe`` when
   that preparation is not trustworthy.
2. **Generate a support.** Non-circular and circular matrices use different generators, described below.
3. **Find an equilibrium candidate.** The selected candidate route solves the equilibrium equations on the support and verifies the
   probability and outside-payoff conditions. A candidate returned publicly is always verified with exact arithmetic.
4. **Prune strict supersets.** Once ``I(x)`` is an exact equilibrium support, no later strict superset can be another candidate.
   The generator therefore skips that subtree. This pruning depends on equilibrium status, not
   on whether ``x`` is an ESS.
5. **Decide stability exactly.** The exact stability path described below labels the candidate and updates the ESS totals.

Candidate methods
*****************

``safe`` directly solves and verifies every support with fraction-free exact integer arithmetic. Exact rational probabilities and
payoff are constructed only for successful candidates that must be exposed as output.

``fast`` first applies a binary64 candidate filter. The complete matrix is converted and equilibrated once, and each reduced
symmetric support system is solved with a pivoted symmetric factorization. A matrix-wide preparation failure or an inconclusive
support solve falls back to exact arithmetic. The remaining floating-point probability and outside-payoff rejections are heuristic,
so ``fast`` can miss a valid candidate and is not a completeness certificate. Candidates that survive the filter are rechecked
exactly, and all final stability decisions are exact.

Exact stability
***************

The first exact stability object is the reduced Hessian on directions that stay inside ``I(x)`` while preserving the sum of the
probabilities. If it is not negative definite, the candidate cannot be a strict local maximizer and is rejected.

When ``J(x) = I(x)``, negative definiteness completes the ESS decision. Otherwise the outside best replies
``K = J(x) \\ I(x)`` introduce additional feasible directions. FracESSA uses the already available exact Hessian factorization to
form an integer-scaled reduced Bomze matrix by a Schur complement. The candidate is an ESS exactly when the required strict
copositivity condition for this smaller matrix holds.

The strict-copositivity path is also exact. It first applies low-dimensional and sign-based decisions, splits the graph of negative
off-diagonal entries into independent connected components, and sends only unresolved components to Hadeler enumeration. Hadeler
checks principal submatrices in increasing cardinality. Each unresolved matrix uses an exact fraction-free determinant, one retained
solve in the nonsingular case, or an exact nullspace in the singular case. Floating-point arithmetic is never accepted as a
stability certificate.

Support generation and pruning
******************************

Non-circular matrices
---------------------

The non-circular generator performs fixed-cardinality binary depth-first search. A forbidden exact-equilibrium support prunes the
recursion as soon as the partial support must contain it, so all strict supersets are skipped before they are materialized.

Circular-symmetric matrices
---------------------------

A circular-symmetric matrix is unchanged after a simultaneous cyclic shift of rows and columns. Symmetry also makes reflected
supports equivalent. Binary supports related by rotation or reflection form a *bracelet*. The circular generator produces one
fixed-cardinality bracelet representative rather than every support in its orbit.

When a representative is an exact equilibrium, the generator forbids the corresponding symmetry orbit. The output row stores a
``multiplier`` equal to the number of distinct rotations and reflections represented by that row. Candidate and ESS totals, as well
as their support-size structures, include this multiplier. Matrices whose repeated entries create further exact affine symmetries
can be reduced by those symmetries as well.

Support representation and limits
*********************************

Dimensions through 64 use one ``uint64_t`` per support. Larger dimensions use a fixed-width vector of 64-bit words selected once for
the matrix. Both representations expose the same generator operations; the one-word route remains specialized so ordinary small
games do not pay for dynamic storage.

Multiword masks remove the former support-storage limit, not the exponential search cost. Running time depends on dimension,
candidate placement, circular symmetry, and how early exact equilibria prune their supersets. One matrix remains single-core; use
Python multiprocessing to analyze independent matrices concurrently.
