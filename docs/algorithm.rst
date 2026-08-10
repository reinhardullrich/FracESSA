Algorithm Overview
##################

For a game with ``n`` pure strategies, FracESSA may need to inspect every nonempty support, up to ``2^n - 1`` possibilities.
It processes each support in three steps:

1. Solve the equilibrium equations on that support. Safe search uses exact fraction-free integer arithmetic; fast search first tries
   a binary64 filter and falls back to the exact solver when the filter cannot decide safely.
2. Verify the equilibrium inequalities and prune every later strict superset of an exact equilibrium support.
3. Test ESS stability exactly. Reduced-Hessian inertia handles the ordinary case. Unresolved outside best replies lead to a scaled
   reduced Bomze matrix and an exact strict-copositivity test.

Non-circular games generate supports by cardinality with binary depth-first search. Circular-symmetric games have the same entries
after a simultaneous cyclic shift of rows and columns. Their rotations and reflections are equivalent, so FracESSA generates one
bracelet representative and records how many supports it represents.

Dimensions through 64 use one machine word per support. Larger dimensions use a multiword mask. This removes the storage limit, not
the exponential running time.
