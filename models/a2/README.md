# Model a2: Safe Full-Simplex Inertia Pruning

Model `a2` extends A1's exact ESS curvature pruning with an optional complete-inertia calculation on the full simplex. It has only
the `safe` method. Binary64 is used only to schedule that optional exact calculation; it never establishes a candidate, an ESS, an
inertia bound, or a pruning decision.

For a support \(S\), let \(H_S\) be the payoff matrix restricted to zero-sum directions on that support, and let

\[
r=\nu_-(H_S)
\]

be its exact number of negative eigenvalues. If \(T\) supports an ESS, then

\[
|T\cap S|\le r+1.
\]

For the full simplex \(S=[n]\), this becomes the global cardinality bound

\[
|T|\le \nu_-(H_{[n]})+1.
\]

## Search order

The normal search is:

1. check every singleton exactly;
2. if no singleton already excludes the full support, estimate its candidate status and inertia in binary64;
3. check the full support exactly when A1 would: it may be an ESS, or binary64 preparation, factorization, inertia, or probabilities are
   uncertain; also check it when A2's reliable estimated cardinality ceiling is at most \(\lceil n/2\rceil\);
4. stop immediately if that exact check proves that the full support is an ESS;
5. when exact full-support inertia was calculated, enumerate only through the certified ceiling
   \(\nu_-(H_{[n]})+1\); otherwise continue ordinary exact enumeration from cardinality two;
6. if the full support survives all smaller exact pruning, check it exactly at the end.

An uncertain floating result follows A1 and requests exact verification because uncertainty cannot safely postpone a possible
full-simplex ESS. A reliable rejection of the full-simplex ESS requests early exact factorization only when the estimated ceiling is at
most \(\lceil n/2\rceil\). For odd \(n\), accepting the ceiling value \((n+1)/2\) still removes roughly the upper third of the unpruned
support space. Skipping that optional inertia calculation only forgoes pruning. A smaller exact curvature ceiling or exact Nash
equilibrium still excludes the full support through the ordinary upward-pruning rules. Once an exact full-support factorization has
been calculated, A2 applies its exact ceiling even if it differs from the floating estimate.

During ordinary enumeration, A2 retains A1's two upward-pruning rules:

- failed exact negative definiteness excludes the support and every superset;
- an exact Nash equilibrium retains its support and excludes every strict superset.

Complete inertia from an ordinary ascending support would add no new pruning. If such a support \(S\) had
\(r\le |S|-3\), one of its smaller \((r+2)\)-subsets would already have failed exact curvature and prevented \(S\) from being generated.
The out-of-order exact full-support check is different because only its singleton subsets have been visited.

The local KKT factorization reports negative inertia after singular as well as nonsingular factorizations. This matters because a
singular full reduced Hessian may produce a strong cardinality bound even though it cannot produce a full-support candidate.

Like A1, A2 preserves ESS output but may omit unstable Nash candidates rejected by the necessary curvature condition. Candidate counts
are therefore not a production-compatibility target.

## Build and check

From the repository root:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S models/a2 -B models/a2/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
cmake --build models/a2/build --target fracessa_core fracessa_a2_core fracessa_a2 --parallel
python3 models/a2/verify.py
taskset -c 2 python3 models/a2/verify.py --target-seconds 0.5
```

`verify.py` loads production safe and A2 in separate persistent processes. It compares ESS counts and support-size structures on direct
edge cases and representative database rows. A2 accepts only `safe`:

```bash
models/a2/build/fracessa_a2 safe '3#0,0,0,0,0,0'
```

The experiment is isolated under `models/a2` and does not alter production, A1, packaging, the database, or the release workflow.
