# Model a1: Exact Curvature Pruning with a Full-Simplex Probe

Model `a1` combines exact curvature failures with one speculative full-simplex probe.

For a support \(S\), let \(Z\) span the zero-sum directions on that support and let

\[
H_S = Z^T A_S Z.
\]

Every ESS whose support contains \(S\) would require \(H_S\) to be strictly negative definite. Therefore, if exact arithmetic proves
that \(H_S\) is singular or has a non-negative direction, neither \(S\) nor any of its supersets can support an ESS. The support
generator can prune that complete upper cone.

The search order is cardinality one, the full simplex, and then cardinalities two through \(n-1\). Singleton supports and all ordinary
supports are checked with exact arithmetic. The full simplex alone is first probed in binary64. If its reduced Hessian is comfortably
negative definite and its candidate probabilities are comfortably positive, or if the floating computation is inconclusive, A1
checks it exactly immediately. An exact full-support ESS ends the complete search.

If binary64 makes the full support look impossible, A1 merely postpones its exact check. Any exact candidate or exact curvature
ceiling found later forbids the full support. Otherwise A1 verifies the full support exactly after enumeration, so floating point can
change scheduling but cannot remove an ESS. There are no outside-payoff inequalities on the full simplex.

During normal enumeration, only the exact integer factorization may add a forbidden support. This includes cardinality two: the
ordinary reduced Hessian is then the scalar \(A_{ii}+A_{jj}-2A_{ij}\), and an exact failure becomes active before cardinality three.
There is deliberately no separate eager all-pairs scan and no per-support floating filter.

Binary64 never authorizes pruning. `fast` and `safe` retain their public names for interface compatibility, but both use this same A1
search order; the requested method changes only fallback reporting and logging in this experiment.

This model reports ESS output, not every unstable Nash candidate. Once exact curvature fails, candidate probabilities and outside
payoffs cannot affect ESS output and are deliberately skipped. Candidate counts can therefore differ from production; ESS counts
and ESS structures are the correctness comparison.

The experiment is isolated under `models/a1`. It reuses unchanged production files and dependencies, but does not alter the main
build, package, CLI, Python module, database, or release workflow.

## Build and check

From the repository root:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S models/a1 -B models/a1/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
cmake --build models/a1/build --target fracessa_core fracessa_a1_core fracessa_a1 --parallel
python3 models/a1/verify.py
taskset -c 2 python3 models/a1/verify.py --target-seconds 0.5
```

`verify.py` loads each native module in a separate persistent process, checks representative database rows plus direct edge cases,
and compares ESS counts and structures. `--target-seconds` repeats each case enough times to report a median native time; these local
results are still experimental rather than a replacement for the maintained production benchmark protocol.

The multiword pair-pruning path has a separate quick smoke test. This dimension-65 zero game should inspect cardinalities one and two
and then stop because every larger support contains a pair rejected exactly during normal enumeration:

```bash
matrix="65#0$(printf ',0%.0s' {1..32})"
models/a1/build/fracessa_a1 fast "$matrix"
```
