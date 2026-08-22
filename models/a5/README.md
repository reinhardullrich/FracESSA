# Model A5: SAT Invader-Interval Pruning

Model A5 is an isolated safe-only experiment based on A4. It replaces fixed-order support generation with one incremental CaDiCaL
2.2.1 solver. Exact numerical analysis remains unchanged.

For a support `S` with exact invaders `V`, A5 learns the clause

```text
not(all strategies in S) or at least one strategy in V.
```

Unlike A4, A5 retains this exact clause at every support size. Curvature failures and equilibria add ordinary upper-cone clauses;
otherwise rejected supports receive a cardinality-expiring exact-support clause. Supports are requested in increasing cardinality.
An unsatisfiable cardinality does not end the search because a larger support may satisfy an invader clause.

Circular games use the same SAT generator. After one orbit member is analyzed, A5 rotates, reflects, and applies every detected
matrix-preserving affine relabeling to both learned clauses and candidate output. The other orbit members are therefore neither
solved nor reported individually; candidate and ESS structures use the same orbit multipliers as A2 and A4.

## Build and check

From the repository root:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S models/a5 -B models/a5/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
cmake --build models/a5/build --parallel
/home/reinhard/.local/bin/python3.12 models/a5/verify.py --target-seconds 0
```

The experiment does not modify production behavior, packaging, the database, or the release workflow.
