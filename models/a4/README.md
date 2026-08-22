# Model A4: Sideways Invader Pruning

Model A4 is an isolated safe-only experiment based on A2. When an exact internally valid support candidate on `S` fails only because
outside strategies in `V` earn strictly more, every later ESS support containing `S` must contain at least one member of `V`.
Equivalently, A4 learns the exact clause

```text
not(all strategies in S) or at least one strategy in V.
```

The retained implementation keeps the broad clauses learned from singleton and pair supports in non-circular games. The fixed-order
support generator indexes each rule by the last required bit, activates it only on branches already containing all of `S`, and checks
only activated rules before emitting a support. Circular games retain A2's unchanged bracelet generator: expanding clauses over
circular symmetry or retaining representative-only clauses measured slower than solving the surviving systems.

## Result

Correctness matched A2 on the maintained ordinary, circular, affine-circular, multiword, and database cases.

The full comparison selected every stored game with a positive production-safe calibration below 100 seconds, excluding dimension
2 as required by the benchmark policy. It used a 0.5-second target, a 30-second timeout, a scheduler on CPU 1, and persistent workers
on CPUs 2 through 9. The unchanged A2 measurements were reused; only A4 was rerun.

- A4 returned the correct ESS count for 1,070 of 1,071 matrices, with no mismatch and one timeout.
- A2 had 1,069 correct results and two timeouts. A4 solved matrix 2453 in 0.169 ms after A2 had timed out at 30 seconds; both timed
  out on matrix 2493.
- Across the 1,069 paired successful matrices, the median A4/A2 ratio was 1.026, or 2.57% slower.
- The non-circular median was 7.11% slower. The circular median was unchanged because circular games use A2's generator.
- Matrix 2424 improved from 3.163 seconds to 0.119 ms, matrix 52 from 1.961 seconds to 0.136 ms, and matrix 2430 from 1.908 seconds
  to 0.141 ms.
- The largest absolute regression was matrix 211, from 0.787 seconds to 1.116 seconds: 41.9% or 0.330 seconds slower.
- The largest ratio regression was matrix 2169 at 77.98 times slower, but its absolute runtime increased only from 0.374 ms to
  29.196 ms.

These are per-matrix native timings from the multi-CPU experiment, not a sum of speedups. A4 is retained because the complete set
shows a small median cost, one removed timeout, and several exceptional improvements; it is not yet a production change.

Earlier complete indexes were discarded:

- a dynamic BDD made representative cases up to roughly 1,163 times slower;
- a general two-watched-literal implementation made matrix 89 roughly 320 times slower;
- expanding circular clauses over rotations, reflections, and affine symmetries made matrix 34 roughly 1,927 times slower.

## Build and check

From the repository root:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S models/a4 -B models/a4/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
cmake --build models/a4/build --parallel
/home/reinhard/.local/bin/python3.12 models/a4/verify.py --target-seconds 0
```

The experiment does not modify production behavior, packaging, the database, or the release workflow.
