# Support-Frontier Generation And Pruning Experiment

Date: 2026-07-29

## Scope

This experiment tested alternatives to eagerly storing every nonempty support
and repeatedly erasing strict supersets with `std::remove_if`. Production C++
source was not changed. All alternative headers, binaries, generated matrices,
and raw results live under `experiments/support_frontier_2026-07-29/`.

Correctness remained the first requirement: every retained generator preserved
cardinality order, increasing numeric mask order, circular representatives, and
candidate IDs.

## Test Corpus And Protocol

- All 44 maintained verification matrices were used for differential checking.
- Thirty-five deterministic non-circular symmetric integer games were generated:
  five each for dimensions 18 through 24. Each upper triangle used Python seed
  `(dimension << 16) + sample`, entries uniformly drawn from `[-9,9]`, and
  contained positive, negative, and zero values. The requested dimensions 19
  through 24 therefore contributed 30 matrices.
- Timed runs parsed each matrix once, constructed a fresh analyzer with candidate
  output enabled per repetition, warmed up, calibrated repetitions, used up to
  seven batches, reported the median, and pinned execution to CPU 2.
- Representative final comparisons ran each binary for approximately three
  seconds and alternated production/adaptive execution order.

## Alternatives Tested

1. Fixed-cardinality Gosper generation with deferred `remove_if` pruning.
2. Directly marking every strict superset in a bitmap.
3. Gosper block jumps that searched the forbidden-support list for a witness.
4. A byte-indexed jump table storing the best lowest-bit jump at each mask.
5. A split high-prefix/low-suffix fixed-cardinality generator.
6. Recursive fixed-cardinality DFS that rejects a branch when it completes a
   forbidden subset.
7. An adaptive combination of split generation, bounded DFS, byte-indexed
   Gosper jumps, and deferred `remove_if`.

The Gosper block jump was checked exhaustively for every single forbidden subset
through dimension 12 and for 2,000 random overlapping forbidden families through
dimension 20. Its emitted sequence always matched filtered ordinary Gosper order.

## Final Paired Timing

Negative change is faster. These are the final paired three-second measurements
after rebuilding the experimental archive correctly from exactly one copy of
each object file.

| ID | Dimension | Production | Adaptive | Change |
|---:|---:|---:|---:|---:|
| 9 | 5 | 10.090 us | 10.128 us | +0.38% |
| 15 | 10 | 258.520 us | 268.719 us | +3.95% |
| 17 | 11 | 46.933 us | 45.674 us | -2.68% |
| 23 | 14 | 1.804 ms | 1.822 ms | +1.03% |
| 27 | 18 | 69.239 ms | 47.621 ms | -31.22% |
| 30 | 20 | 360.342 ms | 227.109 ms | -36.97% |
| 31 | 21 | 761.565 ms | 488.322 ms | -35.88% |
| 32 | 22 | 3.596 s | 1.365 s | -62.04% |
| 33 | 23 | 571.177 ms | 572.619 ms | +0.25% |
| 34 | 24 | 22.686 s | 5.966 s | -73.70% |

The 30 generated dimensions 19-24 used a shorter 0.5-second screening protocol.
Adaptive generation won 18 and lost 12 cases: mean `-3.58%`, median `-2.25%`,
best `-21.24%`, and worst `+2.67%`. Longer reruns confirmed that weakly pruned
cases can regress by approximately 4-5%, while the strongest sampled win remained
approximately 22%.

## Memory

One matrix-34 analyzer run measured with `/usr/bin/time` used:

| Variant | Peak RSS | Elapsed |
|---|---:|---:|
| Production | 108,336 KiB | 21.99 s |
| Adaptive | 55,568 KiB | 5.77 s |
| Byte-indexed Gosper jump | 55,504 KiB | 5.83 s |

The adaptive variant reduced peak RSS by 48.7%. The byte jump table itself still
has `2^n` bytes, so it postpones rather than removes the exponential memory
ceiling.

## Dimension-19 Circular Follow-up

The initial representative table omitted official IDs 28 and 29. Both are
circular-symmetric games at prime dimension 19. A paired three-second follow-up
found materially different behavior:

| ID | Candidate support sizes | Production | Byte-indexed Gosper | Change |
|---:|---|---:|---:|---:|
| 28 | `7:1444, 8:399` | 11.918 ms | 18.124 ms | +52.07% |
| 29 | `7:19, 10:57` | 90.910 ms | 92.997 ms | +2.30% |

For prime dimension 19, every cardinality from 1 through 18 is coprime to the
dimension. The circular search therefore stores only one canonical mask from
each 19-member rotational orbit. ID 28 nevertheless creates all 19 rotated
candidates for each solved representative and invokes superset pruning for every
rotation. Its size-7 and size-8 candidates fall inside the direct-marking window.
The byte implementation consequently attempts approximately 6.73 million
superset-table updates: `1444*(2^12-1) + 399*(2^11-1)`. ID 29 attempts only
`19*(2^12-1) = 77,805` direct updates because its size-10 candidates are outside
that window.

Prime dimension also helps explain ID 33 at dimension 23: production receives
maximum circular representative reduction, while its candidates have sizes 13
and 14 and therefore provide no direct-jump benefit. Byte-indexed Gosper was
2.96% slower there. Prime dimension is thus relevant but not sufficient; the
candidate support sizes and number of rotational copies determine whether the
direct table saves or wastes work.

## Complete Official-Matrix Split

A final paired three-second pass covered all 44 verification matrices that
existed at that point and verified byte-identical candidate output afterward.
The exact per-matrix table is stored in
`experiments/support_frontier_2026-07-29/results/gosper_jump_byte_all44_3s.tsv`.

The 17 ordinary full-matrix fixtures have dimensions only 1 through 5. Their
changes ranged from `-6.15%` to `+0.77%`, but every absolute median was below
11 microseconds, so this is primarily fixed-overhead variation. The larger
generated non-circular corpus remains the relevant evidence for that path.

Among the 27 circular fixtures, the major wins were IDs 27, 30, 31, 32, and 34
at composite dimensions 18, 20, 21, 22, and 24: `-30.31%`, `-36.32%`,
`-35.32%`, `-61.53%`, and `-73.43%`. Prime circular cases exposed the opposite
risk: ID 22 at dimension 13 regressed `13.64%`, ID 26 at dimension 17 regressed
`4.35%`, ID 28 at dimension 19 regressed `52.39%`, and ID 33 at dimension 23
regressed `2.80%`. ID 29 at the same prime dimension 19 regressed only `2.18%`,
confirming that primality explains maximum orbit reduction but does not by
itself predict the amount of direct-marking work.

## High-Dimension Non-Circular Follow-up

The original 17 non-circular fixtures stopped at dimension 5. Seven larger
symmetric matrices were therefore added as verification IDs 48-54; IDs 45-47
remain reserved for the certified-filter counterexamples.

- ID 48 is the dimension-15 Hilbert matrix,
  `H(i,j) = 1/(i+j-1)` in one-based notation. It is a standard poorly
  conditioned test matrix.
- ID 49 is the dimension-16 Sylvester Hadamard matrix. It is symmetric, has
  entries `+1/-1`, and satisfies `H H^T = 16 I`.
- ID 50 is the dimension-18 symmetric Paley conference matrix for `q=17`. It
  has zero diagonal, `+1/-1` off diagonal, and satisfies `C C^T = 17 I`.
- ID 51 is the dimension-20 MINIJ matrix, `A(i,j) = min(i,j)` in one-based
  notation. It is symmetric positive definite.
- ID 52 is the dimension-21 Fiedler matrix for `x=(1,...,21)`,
  `A(i,j) = |x_i-x_j|`. It has one positive eigenvalue and all remaining
  eigenvalues negative.
- ID 53 is a deterministic dimension-23 symmetric matrix with entries sampled
  from `[-9,9]` using seed `(23 << 16) + 101`; it contains negative, zero, and
  positive values.
- ID 54 is a deterministic dimension-24 symmetric positive-entry matrix with
  entries sampled from `[1,9]` using seed `(24 << 16) + 102`.

Definitions and properties were taken from the
[MathWorks Hilbert documentation](https://www.mathworks.com/help/matlab/ref/hilb.html),
[MathWorks test-matrix gallery](https://www.mathworks.com/help/matlab/ref/gallery.html),
[SciPy Hadamard documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hadamard.html),
and the standard
[conference-matrix definition](https://arxiv.org/abs/2004.05829).

The paired three-second result is below. Negative change is faster. A separate
CLI comparison found byte-identical complete candidate output between
production and the byte-indexed Gosper executable for every row; the raw timing
table is
`experiments/support_frontier_2026-07-29/results/high_n_non_circular_3s.tsv`.

| ID | Family | Dimension | Production | Byte-indexed Gosper | Change |
|---:|---|---:|---:|---:|---:|
| 48 | Hilbert | 15 | 87.521 ms | 87.489 ms | -0.04% |
| 49 | Hadamard | 16 | 0.384 ms | 0.399 ms | +3.77% |
| 50 | Paley conference | 18 | 97.313 ms | 98.273 ms | +0.99% |
| 51 | MINIJ | 20 | 4.296 ms | 1.474 ms | -65.70% |
| 52 | Fiedler | 21 | 814.447 ms | 822.040 ms | +0.93% |
| 53 | Random mixed | 23 | 66.393 ms | 70.671 ms | +6.44% |
| 54 | Random positive | 24 | 292.146 ms | 302.989 ms | +3.71% |

This follow-up strengthens both sides of the result. MINIJ provides a large
non-circular win from on-demand generation and early exhaustion; its singleton
candidates are outside the direct byte-marking window, so this case does not
measure a byte-table benefit. The two random matrices regress by 3.7-6.4%,
confirming that the design still needs a cheap gate for weakly pruned games.

## DFS Gate Result

A dimension-only DFS gate does not predict performance. Pure DFS versus split
generation ranged from `-4.66%` at dimension 14 to `+32.29%` at dimension 13
and `+51.38%` at dimension 17. The useful property is the discovered forbidden
family, which is matrix-specific. Adding a more elaborate runtime predictor was
not justified by the measurements.

## Correctness Result

The final adaptive executable produced byte-identical candidate output to the
production executable for all 44 official matrices and all 35 generated
non-circular matrices. Candidate counts and IDs therefore remained unchanged.

## Build Note

An intermediate manual rebuild used `ar qc` repeatedly, creating three archive
members with each object name. Those adaptive gate measurements were discarded
because the linker could select stale members. The final archive was recreated
with `ar rcs` and verified to contain exactly four unique objects. The retained
random screen predates the duplicate appends and used the original clean
archive; the final paired timing and both differential checks were run after the
corrected rebuild.

## Conclusion

The important speedup is real: direct superset marking removes the repeated
large `remove_if` scans, and byte-indexed Gosper jumps add roughly another 2-6%
on the strongly pruned games while halving matrix-34 memory. Pure DFS is not a
reliable standalone replacement.

The 235-line adaptive header is not the best production tradeoff. The simpler
157-line byte-indexed Gosper variant captures the major large-matrix speed and
memory gains, but it is not production-ready unchanged: weak random games have a
roughly 2% median regression and circular prime ID 28 is 52% slower. A retained
design must avoid repeated direct marking when circular representative reduction
already makes the frontier cheap. The current 126-line implementation remains
unchanged pending a simpler circular-aware rule and explicit approval.

A serialized per-dimension framework was not pursued. The reusable structural
state is cheap to construct, while forbidden supports are discovered separately
for each matrix; loading an exponential table from disk would add I/O without
removing the matrix-specific work.

## Pure DFS Non-Circular Follow-up

After circular games received a separate bracelet experiment, production,
byte-indexed Gosper, and pure fixed-cardinality DFS were compared again on all
24 active non-circular matrices. Complete candidate output was byte-identical
for every variant and matrix. Each timing targeted three seconds on performance
CPU 2; dimensions 15-24 use the mean of two independent run medians. Times are
microseconds.

| ID | n | Production | Byte Gosper | Pure DFS |
|---:|---:|---:|---:|---:|
| 1 | 2 | 1.427 | 1.397 | 1.450 |
| 2 | 2 | 1.532 | 1.484 | 1.563 |
| 3 | 3 | 2.445 | 2.400 | 2.454 |
| 4 | 3 | 3.055 | 3.022 | 3.067 |
| 5 | 3 | 3.216 | 3.195 | 3.270 |
| 6 | 3 | 3.317 | 3.275 | 3.342 |
| 7 | 4 | 6.511 | 6.454 | 6.516 |
| 9 | 5 | 9.740 | 9.599 | 9.826 |
| 36 | 4 | 9.558 | 9.440 | 9.486 |
| 37 | 4 | 9.199 | 9.189 | 9.236 |
| 38 | 2 | 2.885 | 2.887 | 2.927 |
| 39 | 2 | 4.410 | 4.333 | 4.370 |
| 40 | 1 | 0.433 | 0.406 | 0.437 |
| 41 | 2 | 1.429 | 1.397 | 1.443 |
| 42 | 2 | 1.246 | 1.205 | 1.247 |
| 43 | 2 | 1.314 | 1.283 | 1.338 |
| 44 | 2 | 1.425 | 1.395 | 1.441 |
| 48 | 15 | 87,256.697 | 87,586.856 | 87,439.174 |
| 49 | 16 | 383.216 | 399.106 | 227.320 |
| 50 | 18 | 98,170.715 | 98,224.456 | 95,800.527 |
| 51 | 20 | 4,257.996 | 1,436.457 | 1,384.051 |
| 52 | 21 | 812,544.532 | 821,446.175 | 834,242.196 |
| 53 | 23 | 66,414.098 | 70,874.092 | 44,370.270 |
| 54 | 24 | 291,950.463 | 303,039.544 | 249,544.128 |

Pure DFS exposes useful non-circular cases that the adaptive dimension/count
gate hid: it improves ID 49 by 40.68%, ID 51 by 67.50%, ID 53 by 33.19%, and
ID 54 by 14.52%. It regresses ID 52 by 2.67%, so it is still not an
unconditional replacement. The machine-readable table is
`experiments/support_frontier_2026-07-29/results/non_circular_current_gosper_dfs_final.tsv`.

## Unique Full-Support ESS Control

The dense non-circular 20-by-20 game

```text
A_ij = i+j                         when i != j
A_ii = -(1 + sum_{j != i} A_ij)   = -(18i+211)
```

uses one-based indices. It has no zero entries, and its off-diagonal entries
span 37 distinct integer values from 3 through 39. This is `A=-I-L`, where `L`
is the Laplacian of the complete graph with positive edge weights `i+j`.
Consequently,

```text
z^T A z = -sum_i z_i^2 - sum_{i<j} (i+j)(z_i-z_j)^2 < 0
```

for every nonzero `z`; hence the simplex quadratic is strictly concave. Every
row sums to `-1`, so the uniform vector is an interior equilibrium. Strict
concavity makes it the only candidate, and negative definiteness makes it ESS.

The current Release executable confirmed one candidate with support mask
`1048575 = 2^20 - 1`, support size 20, `is_ess=1`, and `T_pd_frc`. The matrix,
output, and one-shot timing are saved as
`experiments/support_frontier_2026-07-29/results/n20_unique_full_support_ess_*`.

This game is active verification matrix 55. Its expected uniform candidate was
added to `baseline_candidates.csv`, and its focused CTest passed.

### Frontier timing

Production, byte-indexed Gosper, and pure fixed-cardinality DFS were measured
only on matrix 55. Each run targeted three seconds on pinned CPU 2. The table
reports two independent medians and their mean; negative change is faster.

| Variant | Run 1 | Run 2 | Mean | Change vs production |
|---|---:|---:|---:|---:|
| Production | 430.289 ms | 423.727 ms | 427.008 ms | baseline |
| Byte Gosper | 427.038 ms | 426.428 ms | 426.733 ms | -0.064% |
| Pure DFS | 425.940 ms | 426.110 ms | 426.025 ms | -0.230% |

All three candidate outputs were byte-identical. Neither alternative has a
meaningful speed advantage on this case: the only candidate has cardinality 20,
so no ESS superset pruning is available while all proper supports are examined.
Raw measurements are under
`experiments/support_frontier_2026-07-29/results/matrix_55_frontier_3s/`.

## Table-Free Gosper Candidate-List Scan

The byte-indexed Gosper experiment still allocated one byte for every one of
the `2^n` masks. This follow-up removed that table entirely. It keeps only the
current fixed-cardinality layer and the exact-candidate support masks. The
candidate masks are ordered by descending lowest set bit; while Gosper visits a
mask, the first matching candidate supplies the largest valid forbidden-block
jump. Thus its auxiliary pruning state is proportional to the candidate list,
not `2^n`.

The experimental header is
`experiments/support_frontier_2026-07-29/sources/gosper_candidate_scan/include/fracessa/supports.hpp`.
Production source was not changed. Complete CLI candidate output was
byte-identical to production for all 25 active non-circular matrices, including
candidate IDs, vectors, supports, stability labels, and payoffs.
The four existing `Supports` unit tests passed against the experimental header,
and the exhaustive/random jump-sequence check passed under ASan and UBSan.

Each timing targeted three seconds on pinned performance CPU 2 and alternated
binary order. Dimensions 15 and above report the mean of two independent run
medians, with the repeat run using the reverse binary order. Negative change is
faster; times are microseconds.

| ID | n | Candidates | Production | Candidate scan | Change |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 1 | 1.426 | 1.379 | -3.31% |
| 2 | 2 | 2 | 1.548 | 1.524 | -1.55% |
| 3 | 3 | 1 | 2.431 | 2.396 | -1.41% |
| 4 | 3 | 2 | 3.075 | 3.020 | -1.81% |
| 5 | 3 | 2 | 3.226 | 3.223 | -0.07% |
| 6 | 3 | 2 | 3.311 | 3.311 | -0.01% |
| 7 | 4 | 3 | 6.527 | 6.454 | -1.12% |
| 9 | 5 | 6 | 9.764 | 9.780 | +0.16% |
| 36 | 4 | 2 | 9.481 | 9.606 | +1.32% |
| 37 | 4 | 2 | 9.224 | 9.160 | -0.69% |
| 38 | 2 | 1 | 2.907 | 2.870 | -1.27% |
| 39 | 2 | 1 | 4.440 | 4.392 | -1.07% |
| 40 | 1 | 1 | 0.439 | 0.403 | -8.30% |
| 41 | 2 | 1 | 1.457 | 1.411 | -3.15% |
| 42 | 2 | 1 | 1.266 | 1.226 | -3.11% |
| 43 | 2 | 2 | 1.326 | 1.334 | +0.60% |
| 44 | 2 | 1 | 1.453 | 1.404 | -3.37% |
| 48 | 15 | 1 | 88,590.192 | 88,798.179 | +0.23% |
| 49 | 16 | 8 | 388.538 | 240.037 | -38.22% |
| 50 | 18 | 68 | 99,149.018 | 100,112.463 | +0.97% |
| 51 | 20 | 20 | 4,345.959 | 1,451.370 | -66.60% |
| 52 | 21 | 1 | 823,920.896 | 830,601.216 | +0.81% |
| 53 | 23 | 12 | 67,325.344 | 45,911.895 | -31.81% |
| 54 | 24 | 23 | 295,925.674 | 259,021.275 | -12.47% |
| 55 | 20 | 1 | 426,528.905 | 427,318.300 | +0.19% |

The candidate scan won 18 of 25 rows. Across all rows its geometric-mean
change was `-9.01%`, its median change was `-1.27%`, and the sum of the row
medians changed by `-2.92%`. For the eight dimensions 15-24, the geometric
mean was `-22.82%`; four won and four were effectively neutral or about 1%
slower.

The result is structurally useful but not universal. IDs 49, 51, 53, and 54
benefit because early exact candidates remove large contiguous Gosper ranges.
IDs 48, 52, and 55 have only one late candidate or are dominated by exact
matrix work, so scanning cannot help. ID 50 has 68 candidates and pays list-scan
overhead, but its measured regression remained below 1%. The implementation
therefore removes the exponential jump table and captures the important
non-circular wins without requiring the recursive DFS generator.

Machine-readable and raw results are retained at:

- `experiments/support_frontier_2026-07-29/results/non_circular_gosper_candidate_scan_final.tsv`
- `experiments/support_frontier_2026-07-29/results/non_circular_gosper_candidate_scan_raw/`
- `experiments/support_frontier_2026-07-29/results/gosper_candidate_scan_differential.tsv`

## Streaming Get-Next-Support Follow-up

The next experiment removed even the current-cardinality support vector. The
analyzer keeps the existing outer cardinality loop because the solver and
circular-symmetry rule need that size, but it now asks the generator for one
support at a time:

```cpp
supports_.start_cardinality(i);
while (const bitset64 support = supports_.get_next_support()) {
    search_one_support(support, i, is_cs_and_coprime);
}
```

`get_next_support()` retains only the next Gosper mask and the exact-candidate
list used for block jumps. It returns zero at the end of a cardinality because
zero is never a valid support. Thus the support frontier itself is constant
space; candidate masks, matrix storage, and solver state still remain.

The implementation is isolated under
`experiments/support_frontier_2026-07-29/sources/gosper_candidate_stream/`.
No production source changed. Complete `--candidates` output was byte-identical
to production for all 52 active matrices, including circular games, candidate
IDs, vectors, supports, stability labels, and payoffs.

### Runtime

The same pinned three-second protocol compared production, the table-free
one-layer candidate scan, and streaming. Dimensions 15 and above are means of
two independent medians. Negative percentages are faster.

| ID | n | Production | One layer | Stream | Stream vs production | Stream vs layer |
|---:|---:|---:|---:|---:|---:|---:|
| 48 | 15 | 88.590 ms | 88.798 ms | 88.165 ms | -0.48% | -0.71% |
| 49 | 16 | 0.389 ms | 0.240 ms | 0.240 ms | -38.20% | +0.03% |
| 50 | 18 | 99.149 ms | 100.112 ms | 101.698 ms | +2.57% | +1.58% |
| 51 | 20 | 4.346 ms | 1.451 ms | 1.523 ms | -64.95% | +4.94% |
| 52 | 21 | 823.921 ms | 830.601 ms | 827.546 ms | +0.44% | -0.37% |
| 53 | 23 | 67.325 ms | 45.912 ms | 46.320 ms | -31.20% | +0.89% |
| 54 | 24 | 295.926 ms | 259.021 ms | 262.611 ms | -11.26% | +1.39% |
| 55 | 20 | 426.529 ms | 427.318 ms | 429.466 ms | +0.69% | +0.50% |

Across all 25 non-circular matrices, streaming changed summed medians by
`-2.69%` versus production but by `+0.24%` versus the one-layer version. For
the eight dimensions 15-24, streaming was `1.02%` slower than one layer by
geometric mean and won 2 of 8 comparisons. It therefore preserves the Gosper
candidate-scan gains but is not a speed improvement by itself. Alternating
support generation with solving also loses some locality compared with first
filling one compact layer and then solving it.

### Memory

Peak resident memory was measured separately because ordinary pruning can hide
the retained layer:

| Case | Production | One layer | Stream | Stream vs layer |
|---|---:|---:|---:|---:|
| Active ID 54, strong pruning | 143,120 KiB | 13,072 KiB | 12,608 KiB | -3.55% |
| Constructed n=24, no pruning | 144,144 KiB | 34,576 KiB | 13,632 KiB | -60.57% |

The no-pruning control is the weighted-Laplacian family generalized from matrix
55. Its only candidate is the full-support ESS, so every proper support is
generated and no early candidate can keep the layer small. Elapsed time was
`9.37 s` with one layer and `9.39 s` with streaming, while streaming saved
20,944 KiB of peak resident memory.

Conclusion: the streaming interface is correct, simple, and a useful low-memory
base for a future bracelet generator. It should not be presented as a runtime
optimization; its measured benefit is removal of cardinality-layer storage.

Machine-readable results are retained at:

- `experiments/support_frontier_2026-07-29/results/non_circular_gosper_candidate_stream_final.tsv`
- `experiments/support_frontier_2026-07-29/results/gosper_candidate_stream_all_differential.tsv`
- `experiments/support_frontier_2026-07-29/results/gosper_candidate_stream_memory_id54.tsv`
- `experiments/support_frontier_2026-07-29/results/gosper_candidate_stream_memory_n24_no_pruning.tsv`
