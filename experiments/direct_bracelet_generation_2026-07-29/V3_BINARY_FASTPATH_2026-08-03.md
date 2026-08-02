# V3 binary reverse-comparison check

Date: 2026-08-03

## Changes

- Corrected `median()` for an even number of benchmark batches by averaging the two middle samples.
- Specialized `V3BraceletGenerator::update_reverse_result()` for binary words. The latest recursively placed symbol is always one, while `prefix_density_[mirrored_position]` already identifies whether the mirrored symbol is zero or one. This removes two repeated bit-mask tests without changing the algorithm.

The speed control is the previous V3 with only the median correction, so the timing comparison isolates the reverse-comparison change.

## Correctness

- Release comparison through dimension 30 passed: V3, V1/FKM, and direct BraceletFC produced identical canonical bracelets for every fixed density. Dimension 30 contained 17,920,859 nonempty bracelets.
- ASan and UBSan comparison through dimension 28 passed.
- `git diff --check` passed.

## Timing

The paired benchmark used the same release flags from `build.sh` and a three-second native timing target per implementation and dimension.

| Dimensions | Median runtime change | Mean runtime change | Geometric-mean speedup | Wins |
|---|---:|---:|---:|---:|
| 8-30 measured set | -8.97% | -9.34% | 1.106x | 17 of 18 |
| 19-24 and 30 | -8.43% | -8.49% | 1.093x | 7 of 7 |

At dimension 30, the first paired run measured 507.360 ms for the control and 466.330 ms for the optimized V3, an 8.09% reduction. A separate core-2 confirmation measured 508.351 ms and 472.113 ms respectively, a 7.13% reduction.

Hardware counters at dimension 30 showed:

| Counter | Control | Optimized | Reduction |
|---|---:|---:|---:|
| Cycles | 5,974,609,211 | 5,557,460,475 | 6.98% |
| Instructions | 34,017,817,612 | 30,550,007,956 | 10.19% |
| Branches | 4,190,403,161 | 3,701,619,880 | 11.66% |
| Branch misses | 39,034,211 | 37,415,795 | 4.15% |

Raw paired, confirmation, and `perf stat` outputs are under `results/v3_binary_fastpath/`.
