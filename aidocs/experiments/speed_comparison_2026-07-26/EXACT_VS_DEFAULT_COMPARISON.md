# Default versus `--exact` persistent-process benchmark

> Frozen benchmark snapshot from 2026-07-27. It measures the then-current
> worktree on the Apple M1 Max laptop. Matrix 34 is excluded by request.

## Result

- Across IDs 1-33 and 35, summed primary medians increase from `5.128579 s`
  in default mode to `171.323374 s` in `--exact` mode: `33.41x` slower.
- The independent second pass gives `4.933269 s` versus `160.354808 s`:
  `32.50x` slower. The central conclusion is therefore reproducible.
- IDs 1-12 together are only `1.31x` slower, IDs 13-22 are `4.81x`
  slower, and IDs 23-33 are `31.34x` slower in the primary pass.
- Matrix 33 is the largest primary ratio: `0.313476 s` default versus
  `73.177706 s` exact, or `233.44x` slower.
- Matrix 35 is `0.124641 s` default versus `14.594692 s` exact in the primary
  measurement. Its exact run is a single sample and measured `9.581411 s` in
  the second pass, so its useful conclusion is the broad `81x-117x` range,
  not either exact decimal value.
- Default and exact ESS counts match each other and the expected values for all
  34 included matrices. The detailed second pass also has identical candidate
  counts in both modes for every matrix.

This is not the earlier "exact-PD-only" experiment. The earlier experiment
kept the double candidate filter and changed only the final positive-definite
stability decision, costing roughly 10% on the small matrices. Full `--exact`
skips the double candidate filter, so many supports that default mode rejects
cheaply must instead run the rational candidate solver. The growing ratios are
the direct cost of that additional exact work across an exponential support
space.

## Method

- The frozen 35-matrix experiment input is used, except matrix 34.
- Each process parses its matrix once.
- Every timed repetition constructs and destroys one `fracessa` analyzer.
- Candidate construction is enabled; serialization and terminal output are
  excluded.
- Fast cases receive a 50 ms warmup, adaptive calibration, and up to seven
  measured batches targeting approximately three seconds per mode.
- A single analyzer run is never interrupted when it exceeds three seconds.
- Runs are pinned to CPU 2, an Apple M1 Max performance core.
- Default/exact execution order alternates by matrix ID.
- Both modes use one current binary compiled by GCC 16.1.1 with
  `-O3 -DNDEBUG -flto -funroll-loops -march=native -mcpu=native`.
- The only mode difference is the `exact` constructor argument.

## Reproducibility

The median absolute pass-to-pass timing difference is `3.71%` for default mode
and `2.58%` for exact mode. Most measurements are therefore stable enough for
relative conclusions. The largest percentage differences among tiny default
runs represent very small absolute times. Exact IDs 30, 31, 32, 33, and 35
exceed the target in one invocation and consequently have one exact sample in
the detailed second pass.

This frozen input contains only historical IDs 1-35, with ID 34 excluded. It
does not contain later correctness regressions 38 and 39, where the default
double candidate filter drops a valid ESS support. Agreement here therefore
does not certify the default path; `--exact` remains the correctness reference.

## Grouped results

| IDs | Primary default | Primary exact | Primary factor | Second factor |
|---:|---:|---:|---:|---:|
| 1-12 | 0.000064 s | 0.000084 s | 1.31x | 1.31x |
| 13-22 | 0.003941 s | 0.018970 s | 4.81x | 5.02x |
| 23-33 | 4.999933 s | 156.709627 s | 31.34x | 31.33x |
| 35 | 0.124641 s | 14.594692 s | 117.09x | 81.36x |
| All except 34 | 5.128579 s | 171.323374 s | 33.41x | 32.50x |

## Per-matrix primary medians

Times are milliseconds per analyzer construction.

| ID | n | Default (ms) | Exact (ms) | Exact/default |
|---:|---:|---:|---:|---:|
| 1 | 2 | 0.001322 | 0.001544 | 1.17x |
| 2 | 2 | 0.001382 | 0.001280 | 0.93x |
| 3 | 3 | 0.002098 | 0.003779 | 1.80x |
| 4 | 3 | 0.002879 | 0.003156 | 1.10x |
| 5 | 3 | 0.002908 | 0.003376 | 1.16x |
| 6 | 3 | 0.003048 | 0.003892 | 1.28x |
| 7 | 4 | 0.005845 | 0.006039 | 1.03x |
| 8 | 5 | 0.004026 | 0.005549 | 1.38x |
| 9 | 5 | 0.008954 | 0.011837 | 1.32x |
| 10 | 6 | 0.013323 | 0.016854 | 1.27x |
| 11 | 7 | 0.009456 | 0.013609 | 1.44x |
| 12 | 7 | 0.009160 | 0.013502 | 1.47x |
| 13 | 8 | 0.084019 | 0.170576 | 2.03x |
| 14 | 9 | 0.078442 | 0.118621 | 1.51x |
| 15 | 10 | 0.252864 | 1.046686 | 4.14x |
| 16 | 11 | 0.071254 | 0.330102 | 4.63x |
| 17 | 11 | 0.040453 | 0.083941 | 2.08x |
| 18 | 12 | 0.751408 | 3.429127 | 4.56x |
| 19 | 12 | 0.650215 | 3.556940 | 5.47x |
| 20 | 12 | 0.163210 | 0.259677 | 1.59x |
| 21 | 12 | 1.660420 | 8.452854 | 5.09x |
| 22 | 13 | 0.188481 | 1.521428 | 8.07x |
| 23 | 14 | 1.601033 | 26.874260 | 16.79x |
| 24 | 15 | 3.919407 | 102.980248 | 26.27x |
| 25 | 16 | 12.270484 | 177.212994 | 14.44x |
| 26 | 17 | 6.623400 | 274.112113 | 41.39x |
| 27 | 18 | 62.139109 | 1776.944466 | 28.60x |
| 28 | 19 | 10.638315 | 289.801149 | 27.24x |
| 29 | 19 | 11.059092 | 245.926186 | 22.24x |
| 30 | 20 | 339.141882 | 8021.194004 | 23.65x |
| 31 | 21 | 728.249164 | 18890.924061 | 25.94x |
| 32 | 22 | 3510.814947 | 53725.951373 | 15.30x |
| 33 | 23 | 313.475683 | 73177.706448 | 233.44x |
| 35 | 18 | 124.640873 | 14594.692096 | 117.09x |

## Files

- Primary summary: `experiments/speed_comparison_2026-07-26/results/exact_vs_default_3s_2026-07-27.csv`
- Primary JSON: `experiments/speed_comparison_2026-07-26/results/exact_vs_default_3s_2026-07-27.json`
- Detailed second pass: `experiments/speed_comparison_2026-07-26/results/exact_vs_default_3s_2026-07-27_second_pass.json`
- Harness basis: `experiments/speed_comparison_2026-07-26/microbench/benchmark_one.cpp`
