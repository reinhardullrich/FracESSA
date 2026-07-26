# Persistent-process exact-PD microbenchmark

> Frozen benchmark snapshot from 2026-07-26. It records that experiment, not
> the current worktree.

## Result

- IDs 1-22: summed per-matrix medians increased from `0.003445084 s` to `0.003782342 s` (`+9.79%`).
- IDs 13-22: summed per-matrix medians increased from `0.003387055 s` to `0.003721255 s` (`+9.87%`).
- All 35 matrices: the summed medians increased from `25.610597 s` to `25.861934 s` (`+0.98%`). This aggregate is dominated by matrix 34, which takes about 21 seconds per run.
- The largest measured increase was matrix 13 (`+25.58%`); 9 of 35 matrices increased by more than 10%.
- Exact ESS counts matched the expected baseline for all 35 matrices. Candidate counts matched between both variants for all 35 matrices.

The persistent-process result is the useful result for the small matrices: removing process startup and parsing exposes an approximately 10% exact-PD cost across IDs 1-22. The all-matrix percentage hides that cost because one long-running matrix dominates the sum.

## What is measured

- Each harness process parses its matrix once.
- Every timed repetition constructs and destroys a new `fracessa` analyzer.
- Candidate construction is enabled. Candidate serialization and terminal output are excluded.
- Fast cases use a warmup batch, a calibration batch of up to 20 ms, and up to seven measured batches targeting about three seconds in total.
- The reported value is the median nanoseconds per analyzer run across measured batches.
- The runner pins execution to CPU 2, an Apple M1 Max Firestorm performance core, and alternates which variant runs first by matrix ID.
- Both variants were compiled with GCC 16.1.1 using `-O3 -DNDEBUG -flto -funroll-loops -march=native -mcpu=native`.
- The source trees differ only in `cpp/src/checkstab.cpp`: the experimental variant skips the double positive-definiteness decision and performs the exact rational check directly.

## Interpretation limits

- This measures warm repeated throughput for one matrix, not cold command-line latency. It intentionally answers whether repeating a matrix inside one process improves timing quality.
- Matrices 32 and 34 already exceed the three-second target on their first run, so each has one sample per variant. Their small percentage differences are not precise microbenchmark estimates.
- Small negative changes are measurement and code-layout noise; the exact-only path does not perform less logical work.
- The harness validates ESS and candidate counts, not serialized candidate-row equality. The separate exact-PD experiment already passed the full candidate-level verifier for all 35 matrices.

## Grouped results

| IDs | Current sum | Exact-PD sum | Weighted change | Median row change |
|---:|---:|---:|---:|---:|
| 1-12 | 0.000058029 s | 0.000061087 s | +5.27% | +6.20% |
| 13-22 | 0.003387055 s | 0.003721255 s | +9.87% | +12.39% |
| 1-22 | 0.003445084 s | 0.003782342 s | +9.79% | +7.96% |
| 23-35 | 25.607152 s | 25.858152 s | +0.98% | +4.25% |
| 1-35 | 25.610597 s | 25.861934 s | +0.98% | +5.70% |

## Per-matrix medians

| ID | n | Current (us) | Exact PD (us) | Change | Repetitions C/E | Batches C/E |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 1.128 | 1.252 | +11.06% | 2512753/2016825 | 7/7 |
| 2 | 2 | 1.318 | 1.315 | -0.28% | 1933909/1973615 | 7/7 |
| 3 | 3 | 1.877 | 2.009 | +7.07% | 1445291/1343219 | 7/7 |
| 4 | 3 | 2.649 | 2.746 | +3.67% | 1113371/993466 | 7/7 |
| 5 | 3 | 2.869 | 2.768 | -3.50% | 933329/961412 | 7/7 |
| 6 | 3 | 2.862 | 2.905 | +1.51% | 906117/920990 | 7/7 |
| 7 | 4 | 5.763 | 5.592 | -2.97% | 476959/496253 | 7/7 |
| 8 | 5 | 3.562 | 3.888 | +9.15% | 781208/711308 | 7/7 |
| 9 | 5 | 7.817 | 8.498 | +8.72% | 384097/335683 | 7/7 |
| 10 | 6 | 11.946 | 12.583 | +5.33% | 242671/227475 | 7/7 |
| 11 | 7 | 8.125 | 8.780 | +8.06% | 357757/319286 | 7/7 |
| 12 | 7 | 8.113 | 8.751 | +7.86% | 346390/341878 | 7/7 |
| 13 | 8 | 60.836 | 76.396 | +25.58% | 48476/39302 | 7/7 |
| 14 | 9 | 66.045 | 76.509 | +15.84% | 44894/39175 | 7/7 |
| 15 | 10 | 208.220 | 245.427 | +17.87% | 14320/12042 | 7/7 |
| 16 | 11 | 59.842 | 68.177 | +13.93% | 49926/43809 | 7/7 |
| 17 | 11 | 36.507 | 38.517 | +5.51% | 79063/77061 | 7/7 |
| 18 | 12 | 618.397 | 716.615 | +15.88% | 4804/4122 | 7/7 |
| 19 | 12 | 619.075 | 630.386 | +1.83% | 4732/4692 | 7/7 |
| 20 | 12 | 149.938 | 155.532 | +3.73% | 19886/19144 | 7/7 |
| 21 | 12 | 1404.756 | 1532.512 | +9.09% | 2134/1955 | 7/7 |
| 22 | 13 | 163.440 | 181.186 | +10.86% | 18087/16225 | 7/7 |
| 23 | 14 | 1503.819 | 1524.760 | +1.39% | 1990/1967 | 7/7 |
| 24 | 15 | 3647.889 | 3748.958 | +2.77% | 823/795 | 7/7 |
| 25 | 16 | 10450.305 | 11805.638 | +12.97% | 288/254 | 7/7 |
| 26 | 17 | 5811.420 | 6371.836 | +9.64% | 516/472 | 7/7 |
| 27 | 18 | 59741.604 | 59863.890 | +0.20% | 51/51 | 7/7 |
| 28 | 19 | 9583.403 | 9990.839 | +4.25% | 313/302 | 7/7 |
| 29 | 19 | 10642.363 | 10358.701 | -2.67% | 285/290 | 7/7 |
| 30 | 20 | 306596.241 | 323463.227 | +5.50% | 10/10 | 7/7 |
| 31 | 21 | 650995.409 | 694427.998 | +6.67% | 5/5 | 5/5 |
| 32 | 22 | 3340335.892 | 3313418.124 | -0.81% | 1/1 | 1/1 |
| 33 | 23 | 268598.899 | 283903.177 | +5.70% | 12/11 | 7/7 |
| 34 | 24 | 20835187.517 | 21019492.279 | +0.88% | 1/1 | 1/1 |
| 35 | 18 | 104057.337 | 119782.076 | +15.11% | 29/26 | 7/7 |

## Files

- Harness: `experiments/speed_comparison_2026-07-26/microbench/benchmark_one.cpp`
- Build script: `experiments/speed_comparison_2026-07-26/microbench/build.sh`
- 35-matrix runner: `experiments/speed_comparison_2026-07-26/microbench/run_all.py`
- Raw measurements: `experiments/speed_comparison_2026-07-26/results/microbench_exact_pd.json`
- Earlier command-line benchmark: `EXACT_PD_ONLY_COMPARISON.md`

## Reproduce

From `experiments/speed_comparison_2026-07-26`:

```bash
./microbench/build.sh
./microbench/run_all.py --target-seconds 3
```

The harness parses each matrix once, warms short cases for 50 ms, runs up to
seven adaptive batches targeting about three measured seconds, and reports the
median nanoseconds per analyzer invocation.
