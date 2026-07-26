# Exact rational PD experiment

> Frozen benchmark snapshot from 2026-07-26. It records that experiment, not
> the current worktree.

## Result

- Prior internal timing sum: `25.855308 s`.
- Exact-PD-only timing sum: `26.002413 s`.
- Difference: `+0.147105 s` (`+0.569%`).
- ESS counts matched for all 35 matrices.
- Candidate-level verification passed 35/35.
- Candidate-row counts were unchanged for all 35 matrices.
- The two locally constructed floating-PD false positives return zero ESS with this variant.

## Isolated change

`check_stability()` no longer constructs or tests `Bee_dbl`. It constructs `Bee` in exact rationals and calls only `matrix_frc::is_positive_definite()` before the existing exact copositivity path.

- Source: `experiments/speed_comparison_2026-07-26/sources/exact_pd_only/`
- Executable: `experiments/speed_comparison_2026-07-26/builds/exact_pd_only/fracessa`
- Source diff: `experiments/speed_comparison_2026-07-26/metadata/exact_pd_only.patch`
- Prior raw result: `experiments/speed_comparison_2026-07-26/results/current_worktree.json`
- Exact-only raw result: `experiments/speed_comparison_2026-07-26/results/exact_pd_only.json`

Both builds use GCC 16.1.1 with `-O3`, LTO, `-funroll-loops`, `-march=native`, and `-mcpu=native`. Benchmarks use internal C++ timing, candidates enabled, safe parsing, CPU 2 affinity, and 20/5/1 repetitions for dimensions below 10/below 15/15 or greater.

## Focused small-matrix benchmark

IDs 13-22 were rerun with 5 warmups and 200 alternating repetitions per binary. Alternating the prior and exact-only executables reduces frequency and run-order bias. The table reports medians.

| ID | n | Prior (s) | Exact PD (s) | Delta | Change |
|---:|---:|---:|---:|---:|---:|
| 13 | 8 | 0.000115 | 0.000133 | +18.0 us | +15.7% |
| 14 | 9 | 0.000120 | 0.000135 | +15.5 us | +13.0% |
| 15 | 10 | 0.000353 | 0.000398 | +44.0 us | +12.4% |
| 16 | 11 | 0.000121 | 0.000133 | +12.0 us | +9.9% |
| 17 | 11 | 0.000078 | 0.000082 | +4.0 us | +5.1% |
| 18 | 12 | 0.000946 | 0.001061 | +114.5 us | +12.1% |
| 19 | 12 | 0.001289 | 0.001351 | +62.0 us | +4.8% |
| 20 | 12 | 0.000248 | 0.000265 | +17.0 us | +6.9% |
| 21 | 12 | 0.002420 | 0.002562 | +142.0 us | +5.9% |
| 22 | 13 | 0.000300 | 0.000325 | +24.5 us | +8.2% |

The sum of these ten medians increased from `0.0059905 s` to `0.0064440 s`, or `7.57%`. Raw samples are in `experiments/speed_comparison_2026-07-26/results/paired_small_exact_pd.json`. Matrices below ID 13 execute in roughly 8-30 microseconds, where the executable's microsecond timing precision makes percentage comparisons unstable.

## Per-matrix timings

| ID | n | Prior (s) | Exact PD (s) | Delta (s) | Change |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.000008 | 0.000009 | +0.000001 | +12.50% |
| 2 | 2 | 0.000010 | 0.000009 | -0.000001 | -10.00% |
| 3 | 3 | 0.000010 | 0.000011 | +0.000001 | +10.00% |
| 4 | 3 | 0.000014 | 0.000014 | +0.000000 | +0.00% |
| 5 | 3 | 0.000012 | 0.000013 | +0.000001 | +8.33% |
| 6 | 3 | 0.000012 | 0.000014 | +0.000002 | +16.67% |
| 7 | 4 | 0.000021 | 0.000018 | -0.000003 | -14.29% |
| 8 | 5 | 0.000013 | 0.000013 | +0.000000 | +0.00% |
| 9 | 5 | 0.000021 | 0.000018 | -0.000003 | -14.29% |
| 10 | 6 | 0.000021 | 0.000021 | +0.000000 | +0.00% |
| 11 | 7 | 0.000019 | 0.000018 | -0.000001 | -5.26% |
| 12 | 7 | 0.000021 | 0.000018 | -0.000003 | -14.29% |
| 13 | 8 | 0.000089 | 0.000105 | +0.000016 | +17.98% |
| 14 | 9 | 0.000092 | 0.000097 | +0.000005 | +5.43% |
| 15 | 10 | 0.000268 | 0.000363 | +0.000095 | +35.45% |
| 16 | 11 | 0.000088 | 0.000100 | +0.000012 | +13.64% |
| 17 | 11 | 0.000061 | 0.000079 | +0.000018 | +29.51% |
| 18 | 12 | 0.000765 | 0.001006 | +0.000241 | +31.50% |
| 19 | 12 | 0.001104 | 0.001203 | +0.000099 | +8.97% |
| 20 | 12 | 0.000200 | 0.000246 | +0.000046 | +23.00% |
| 21 | 12 | 0.002296 | 0.002356 | +0.000060 | +2.61% |
| 22 | 13 | 0.000222 | 0.000278 | +0.000056 | +25.23% |
| 23 | 14 | 0.002101 | 0.001597 | -0.000504 | -23.99% |
| 24 | 15 | 0.004762 | 0.003830 | -0.000932 | -19.57% |
| 25 | 16 | 0.013577 | 0.012065 | -0.001512 | -11.14% |
| 26 | 17 | 0.007014 | 0.006981 | -0.000033 | -0.47% |
| 27 | 18 | 0.066413 | 0.075377 | +0.008964 | +13.50% |
| 28 | 19 | 0.014944 | 0.013845 | -0.001099 | -7.35% |
| 29 | 19 | 0.012586 | 0.013085 | +0.000499 | +3.96% |
| 30 | 20 | 0.328231 | 0.329471 | +0.001240 | +0.38% |
| 31 | 21 | 0.656870 | 0.696646 | +0.039776 | +6.06% |
| 32 | 22 | 3.325374 | 3.323353 | -0.002021 | -0.06% |
| 33 | 23 | 0.273579 | 0.288519 | +0.014940 | +5.46% |
| 34 | 24 | 21.039310 | 21.107542 | +0.068232 | +0.32% |
| 35 | 18 | 0.105180 | 0.124093 | +0.018913 | +17.98% |

For very short matrices, percentage changes amplify microsecond-level noise. Dimensions 15 and above were sampled once under the existing protocol, so small differences there should also be treated as run-to-run noise rather than exact costs.
