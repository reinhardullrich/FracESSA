# Exact Integer Solver Comparison

Date: 2026-07-31

## Question

Can one fraction-free symmetric factorization solve FracESSA's reduced candidate equations and provide the reduced-Hessian inertia needed by stability?

## Compared procedures

1. **Old bordered Gaussian:** the exact bordered rational solver used immediately before the reduced-Hessian change.
2. **Current reduced-Hessian LDLT:** the current working-tree implementation. It solves the candidate equations and records exact reduced-Hessian inertia.
3. **Integer bordered FFLU:** clear every game row's denominators once, build an integer bordered system, and call `fmpz_mat_solve_fflu`. This is candidate-only.
4. **Conservative integer hybrid:** use integer FFLU for every support, then call the unchanged current LDLT finder only for a successful candidate. It therefore supplies and verifies the same inertia as current LDLT.
5. **FFLDLT candidate solver:** clear the whole symmetric game's denominators once, build the integer border-reduced system $H y=r$, and solve it with one symmetric fraction-free LDLT-style factorization. The same factorization records exact inertia. It is a general symmetric-indefinite solver: pivot signs never reject a candidate.
6. **Rational fraction-free:** build the rational reduced system and call `fmpq_mat_solve_fraction_free`. This is candidate-only.

The conservative hybrid deliberately repeats the solve and outside-support checks for successful candidates. FFLDLT instead performs the candidate solve and inertia classification together. Its immediate-integer arithmetic is adapted from FLINT's FFLU implementation and switches to ordinary arbitrary-precision `fmpz` operations as soon as an intermediate result no longer fits FLINT's immediate representation.

Production source is unchanged. The isolated executable uses the production parser, support generators, pruning rule, and current exact candidate finder.

## Build and run

```bash
cmake -S experiments/exact_integer_solver_comparison_2026-07-31 \
      -B experiments/exact_integer_solver_comparison_2026-07-31/builds/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
cmake --build experiments/exact_integer_solver_comparison_2026-07-31/builds/release -j2
```

One matrix can be measured on CPU 2 with the adaptive one-second protocol:

```bash
taskset -c 2 experiments/exact_integer_solver_comparison_2026-07-31/builds/release/exact_integer_solver_benchmark \
    '3#0,1,0,0,1,0' auto 1000
```

After the cross-procedure checks, the new variant alone can be timed without repeating the five stored baselines:

```bash
taskset -c 2 experiments/exact_integer_solver_comparison_2026-07-31/builds/release/exact_integer_solver_benchmark \
    '3#0,1,0,0,1,0' ffldlt_only 1000
```

## Method

- Apple M1 Max performance CPU 2, Fedora Linux Asahi Remix.
- GCC 16.1.1, C++17 Release, `-march=native`, loop unrolling, and LTO.
- FLINT 3.4.0 and GMP 6.3.0.
- Timing starts inside C++ after process startup and parsing; process-launch overhead is excluded.
- Each timed search includes solver construction, support generation, and one-time denominator preprocessing.
- The first measured duration determines `ceil(1 second / duration)` samples for that variant; the first sample is retained and the median is reported.
- Long variants naturally have one sample.
- Stability after the stored inertia decision is excluded. Current LDLT, the conservative hybrid, and FFLDLT all compute exact reduced-Hessian inertia.
- IDs 45, 47, 65, 66, and 90 were excluded because the stored current reduced-Hessian exact time exceeds two minutes.
- IDs 33 and 34 have no special rule and are included normally.
- The raw nanosecond values and sample counts are in [all_matrices_cpu2.tsv](all_matrices_cpu2.tsv).

## Result summary

FFLDLT alone solves every nonsingular reduced candidate system, reconstructs and validates the candidate, and records inertia. Indefinite systems continue through the complete candidate test; only singularity, a nonpositive support probability, or an outside payoff violation rejects the support.

The Old Gaussian summed baseline is **458,322,358.311 microseconds**. Every percentage below uses Old Gaussian as its baseline; a negative value means the named alternative is faster.

| Comparison | Geometric-mean change | Median change | Change in summed medians | Alternative total (µs) |
|---|---:|---:|---:|---:|
| Old Gaussian vs Current LDLT | -32.52% | -29.87% | -28.82% | 326,256,174.295 |
| Old Gaussian vs Integer FFLU | -55.70% | -61.21% | -81.90% | 82,952,181.061 |
| Old Gaussian vs Conservative hybrid | -50.57% | -59.47% | -81.84% | 83,230,364.284 |
| Old Gaussian vs FFLDLT | -73.94% | -77.08% | -87.23% | 58,543,326.382 |
| Old Gaussian vs Rational FFLU | -55.16% | -53.91% | -70.72% | 134,188,201.401 |

Relative to Old Gaussian, FFLDLT reduces the geometric mean by 73.94%, the median by 77.08%, and the summed medians by 87.23%. Its 58,543,326.382-microsecond total is also below candidate-only Integer FFLU's 82,952,181.061 microseconds while FFLDLT additionally supplies exact inertia.

Summed medians are dominated by the longest matrices, so the dimension groups show how the Old Gaussian comparison changes with dimension:

| Dimension group | Old Gaussian total (µs) | FFLDLT total (µs) | Old Gaussian vs FFLDLT |
|---|---:|---:|---:|
| 2-8 | 240.708 | 79.710 | -66.89% |
| 9-16 | 988,418.058 | 210,007.863 | -78.75% |
| 17-25 | 457,333,699.545 | 58,333,238.809 | -87.24% |

The improvement over Old Gaussian grows strongly with dimension. These results support a production experiment replacing the rational reduced-system solver with the specialized FFLDLT candidate solver; production source remains unchanged here.

## Complete 82-matrix table

Times are microseconds. Every percentage uses Old Gaussian as its baseline; a negative value means the named alternative is faster. Rows are sorted by dimension, then non-circular before circular, then matrix ID.

| ID | n | CS | Visited | Candidates | Old Gaussian µs | Current LDLT µs | Old Gaussian vs Current LDLT | Integer FFLU µs | Old Gaussian vs Integer FFLU | Conservative hybrid µs | Old Gaussian vs Conservative hybrid | FFLDLT µs | Old Gaussian vs FFLDLT | Rational FFLU µs | Old Gaussian vs Rational FFLU |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | no | 3 | 1 | 0.917 | 0.625 | -31.84% | 1.333 | +45.37% | 1.625 | +77.21% | 0.584 | -36.31% | 0.625 | -31.84% |
| 2 | 2 | no | 2 | 2 | 0.833 | 0.500 | -39.98% | 1.042 | +25.09% | 1.375 | +65.07% | 0.584 | -29.89% | 0.500 | -39.98% |
| 38 | 2 | no | 3 | 1 | 1.917 | 1.167 | -39.12% | 2.125 | +10.85% | 2.750 | +43.45% | 0.917 | -52.16% | 1.125 | -41.31% |
| 39 | 2 | no | 3 | 1 | 4.750 | 1.000 | -78.95% | 2.459 | -48.23% | 2.917 | -38.59% | 0.875 | -81.58% | 1.042 | -78.06% |
| 41 | 2 | no | 3 | 1 | 1.042 | 0.584 | -43.95% | 1.292 | +23.99% | 1.583 | +51.92% | 0.709 | -31.96% | 0.625 | -40.02% |
| 42 | 2 | no | 2 | 1 | 0.625 | 0.375 | -40.00% | 0.917 | +46.72% | 1.083 | +73.28% | 0.458 | -26.72% | 0.416 | -33.44% |
| 43 | 2 | no | 2 | 2 | 0.666 | 0.458 | -31.23% | 0.959 | +43.99% | 1.250 | +87.69% | 0.500 | -24.92% | 0.458 | -31.23% |
| 44 | 2 | no | 3 | 1 | 1.042 | 0.583 | -44.05% | 1.292 | +23.99% | 1.583 | +51.92% | 0.584 | -43.95% | 0.625 | -40.02% |
| 67 | 2 | yes | 1 | 1 | 0.458 | 0.334 | -27.07% | 0.708 | +54.59% | 0.875 | +91.05% | 0.417 | -8.95% | 0.375 | -18.12% |
| 3 | 3 | no | 6 | 1 | 2.625 | 1.250 | -52.38% | 2.500 | -4.76% | 2.875 | +9.52% | 1.125 | -57.14% | 1.459 | -44.42% |
| 4 | 3 | no | 4 | 2 | 1.583 | 0.875 | -44.73% | 1.750 | +10.55% | 2.333 | +47.38% | 0.834 | -47.32% | 0.917 | -42.07% |
| 5 | 3 | no | 6 | 2 | 1.875 | 1.333 | -28.91% | 2.417 | +28.91% | 3.042 | +62.24% | 1.167 | -37.76% | 1.416 | -24.48% |
| 6 | 3 | no | 6 | 2 | 2.292 | 1.292 | -43.63% | 2.625 | +14.53% | 3.291 | +43.59% | 1.208 | -47.29% | 1.542 | -32.72% |
| 68 | 3 | yes | 1 | 1 | 0.542 | 0.459 | -15.31% | 0.834 | +53.87% | 1.042 | +92.25% | 0.625 | +15.31% | 0.458 | -15.50% |
| 7 | 4 | no | 5 | 3 | 1.916 | 1.208 | -36.95% | 2.292 | +19.62% | 3.167 | +65.29% | 1.042 | -45.62% | 1.292 | -32.57% |
| 36 | 4 | no | 8 | 2 | 4.042 | 2.125 | -47.43% | 3.542 | -12.37% | 4.458 | +10.29% | 1.625 | -59.80% | 2.334 | -42.26% |
| 37 | 4 | no | 8 | 2 | 3.917 | 2.042 | -47.87% | 3.542 | -9.57% | 4.417 | +12.76% | 1.667 | -57.44% | 2.292 | -41.49% |
| 69 | 4 | yes | 3 | 1 | 1.500 | 1.000 | -33.33% | 1.875 | +25.00% | 2.291 | +52.73% | 1.125 | -25.00% | 1.208 | -19.47% |
| 9 | 5 | no | 16 | 6 | 8.083 | 4.792 | -40.72% | 7.167 | -11.33% | 9.708 | +20.10% | 3.750 | -53.61% | 5.459 | -32.46% |
| 8 | 5 | yes | 5 | 1 | 3.208 | 2.167 | -32.45% | 3.291 | +2.59% | 3.959 | +23.41% | 2.000 | -37.66% | 2.333 | -27.28% |
| 71 | 6 | no | 37 | 2 | 26.708 | 13.375 | -49.92% | 18.125 | -32.14% | 19.208 | -28.08% | 8.458 | -68.33% | 13.834 | -48.20% |
| 10 | 6 | yes | 5 | 2 | 3.459 | 2.541 | -26.54% | 3.667 | +6.01% | 4.667 | +34.92% | 2.416 | -30.15% | 2.750 | -20.50% |
| 72 | 7 | no | 23 | 4 | 14.625 | 7.917 | -45.87% | 11.541 | -21.09% | 13.500 | -7.69% | 5.208 | -64.39% | 8.583 | -41.31% |
| 11 | 7 | yes | 9 | 2 | 7.458 | 5.209 | -30.16% | 6.542 | -12.28% | 8.375 | +12.30% | 4.208 | -43.58% | 5.375 | -27.93% |
| 73 | 8 | no | 90 | 3 | 119.209 | 55.334 | -53.58% | 55.000 | -53.86% | 56.833 | -52.32% | 27.916 | -76.58% | 48.625 | -59.21% |
| 13 | 8 | yes | 20 | 3 | 25.416 | 16.500 | -35.08% | 15.833 | -37.70% | 20.125 | -20.82% | 9.708 | -61.80% | 15.291 | -39.84% |
| 74 | 9 | no | 257 | 4 | 498.124 | 283.291 | -43.13% | 215.208 | -56.80% | 217.333 | -56.37% | 101.083 | -79.71% | 206.042 | -58.64% |
| 14 | 9 | yes | 16 | 3 | 16.792 | 11.666 | -30.53% | 13.000 | -22.58% | 16.292 | -2.98% | 8.708 | -48.14% | 11.667 | -30.52% |
| 80 | 9 | yes | 35 | 3 | 64.667 | 40.375 | -37.56% | 31.042 | -52.00% | 37.833 | -41.50% | 19.916 | -69.20% | 36.041 | -44.27% |
| 75 | 10 | no | 346 | 5 | 793.000 | 420.583 | -46.96% | 319.375 | -59.73% | 319.374 | -59.73% | 145.208 | -81.69% | 296.834 | -62.57% |
| 15 | 10 | yes | 40 | 4 | 77.500 | 46.208 | -40.38% | 36.042 | -53.49% | 43.542 | -43.82% | 22.708 | -70.70% | 37.834 | -51.18% |
| 81 | 10 | yes | 66 | 4 | 227.167 | 143.125 | -37.00% | 72.625 | -68.03% | 89.250 | -60.71% | 43.458 | -80.87% | 94.000 | -58.62% |
| 46 | 11 | no | 2,047 | 1 | 25,237.860 | 28,959.567 | +14.75% | 12,594.493 | -50.10% | 12,655.242 | -49.86% | 8,176.309 | -67.60% | 18,777.156 | -25.60% |
| 76 | 11 | no | 578 | 3 | 1,747.020 | 1,020.437 | -41.59% | 680.478 | -61.05% | 681.042 | -61.02% | 307.292 | -82.41% | 646.229 | -63.01% |
| 82 | 11 | no | 884 | 70 | 2,243.082 | 1,372.416 | -38.82% | 786.333 | -64.94% | 953.978 | -57.47% | 433.083 | -80.69% | 909.041 | -59.47% |
| 16 | 11 | yes | 83 | 5 | 184.334 | 131.500 | -28.66% | 87.125 | -52.74% | 99.666 | -45.93% | 57.375 | -68.87% | 104.874 | -43.11% |
| 17 | 11 | yes | 27 | 4 | 41.625 | 31.417 | -24.52% | 27.833 | -33.13% | 34.250 | -17.72% | 18.750 | -54.95% | 28.208 | -32.23% |
| 77 | 12 | no | 1,364 | 5 | 5,033.039 | 3,432.019 | -31.81% | 1,944.436 | -61.37% | 1,948.791 | -61.28% | 906.791 | -81.98% | 2,052.686 | -59.22% |
| 18 | 12 | yes | 80 | 8 | 234.791 | 162.500 | -30.79% | 91.792 | -60.90% | 112.208 | -52.21% | 54.583 | -76.75% | 104.041 | -55.69% |
| 19 | 12 | yes | 103 | 7 | 309.958 | 177.084 | -42.87% | 104.417 | -66.31% | 127.542 | -58.85% | 64.584 | -79.16% | 127.000 | -59.03% |
| 20 | 12 | yes | 15 | 3 | 20.041 | 15.125 | -24.53% | 14.833 | -25.99% | 17.583 | -12.26% | 10.792 | -46.15% | 14.833 | -25.99% |
| 78 | 13 | no | 1,538 | 3 | 7,451.038 | 4,502.164 | -39.58% | 2,450.644 | -67.11% | 2,454.062 | -67.06% | 1,180.583 | -84.16% | 2,594.936 | -65.17% |
| 83 | 13 | no | 4,309 | 157 | 18,863.781 | 12,635.410 | -33.02% | 5,482.622 | -70.94% | 6,014.455 | -68.12% | 3,325.957 | -82.37% | 7,254.704 | -61.54% |
| 22 | 13 | yes | 191 | 7 | 807.958 | 521.020 | -35.51% | 254.125 | -68.55% | 284.750 | -64.76% | 163.875 | -79.72% | 317.750 | -60.67% |
| 79 | 14 | no | 1,769 | 6 | 6,729.580 | 4,738.039 | -29.59% | 2,578.749 | -61.68% | 2,590.165 | -61.51% | 1,243.625 | -81.52% | 2,767.728 | -58.87% |
| 84 | 14 | no | 10,903 | 233 | 91,062.948 | 58,904.216 | -35.31% | 18,185.365 | -80.03% | 19,676.739 | -78.39% | 10,324.995 | -88.66% | 25,480.840 | -72.02% |
| 23 | 14 | yes | 324 | 10 | 1,456.250 | 1,091.416 | -25.05% | 533.083 | -63.39% | 574.458 | -60.55% | 342.542 | -76.48% | 696.666 | -52.16% |
| 48 | 15 | no | 16,384 | 1 | 172,327.984 | 157,820.325 | -8.42% | 104,225.690 | -39.52% | 104,181.294 | -39.54% | 78,331.797 | -54.54% | 92,931.155 | -46.07% |
| 56 | 15 | no | 17,548 | 216 | 88,158.657 | 63,763.068 | -27.67% | 27,672.400 | -68.61% | 28,231.066 | -67.98% | 20,991.949 | -76.19% | 48,018.223 | -45.53% |
| 24 | 15 | yes | 921 | 16 | 7,338.392 | 5,619.038 | -23.43% | 2,147.749 | -70.73% | 2,293.749 | -68.74% | 1,335.750 | -81.80% | 2,885.644 | -60.68% |
| 49 | 16 | no | 263 | 8 | 240.750 | 183.666 | -23.71% | 164.083 | -31.85% | 170.833 | -29.04% | 103.000 | -57.22% | 179.625 | -25.39% |
| 57 | 16 | no | 29,844 | 324 | 164,774.363 | 129,816.196 | -21.22% | 51,715.345 | -68.61% | 52,378.490 | -68.21% | 39,134.190 | -76.25% | 93,015.821 | -43.55% |
| 85 | 16 | no | 35,318 | 536 | 381,675.238 | 255,454.186 | -33.07% | 70,428.376 | -81.55% | 73,745.000 | -80.68% | 41,211.648 | -89.20% | 97,798.693 | -74.38% |
| 25 | 16 | yes | 1,253 | 20 | 10,802.119 | 7,929.078 | -26.60% | 2,928.269 | -72.89% | 3,056.956 | -71.70% | 1,947.312 | -81.97% | 4,072.998 | -62.29% |
| 58 | 17 | no | 81,136 | 486 | 601,724.048 | 425,371.171 | -29.31% | 164,930.613 | -72.59% | 166,371.862 | -72.35% | 122,959.922 | -79.57% | 291,111.664 | -51.62% |
| 86 | 17 | no | 62,654 | 784 | 758,374.582 | 511,696.100 | -32.53% | 135,848.359 | -82.09% | 140,484.376 | -81.48% | 80,503.484 | -89.38% | 188,696.432 | -75.12% |
| 26 | 17 | yes | 3,677 | 24 | 121,354.721 | 97,142.798 | -19.95% | 21,252.863 | -82.49% | 22,556.487 | -81.41% | 10,588.224 | -91.27% | 25,053.444 | -79.36% |
| 50 | 18 | no | 139,130 | 68 | 1,242,790.657 | 1,134,865.969 | -8.68% | 506,177.249 | -59.27% | 506,944.456 | -59.21% | 290,816.426 | -76.60% | 581,502.998 | -53.21% |
| 59 | 18 | no | 154,746 | 648 | 1,319,931.237 | 952,548.741 | -27.83% | 353,303.712 | -73.23% | 355,970.128 | -73.03% | 267,002.584 | -79.77% | 628,745.408 | -52.37% |
| 87 | 18 | no | 164,748 | 1,164 | 2,157,345.544 | 1,645,840.423 | -23.71% | 456,916.527 | -78.82% | 467,424.855 | -78.33% | 294,350.549 | -86.36% | 640,168.026 | -70.33% |
| 27 | 18 | yes | 4,835 | 36 | 62,116.172 | 48,489.784 | -21.94% | 15,561.741 | -74.95% | 15,953.283 | -74.32% | 10,895.286 | -82.46% | 21,690.279 | -65.08% |
| 35 | 18 | yes | 7,591 | 11 | 391,300.148 | 284,328.919 | -27.34% | 50,016.908 | -87.22% | 51,187.804 | -86.92% | 24,557.342 | -93.72% | 67,058.336 | -82.86% |
| 60 | 19 | no | 273,154 | 972 | 2,502,754.178 | 1,922,023.306 | -23.20% | 681,754.314 | -72.76% | 685,470.208 | -72.61% | 531,205.668 | -78.78% | 1,257,747.606 | -49.75% |
| 88 | 19 | no | 304,102 | 1,694 | 4,577,673.396 | 3,502,734.225 | -23.48% | 934,614.752 | -79.58% | 948,095.556 | -79.29% | 608,580.174 | -86.71% | 1,302,683.205 | -71.54% |
| 28 | 19 | yes | 8,042 | 60 | 144,288.833 | 107,305.084 | -25.63% | 29,890.295 | -79.28% | 30,780.149 | -78.67% | 21,446.010 | -85.14% | 36,658.208 | -74.59% |
| 29 | 19 | yes | 12,776 | 4 | 120,537.056 | 92,323.405 | -23.41% | 37,888.061 | -68.57% | 37,824.437 | -68.62% | 27,552.216 | -77.14% | 50,711.804 | -57.93% |
| 51 | 20 | no | 20 | 20 | 25.083 | 21.792 | -13.12% | 22.084 | -11.96% | 43.000 | +71.43% | 9.416 | -62.46% | 21.792 | -13.12% |
| 55 | 20 | no | 1,048,575 | 1 | 102,021,095.861 | 61,490,146.193 | -39.73% | 7,645,757.804 | -92.51% | 7,659,938.551 | -92.49% | 4,279,167.371 | -95.81% | 16,316,634.585 | -84.01% |
| 61 | 20 | no | 697,086 | 1,458 | 8,290,657.874 | 5,870,352.648 | -29.19% | 2,060,713.934 | -75.14% | 2,067,285.430 | -75.06% | 1,570,267.598 | -81.06% | 3,624,898.363 | -56.28% |
| 30 | 20 | yes | 18,801 | 72 | 358,561.167 | 297,706.598 | -16.97% | 86,035.388 | -76.01% | 87,019.783 | -75.73% | 59,632.660 | -83.37% | 116,975.932 | -67.38% |
| 52 | 21 | no | 1,572,864 | 1 | 19,347,458.190 | 10,242,418.538 | -47.06% | 5,389,934.801 | -72.14% | 5,395,540.422 | -72.11% | 3,813,422.798 | -80.29% | 6,966,168.307 | -63.99% |
| 62 | 21 | no | 1,342,780 | 1,944 | 18,231,955.002 | 13,110,048.921 | -28.09% | 4,430,005.898 | -75.70% | 4,459,200.672 | -75.54% | 3,408,102.446 | -81.31% | 7,845,274.339 | -56.97% |
| 31 | 21 | yes | 38,909 | 123 | 985,539.264 | 799,231.767 | -18.90% | 223,979.829 | -77.27% | 226,428.827 | -77.02% | 159,852.967 | -83.78% | 297,316.390 | -69.83% |
| 63 | 22 | no | 2,432,484 | 2,916 | 35,214,308.184 | 26,759,302.906 | -24.01% | 8,776,647.885 | -75.08% | 8,789,783.294 | -75.04% | 6,909,950.014 | -80.38% | 15,980,903.595 | -54.62% |
| 89 | 22 | no | 3,095,388 | 6,300 | 90,000,177.444 | 71,196,569.045 | -20.89% | 15,537,235.283 | -82.74% | 15,639,881.830 | -82.62% | 9,938,538.388 | -88.96% | 22,064,749.198 | -75.48% |
| 32 | 22 | yes | 72,138 | 156 | 1,601,594.199 | 1,290,718.545 | -19.41% | 473,906.518 | -70.41% | 475,849.100 | -70.29% | 368,168.786 | -77.01% | 633,743.759 | -60.43% |
| 53 | 23 | no | 107,275 | 12 | 1,617,821.731 | 1,216,776.089 | -24.79% | 356,732.419 | -77.95% | 356,421.502 | -77.97% | 194,854.554 | -87.96% | 418,220.675 | -74.15% |
| 64 | 23 | no | 5,922,352 | 4,374 | 107,894,178.079 | 76,002,069.081 | -29.56% | 24,892,561.047 | -76.93% | 24,920,362.833 | -76.90% | 19,294,581.090 | -82.12% | 43,444,342.041 | -59.73% |
| 33 | 23 | yes | 179,397 | 95 | 34,960,738.373 | 28,070,901.897 | -19.71% | 3,775,946.317 | -89.20% | 3,800,327.344 | -89.13% | 2,055,311.913 | -94.12% | 4,068,153.515 | -88.36% |
| 54 | 24 | no | 520,742 | 23 | 8,512,539.454 | 6,566,341.454 | -22.86% | 1,949,961.456 | -77.09% | 1,943,587.919 | -77.17% | 1,174,758.283 | -86.20% | 2,415,328.395 | -71.63% |
| 34 | 24 | yes | 273,666 | 345 | 11,592,848.756 | 9,554,646.060 | -17.58% | 3,039,419.618 | -73.78% | 3,047,124.905 | -73.72% | 2,444,065.814 | -78.92% | 3,775,172.359 | -67.44% |
| 70 | 25 | yes | 110,976 | 8 | 2,704,010.312 | 2,324,902.656 | -14.02% | 619,225.184 | -77.10% | 619,320.996 | -77.10% | 372,096.826 | -86.24% | 726,888.518 | -73.12% |

## Verification

- Detailed runs on IDs 38, 39, 46, 48, and 51 compared every exact candidate vector and payoff across all six procedures.
- Those runs also matched every per-candidate negative-definiteness decision. In particular, the nonsingular indefinite case $H=\begin{bmatrix}0&1\\1&0\end{bmatrix}$ returned its candidate and recorded non-negative-definite inertia instead of rejecting it.
- Another 357 deterministic symmetric games compared complete candidates and inertia against current exact: dimensions 2 through 12, ordinary fractions, and 72 forced arbitrary-precision cases. Four structured systems forced zero-diagonal pivots both initially and after earlier elimination steps.
- Another 105 adversarial indefinite systems used randomly ordered $2\times2$ zero-diagonal blocks, forcing repeated coordinate repairs while still returning the exact candidates needed for superset pruning.
- The complete FFLDLT-only timing sweep matched every stored visited-support and candidate count. A separate sweep matched the ordered support/extended-support checksum against every retained exact database candidate, and repeated samples matched support, extended-support, and inertia checksums against their first sample.
- AddressSanitizer, UndefinedBehaviorSanitizer, and leak detection passed 48 additional deterministic integer, fractional, and arbitrary-precision games.
- `git diff --check` passes.
