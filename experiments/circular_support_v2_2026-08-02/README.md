# Circular Support Generator V1 Versus V2

Date: 2026-08-02

## Question

Compare the then-production `CircularSupportGenerator` (V1), which expands each forbidden candidate into its distinct
rotation/reflection masks once, with the experimental `CircularSupportGeneratorV2`, which stores one representative and tests all
alignments with bit rotations during every pruning query.

## Correctness

V1 and V2 produced identical results for all 81 matrices tagged `quick_test`: candidates, ESS classifications, multipliers, fallback results, and output order all matched.

The generator logic was also checked independently with:

- every arbitrary nonempty forbidden support through dimension 10;
- every generated representative through dimension 10 and sampled representatives through dimension 18;
- explicit dihedral-orbit multiplier calculation;
- the dimension-63 boundary;
- Release and ASan/UBSan builds.

All checks passed.

## Benchmark Method

- Current fast search, changing only `CircularSupportGenerator` to `CircularSupportGeneratorV2` in the V2 build.
- GCC 16.1.1, CMake Release, `-O3 -DNDEBUG`, `-funroll-loops`, `-march=native`, and Release IPO/LTO.
- CPU 2 only.
- Both Python/native processes stayed open for the complete run; process startup was excluded.
- Measurements used the native C++ `elapsed_ns` returned by each call.
- Runs below 0.5 seconds were repeated to target approximately 0.5 seconds, and the median individual native duration was retained. Longer runs were measured once.
- Variant order alternated by matrix.
- All 81 quick-test matrices were used for correctness. Only the 33 circular matrices were timed because non-circular search calls neither V1 nor V2.

## Results

`Change` is `100 * (V2 / V1 - 1)`, so a positive value means V2 is slower.

| Matrix ID | Dimension | V1 (microseconds) | V2 (microseconds) | Change | V2 / V1 |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 0.583 | 0.583 | +0.00% | 1.000x |
| 67 | 2 | 0.500 | 0.459 | -8.20% | 0.918x |
| 68 | 3 | 0.666 | 0.625 | -6.16% | 0.938x |
| 69 | 4 | 1.125 | 1.125 | +0.00% | 1.000x |
| 8 | 5 | 1.708 | 1.709 | +0.06% | 1.001x |
| 10 | 6 | 2.250 | 2.250 | +0.00% | 1.000x |
| 11 | 7 | 3.416 | 3.750 | +9.78% | 1.098x |
| 13 | 8 | 5.667 | 6.666 | +17.63% | 1.176x |
| 14 | 9 | 6.375 | 8.500 | +33.33% | 1.333x |
| 80 | 9 | 8.750 | 10.416 | +19.04% | 1.190x |
| 15 | 10 | 11.167 | 15.500 | +38.80% | 1.388x |
| 81 | 10 | 15.709 | 19.500 | +24.13% | 1.241x |
| 16 | 11 | 22.708 | 40.042 | +76.33% | 1.763x |
| 17 | 11 | 12.083 | 16.667 | +37.94% | 1.379x |
| 18 | 12 | 25.667 | 51.000 | +98.70% | 1.987x |
| 19 | 12 | 73.541 | 106.584 | +44.93% | 1.449x |
| 20 | 12 | 12.959 | 14.708 | +13.50% | 1.135x |
| 22 | 13 | 58.583 | 123.417 | +110.67% | 2.107x |
| 23 | 14 | 114.750 | 319.667 | +178.58% | 2.786x |
| 24 | 15 | 368.375 | 1,001.708 | +171.93% | 2.719x |
| 25 | 16 | 761.021 | 3,455.188 | +354.02% | 4.540x |
| 26 | 17 | 2,180.459 | 5,159.875 | +136.64% | 2.366x |
| 27 | 18 | 3,900.000 | 17,543.251 | +349.83% | 4.498x |
| 35 | 18 | 3,316.625 | 4,692.417 | +41.48% | 1.415x |
| 28 | 19 | 9,234.917 | 49,333.546 | +434.21% | 5.342x |
| 29 | 19 | 12,993.543 | 16,164.209 | +24.40% | 1.244x |
| 30 | 20 | 40,485.191 | 168,566.600 | +316.37% | 4.164x |
| 31 | 21 | 78,632.466 | 489,680.484 | +522.75% | 6.227x |
| 32 | 22 | 194,438.352 | 1,253,605.453 | +544.73% | 6.447x |
| 33 | 23 | 168,699.453 | 595,805.557 | +253.18% | 3.532x |
| 34 | 24 | 1,585,049.069 | 10,816,588.791 | +582.41% | 6.824x |
| 70 | 25 | 91,493.863 | 259,561.775 | +183.69% | 2.837x |
| 2203 | 50 | 130.375 | 130.959 | +0.45% | 1.004x |

Summary across the 33 circular quick-test matrices:

- median V2 slowdown: **41.48%** (`1.415x`);
- geometric-mean V2 slowdown: **91.40%** (`1.914x`);
- arithmetic-mean V2 slowdown: **139.55%**;
- V1 faster on 28 matrices, equal on 3, V2 faster only on two sub-microsecond dimension-2/3 cases.

## Decision

At the time, the result was to keep V1 in production. V1 expands each orbit once and then performs cheap subset tests in the
recursive hot path. V2 saves a modest number of stored `uint64_t` masks but repeatedly scans representatives and performs rotations
for every pruning query. The compact representation is correct, but its speed tradeoff is decisively worse for this workload.

V2 never entered production. V3 later superseded V1 and is also faster than V2. The V2 source remains deliberately preserved as a
record of this failed experiment so the same design is not tested again without new evidence.
