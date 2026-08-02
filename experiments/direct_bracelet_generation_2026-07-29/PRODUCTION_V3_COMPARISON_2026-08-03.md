# Production V3 Comparison

Date: 2026-08-03

## Change

FracESSA now uses `CircularSupportGeneratorV3` for circular matrices. V3 is the binary specialization of the direct fixed-density
`BraceFD` recursion from Karim, Alamgir, and Husnine (2014). The previous `CircularSupportGenerator` remains as V1, and
`CircularSupportGeneratorV2` also remains available. V3 retains V1's callback, expanded forbidden-orbit pruning, and orbit
multiplier.

## Correctness

- V3 and V1 emitted identical ordered representatives through dimension 24.
- A dynamic forbidden-support run produced the same remaining sequence and multiplier.
- V3's dimension-63 boundary passed.
- All 81 quick-test matrices produced byte-identical complete candidate output in both fast and safe modes.
- All 10 Release CTests and all 10 ASan/UBSan CTests passed.
- All 63 Python wrapper tests passed.

## Runtime method

Both revisions were built immediately before and after the production switch with GCC 16.1.1 and the canonical Release/native/LTO
CMake configuration. Each build ran in one persistent Pybind process pinned to CPU 2. The panel contains the 31 dimension-3-or-
larger circular quick-test matrices and uses each matrix's retained calibration with a 0.5-second native-duration target.

The forward-order panel measured V1 first and V3 second. The reverse-order panel measured V3 first and V1 second.

| Panel | V3 median change | V3 geometric-mean change | Wins |
|---|---:|---:|---:|
| Forward order | -23.68% | -24.64% | 31/31 |
| Reverse order | -19.90% | -21.04% | 30/31 plus one tie |
| Forward, dimensions 19+ | -34.77% | -30.29% | 9/9 |
| Reverse, dimensions 19+ | -34.70% | -29.61% | 9/9 |

The complete conservative reverse-order panel is:

| ID | n | V1 (ns) | V3 (ns) | Change |
|---:|---:|---:|---:|---:|
| 68 | 3 | 666 | 666 | 0.00% |
| 69 | 4 | 1,083 | 1,042 | -3.79% |
| 8 | 5 | 1,667 | 1,583 | -5.04% |
| 10 | 6 | 2,250 | 2,084 | -7.38% |
| 11 | 7 | 3,375 | 3,041 | -9.90% |
| 13 | 8 | 5,666 | 4,958 | -12.50% |
| 14 | 9 | 6,375 | 5,250 | -17.65% |
| 80 | 9 | 8,792 | 7,417 | -15.64% |
| 15 | 10 | 11,125 | 9,000 | -19.10% |
| 81 | 10 | 15,709 | 13,250 | -15.65% |
| 16 | 11 | 22,750 | 17,709 | -22.16% |
| 17 | 11 | 12,042 | 9,625 | -20.07% |
| 18 | 12 | 25,667 | 19,584 | -23.70% |
| 19 | 12 | 83,541 | 66,917 | -19.90% |
| 20 | 12 | 12,917 | 11,500 | -10.97% |
| 22 | 13 | 58,625 | 42,667 | -27.22% |
| 23 | 14 | 114,917 | 80,625 | -29.84% |
| 24 | 15 | 368,125 | 282,667 | -23.21% |
| 25 | 16 | 658,292 | 472,459 | -28.23% |
| 26 | 17 | 1,854,915 | 1,542,730 | -16.83% |
| 27 | 18 | 3,881,664 | 2,726,772 | -29.75% |
| 35 | 18 | 3,249,415 | 2,871,209 | -11.64% |
| 28 | 19 | 8,944,660 | 5,841,209 | -34.70% |
| 29 | 19 | 11,878,366 | 10,341,522 | -12.94% |
| 30 | 20 | 24,843,626 | 16,009,086 | -35.56% |
| 31 | 21 | 71,566,402 | 44,951,548 | -37.19% |
| 32 | 22 | 192,587,513 | 111,994,100 | -41.85% |
| 33 | 23 | 166,415,201 | 131,979,686 | -20.69% |
| 34 | 24 | 1,569,960,252 | 861,487,587 | -45.13% |
| 70 | 25 | 91,313,698 | 67,941,404 | -25.60% |
| 2203 | 50 | 130,375 | 127,667 | -2.08% |
