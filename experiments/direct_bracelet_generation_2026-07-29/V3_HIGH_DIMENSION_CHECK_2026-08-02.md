# V3 Higher-Dimension Correctness and Speed Check

Date: 2026-08-02  
Source: `compare_bracelet_generators.cpp` SHA-256 `177a5e232864eeff26f08f628ec5cd658b50d1cbeb8f63fd56646772b9688ed6`

## Scope

This is a read-only review and benchmark of experimental `V3BraceletGenerator`. No C++ source was changed during this check.

The build uses the experiment's existing C++17 settings: GCC 16.1.1, `-O3`, `-DNDEBUG`, LTO, loop unrolling, `-march=native`, and `-mcpu=native`. Every run is pinned to performance CPU 2 on the Apple M1 Max Fedora Asahi machine.

## Implementation audit

The implementation matches the binary specialization of Karim, Alamgir, and Husnine's `BraceFD` pseudocode:

- recursion places one one-bit at a time and skips each intervening zero run;
- `period_density` and the position recurrence generate only fixed-density necklaces;
- `emit_final()` places the final one-bit directly at position `n` and implements the binary cases of the paper's `PrintD` completion test;
- `check_reverse()` performs the paper's prefix-versus-reversal pruning;
- `update_reverse_result()` carries the final suffix comparison in constant work per recursion step and delays a zero-run comparison when the mirrored run is not known yet;
- the output path performs no rotation scan or final canonicalization;
- recursion performs no heap allocation; only the experiment's output vector grows, and its capacity is reused between support sizes;
- fixed arrays cover all supported dimensions through 63, and `uint8_t` safely stores every possible density;
- support sizes zero, one, and full support have explicit safe paths.

No correctness defect or clear low-risk speed fix was found. The profiler attributes 98.68% of V3 cycles to the recursive generator itself. There is no separate canonicalization, allocation, or setup hotspot to remove.

## Correctness

The optimized three-way comparison checks V3 against both V1 (fixed-content FKM plus reflection) and the independent direct fixed-content `BraceletFC` adaptation. It compares every support cardinality, rejects duplicates, and verifies that V3's raw representatives are canonical and ascending.

| Check | Highest dimension | Result | Wall time | Maximum RSS |
|---|---:|---|---:|---:|
| Optimized three-way equality | 30 | Passed | 14.61 s | 168.7 MiB |
| ASan + UBSan three-way equality | 28 | Passed | 27.14 s | 341.0 MiB |

At dimension 30 all three generators produced exactly 17,920,859 nonempty bracelet representatives.
The complete output is saved in `results/verification_v3_release_n30.txt` and `results/verification_v3_sanitized_n28.txt`.

## Dimension-30 timing

Each method was measured twice from the same binary in both orders. The table reports the mean of the two native medians.

| Generator | Run 1 | Run 2 | Mean | Nanoseconds per bracelet |
|---|---:|---:|---:|---:|
| V1: FKM plus reflection | 2,955.702 ms | 2,954.379 ms | 2,955.040 ms | 164.89 ns |
| Direct fixed-content | 1,339.138 ms | 1,338.578 ms | 1,338.858 ms | 74.71 ns |
| V3: direct fixed-density | 500.435 ms | 500.261 ms | 500.348 ms | 27.92 ns |

V3 is **5.906 times as fast as V1**, an 83.07% time reduction. It is **2.676 times as fast as the earlier direct fixed-content algorithm**, a 62.63% time reduction.

The repetitions agree within 0.04%, so the result is not timing noise.

## Hardware counters at dimension 30

Perf measured the same four complete generation passes for V1 and V3. Cache-miss events are not supported by this Apple PMU kernel interface.

| Counter | V1 | V3 | V1 / V3 |
|---|---:|---:|---:|
| Cycles | 35,443,257,971 | 5,968,871,055 | 5.94x |
| Instructions | 125,071,173,744 | 34,017,817,816 | 3.68x |
| Branches | 14,691,953,158 | 4,190,403,138 | 3.51x |
| Branch misses | 519,032,200 | 39,042,343 | 13.29x |

The speedup therefore comes from substantially less work and far fewer branch mispredictions, not merely a favorable clock-frequency sample.

## Dimension-35 feasibility

One V3 generation at dimension 35 produced 490,984,487 nonempty bracelets in 12.411 seconds, or 25.28 ns per bracelet. The count independently matches Burnside's formula: 490,984,488 binary bracelets including the all-zero bracelet.

The existing benchmark driver necessarily executes a one-repetition case four times for pilot, calibration, measurement, and checksum. The complete V3 command therefore took 49.64 seconds and peaked at 515.8 MiB. A full V1/direct/V3 dimension-35 comparison would repeat hundreds of millions of stored representatives several times and add minutes without improving the already stable dimension-30 conclusion, so it was not run.

The large memory figure belongs to the equal-shape experiment, which materializes one complete cardinality layer in a vector. Production support generation uses a callback and does not retain a layer.

## Conclusion

V3 is correct over a materially larger exhaustive range and its high-dimension speed gain is substantial. The next meaningful test is production integration with dynamic forbidden-support pruning; further isolated micro-optimization is not justified by this profile.
