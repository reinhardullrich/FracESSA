# Direct Bracelet Generation Experiments

Date: 2026-07-29

## Question

Is the direct fixed-content bracelet generator from Karim, Sawada, Alamgir,
and Husnine faster than the current experimental method?

The two compared methods are:

1. **FKM plus reflection:** generate every fixed-weight binary necklace, find
   the smallest rotation of its reflection, and retain one necklace per
   bracelet. This is the algorithm used by the existing bracelet-pruning
   experiment.
2. **Direct BraceletFC:** maintain reversal order during fixed-content
   necklace recursion and reject impossible branches before completing the
   mirrored necklace. This is an independent C++17 adaptation of the optimized
   algorithm in [Generating Bracelets with Fixed Content](https://www.socs.uoguelph.ca/~sawada/papers/fix-brace.pdf).

Both implementations remain in
`compare_bracelet_generators.cpp`; neither production code nor the earlier
bracelet experiments were changed.

## Correctness

For every dimension `2..24` and every support size `1..n`, the two generators'
outputs were reduced to their smallest rotation/reflection representative only
for untimed comparison. Every orbit matched and the direct output contained no
duplicates. At dimension 24 both methods generated exactly `352,697` bracelet
orbits across all nonempty support sizes.

The complete comparison passed in both the optimized build and an
AddressSanitizer/UndefinedBehaviorSanitizer build with leak detection. Saved
outputs are in `results/verification_release.txt` and
`results/verification_sanitized.txt`.

The direct generator labels the less frequent category as the lower symbol,
as required by the paper's constant-amortized-time condition. Therefore its
raw representative can have a different orientation or complemented symbol
ordering from the FKM representative. This does not change the represented
rotation/reflection orbit.

## Benchmark Method

- Hardware: Apple M1 Max, Fedora Linux Asahi, performance CPU 2.
- Build: C++17, `-O3`, `-DNDEBUG`, LTO, loop unrolling, `-march=native`, and
  `-mcpu=native`.
- Each measured call constructs one generator and creates bracelets for every
  support size `1..n` into the same `std::vector<uint64_t>` output shape.
- Canonicalization used by the correctness comparison is not timed.
- No support pruning, matrix construction, solving, or candidate reconstruction
  is timed.
- Each algorithm/dimension measurement targets three seconds and uses the
  median of seven batches. Algorithm order alternates by dimension.
- Negative change means the direct paper algorithm is faster.

| n | Bracelets | FKM + reflection (us) | Direct (us) | Change |
|---:|---:|---:|---:|---:|
| 2 | 2 | 0.070 | 0.126 | +80.59% |
| 3 | 3 | 0.097 | 0.146 | +50.79% |
| 4 | 5 | 0.161 | 0.104 | -35.20% |
| 5 | 7 | 0.233 | 0.156 | -33.21% |
| 6 | 12 | 0.395 | 0.329 | -16.63% |
| 7 | 17 | 0.630 | 0.604 | -4.14% |
| 8 | 29 | 1.127 | 1.204 | +6.89% |
| 9 | 45 | 1.973 | 2.273 | +15.24% |
| 10 | 77 | 3.519 | 4.175 | +18.64% |
| 11 | 125 | 6.321 | 7.749 | +22.61% |
| 12 | 223 | 11.694 | 13.974 | +19.50% |
| 13 | 379 | 21.453 | 25.851 | +20.50% |
| 14 | 686 | 41.971 | 49.832 | +18.73% |
| 15 | 1,223 | 138.197 | 94.894 | -31.33% |
| 16 | 2,249 | 175.073 | 185.969 | +6.22% |
| 17 | 4,111 | 358.402 | 377.499 | +5.33% |
| 18 | 7,684 | 721.960 | 741.453 | +2.70% |
| 19 | 14,309 | 1,419.753 | 1,384.218 | -2.50% |
| 20 | 27,011 | 2,682.011 | 2,538.333 | -5.36% |
| 21 | 50,963 | 5,131.472 | 4,642.459 | -9.53% |
| 22 | 96,908 | 9,894.982 | 8,556.444 | -13.53% |
| 23 | 184,409 | 19,252.769 | 15,739.138 | -18.25% |
| 24 | 352,697 | 37,712.679 | 29,300.653 | -22.31% |

Raw batch results are stored under `results/raw/`; the paired table is
`results/paired_3s.tsv`.

## Result

The direct algorithm is not a universally faster replacement:

- At dimensions `8..14`, it is `17.34%` slower by geometric mean.
- At dimensions `15..18`, neither method consistently wins; dimension 15 is an
  isolated direct-generator win, while dimensions 16-18 are small regressions.
- At dimensions `19..24`, direct generation is `12.19%` faster by geometric
  mean.
- At dimensions `21..24`, the advantage grows monotonically from `9.53%` to
  `22.31%`.

For the current project, where small matrices matter heavily, the paper
algorithm should not replace FKM unconditionally. It is a credible large-
dimension generator, but its end-to-end value must still be measured inside
the bracelet-pruning search before any production choice.

## V3: direct fixed-density generation

The follow-up implementation called **V3** is an independent binary C++17 adaptation of `BraceFD` from Karim, Alamgir, and Husnine's 2014 paper *Generating fixed density bracelets of arbitrary base*. Unlike the fixed-content generator above, V3 places one one-bit per recursive call and skips the intervening zero run. The final one-bit is placed directly at position `n`, as required by the paper's fixed-density necklace completion test.

V3 was compared with both existing generators for every dimension `1..24` and every support size. It produced exactly the same `352,697` orbits at dimension 24, with no duplicate, noncanonical, or out-of-order representatives. The optimized build and the AddressSanitizer/UndefinedBehaviorSanitizer build both passed through dimension 24. Saved verification output is in `results/verification_v3_release.txt` and `results/verification_v3_sanitized.txt`.

Only V3 was newly timed; the unchanged V1 and direct fixed-content columns below come from the saved three-second run. Negative change means V3 is faster.

| n | Bracelets | V1: FKM + reflection (us) | V3 (us) | V3 vs V1 |
|---:|---:|---:|---:|---:|
| 5 | 7 | 0.233 | 0.151 | -35.31% |
| 6 | 12 | 0.395 | 0.569 | +44.00% |
| 7 | 17 | 0.630 | 0.590 | -6.32% |
| 8 | 29 | 1.127 | 0.700 | -37.88% |
| 9 | 45 | 1.973 | 1.156 | -41.41% |
| 10 | 77 | 3.519 | 1.764 | -49.88% |
| 11 | 125 | 6.321 | 2.686 | -57.50% |
| 12 | 223 | 11.694 | 4.609 | -60.59% |
| 13 | 379 | 21.453 | 8.278 | -61.41% |
| 14 | 686 | 41.971 | 14.712 | -64.95% |
| 15 | 1,223 | 138.197 | 27.043 | -80.43% |
| 16 | 2,249 | 175.073 | 51.306 | -70.69% |
| 17 | 4,111 | 358.402 | 165.041 | -53.95% |
| 18 | 7,684 | 721.960 | 322.803 | -55.29% |
| 19 | 14,309 | 1,419.753 | 682.770 | -51.91% |
| 20 | 27,011 | 2,682.011 | 1,293.187 | -51.78% |
| 21 | 50,963 | 5,131.472 | 2,293.634 | -55.30% |
| 22 | 96,908 | 9,894.982 | 4,089.496 | -58.67% |
| 23 | 184,409 | 19,252.769 | 7,321.643 | -61.97% |
| 24 | 352,697 | 37,712.679 | 13,183.846 | -65.04% |

Across dimensions `5..24`, V3 is `55.30%` faster at the median and `53.31%` faster by geometric mean. At dimension 24 it is `2.86x` as fast as V1. Full dimensions `2..24`, both comparator columns, and raw V3 batches are saved in `results/paired_v3_3s.tsv` and `results/raw/n*_v3.txt`.

V3 subsequently passed its production gate while preserving dynamic forbidden-support pruning and candidate multipliers. V1 and
V2 remain available for comparison, while FracESSA now selects V3 for circular matrices. The end-to-end correctness and timing
record is in `PRODUCTION_V3_COMPARISON_2026-08-03.md`.

## Reproduce

```bash
cd experiments/direct_bracelet_generation_2026-07-29
./build.sh
./builds/compare_bracelet_generators_san verify 24
./benchmark.sh 3
taskset -c 2 ./builds/compare_bracelet_generators benchmark v3 24 3
```

## Primary Sources

- S. Karim, J. Sawada, Z. Alamgir, and S. M. Husnine,
  [Generating Bracelets with Fixed Content](https://www.socs.uoguelph.ca/~sawada/papers/fix-brace.pdf),
  2011. The optimized `BraceletFC` algorithm and complete C reference are in
  Section 2.2 and the appendix.
- J. Sawada,
  [Generating Bracelets in Constant Amortized Time](https://www.cis.uoguelph.ca/~sawada/papers/brace.pdf),
  SIAM Journal on Computing 31(1), 2001.
- S. Karim, Z. Alamgir, and S. M. Husnine,
  [Generating Fixed-Density Bracelets of Arbitrary Base](https://doi.org/10.1080/00207160.2013.805753),
  *International Journal of Computer Mathematics* 91(3), 2014.
