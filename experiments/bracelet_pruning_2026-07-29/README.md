# Circular Bracelet-Pruning Experiment

Date: 2026-07-29

## Result

The experiment is successful as a circular-support search algorithm:

- All tested mathematical results match production.
- 24 of 27 circular verification matrices are faster.
- The median per-matrix change is `-71.99%`.
- The geometric-mean speedup is `4.28x`.
- Verification matrix 34 (`n=24`) improves from `22.132 s` to `0.650 s`.

It is not yet a byte-compatible production replacement. Grouping candidates by
bracelet changes candidate order, `candidate_id`, and `shift_reference` for 23
of the 27 circular fixtures, although the complete mathematical candidate sets
and ESS counts are identical.

## Algorithm

For a circular-symmetric game, rotating or mirroring a support produces an
equivalent case. All supports related by these operations form a binary
bracelet.

The experimental path:

1. Generates fixed-cardinality binary necklaces directly with the FKM
   recursion, rather than scanning all `C(n,k)` masks.
2. Compares each necklace with the smallest rotation of its mirror and keeps
   one representative per bracelet.
3. Solves that representative once.
4. After an exact equilibrium is found, records all distinct rotations and
   reflections of its support as pruning rules. A later support is skipped if
   it contains one of those rules.
5. Reconstructs every equivalent candidate and ESS count from the solved
   representative.

`remove_supersets()` still physically erases supersets in the unchanged
non-circular path. In the circular path, future supports do not exist yet, so it
stores subset rules instead. Its returned `SupportOrbit` is information already
found while registering the rules: the number of distinct rotations and
whether mirroring creates a second rotation family. The caller reuses it for
candidate reconstruction instead of traversing the orbit again.

## Correctness Checks

- All 27 circular verification fixtures: identical ESS counts and identical
  candidate vectors, supports, extended supports, stability results, and
  payoffs after ignoring only enumeration metadata and row order.
- All 24 non-circular verification fixtures: byte-identical complete output.
- 104 additional deterministic circular games at dimensions 2 through 14:
  identical mathematical candidate sets and ESS counts.
- ASan, UBSan, and leak checks over all 27 circular fixtures: no errors.
- The separate generator experiment exhaustively covers every nonempty support
  at dimensions 8 and 10 with no missing or duplicate bracelet orbits.

The normal unsafe numerical filter was used for the paired benchmark and main
differential run. Exact high-dimensional runs were not started; in particular,
matrix 33 remains excluded from exact reruns as instructed.

## Benchmark Method

- Hardware: Apple M1 Max, Fedora Linux Asahi, performance CPU 2.
- Both variants: C++17, `-O3`, `-DNDEBUG`, LTO, loop unrolling, and native CPU
  tuning.
- The production and bracelet binaries use the same parser, candidate solver,
  exact solver, stability code, FLINT, and spdlog build. Only the experimental
  `Supports` header and `fracessa.cpp` differ.
- Existing in-process benchmark harness: matrix parsing occurs once, analyzer
  construction is repeated, and candidates are materialized.
- Each measurement targets about three seconds and reports the median of up to
  seven batches. Variant order alternates by matrix. A call already exceeding
  three seconds is measured once.
- Negative change means the bracelet implementation is faster.

| ID | n | Production (us) | Bracelet (us) | Change |
|---:|---:|---:|---:|---:|
| 8 | 5 | 5.274 | 5.536 | +4.97% |
| 10 | 6 | 14.485 | 7.993 | -44.82% |
| 11 | 7 | 11.849 | 12.218 | +3.11% |
| 12 | 7 | 12.100 | 12.472 | +3.07% |
| 13 | 8 | 77.671 | 21.756 | -71.99% |
| 14 | 9 | 79.159 | 22.627 | -71.42% |
| 15 | 10 | 247.296 | 40.765 | -83.52% |
| 16 | 11 | 75.825 | 66.923 | -11.74% |
| 17 | 11 | 45.153 | 42.027 | -6.92% |
| 18 | 12 | 735.168 | 94.993 | -87.08% |
| 19 | 12 | 2369.221 | 325.777 | -86.25% |
| 20 | 12 | 160.527 | 45.851 | -71.44% |
| 21 | 12 | 1875.129 | 213.036 | -88.64% |
| 22 | 13 | 209.876 | 177.631 | -15.36% |
| 23 | 14 | 1753.079 | 292.898 | -83.29% |
| 24 | 15 | 5117.565 | 828.652 | -83.81% |
| 25 | 16 | 13461.029 | 1171.753 | -91.30% |
| 26 | 17 | 7557.945 | 4942.622 | -34.60% |
| 27 | 18 | 67418.219 | 4046.047 | -94.00% |
| 28 | 19 | 11809.369 | 8748.135 | -25.92% |
| 29 | 19 | 90872.127 | 47926.886 | -47.26% |
| 30 | 20 | 347756.162 | 17867.468 | -94.86% |
| 31 | 21 | 735827.023 | 44748.827 | -93.92% |
| 32 | 22 | 3597337.971 | 105518.776 | -97.07% |
| 33 | 23 | 563055.150 | 322239.502 | -42.77% |
| 34 | 24 | 22132065.526 | 650213.754 | -97.06% |
| 35 | 18 | 140235.215 | 7118.475 | -94.92% |

## Dimension-24 Cross-Check

The large matrix-34 result is consistent with the amount of support work:

- Raw nonempty supports: `16,777,215`.
- Supports stored by the current coprime-only circular reduction: `11,248,061`.
- Bracelet representatives before candidate pruning: `352,697`.
- Representative reduction: `31.89x`.
- Measured in-process speedup: `34.04x`.

An independent one-shot CLI run measured `21.38 s` versus `0.63 s` of CPU
time and returned the same 15,120 ESS and 15,120 mathematical candidates. Peak
resident memory fell from `99,520 KiB` to `13,008 KiB`.

## Decision Still Required

Before production integration, decide whether candidate IDs, row order, and
the historical meaning of `shift_reference` must remain byte-compatible. If
they do, result ordering needs a separate design and benchmark. If a new
deterministic order is acceptable for a new program version, the mathematical
result already passes the experiment's checks.
