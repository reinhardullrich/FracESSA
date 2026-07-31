# Copositivity Subset Generation

Date: 2026-07-31

## Question

For the exact strict-copositivity test, which cardinality-first enumeration is
the cheapest way to visit all `2^n - 1` nonempty principal subsets?

The standalone benchmark compares:

1. **Eager buckets:** create every subset, group by cardinality, then iterate;
2. **Streaming Gosper:** generate complete fixed-cardinality masks one at a time;
3. **Bitwise DFS:** recursively exclude/include bits from high to low, matching
   the non-circular production generator without its forbidden-support pruning.

No matrix construction, rational arithmetic, Hadeler test, or early rejection
is included. This isolates only subset generation plus one identical count/sum
operation per emitted mask.

## Method

- Dimensions: 2 through 30.
- Hardware: Apple M1 Max, Fedora Linux Asahi, performance CPU 2.
- Build: C++17, `-O3`, `-DNDEBUG`, LTO, loop unrolling, and native CPU flags.
- One process handles every dimension.
- At each dimension, every algorithm runs once untimed as a warm-up.
- Each algorithm is then timed once inside the process; order rotates by
  dimension.
- All three runs must match the exact subset count and checksum before the next
  dimension begins.

This first pass intentionally has no repeated batches or median. Small-dimension
times are therefore only orientation data; repeat measurements are needed before
choosing production code.

## Hot Run Through Dimension 30

The complete output is in `results/one_hot_run_cpu2.tsv`. Selected larger
dimensions are shown below; change is relative to eager buckets.

| n | Subsets | Eager (us) | Gosper (us) | Change | DFS (us) | Change |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 33,554,431 | 69,376.233 | 157,732.294 | +127.36% | 129,309.843 | +86.39% |
| 26 | 67,108,863 | 140,674.590 | 305,683.258 | +117.30% | 257,484.687 | +83.04% |
| 27 | 134,217,727 | 275,756.015 | 617,294.931 | +123.86% | 512,549.332 | +85.87% |
| 28 | 268,435,455 | 548,902.447 | 1,229,734.946 | +124.04% | 1,032,149.370 | +88.04% |
| 29 | 536,870,911 | 1,256,224.523 | 2,499,615.717 | +98.98% | 2,069,081.781 | +64.71% |
| 30 | 1,073,741,823 | 2,676,570.089 | 5,041,550.505 | +88.36% | 4,191,859.049 | +56.61% |

For a complete traversal, eager buckets were fastest from dimension 7 onward
in this run. At dimension 30, Gosper took 1.88 times the eager time and DFS took
1.57 times the eager time; DFS was 16.9% faster than Gosper.

This does not make eager storage suitable for copositivity. Its stored masks
alone occupy about 8 GiB at dimension 30, and its memory requirement continues
to double with every dimension. It also must construct the complete frontier
before discovering an immediately failing 1-by-1 principal submatrix. Gosper
and DFS use constant generator memory and can stop at the first rejection.

The result answers only the full-enumeration generator question. It does not yet
measure the generators around the real exact principal-submatrix work or compare
their early-rejection behavior.

## Eager Phase Split At Dimension 30

A separate hot run divided the eager call into bucket creation, traversal, and
cleanup. Cleanup is reported because the earlier complete eager time included
destruction and deallocation before the function returned.

| Phase | Time | Share of complete eager call |
|---|---:|---:|
| Create and fill cardinality buckets | 2.876875 s | 83.85% |
| Iterate through all buckets | 0.394295 s | 11.49% |
| Destroy and free buckets | 0.159746 s | 4.66% |
| **Total** | **3.430916 s** | **100.00%** |

Considering only the two active steps, creation took 87.95% and the traversal
loop took 12.05%. The absolute total differs from the earlier one-shot eager
measurement because this was a separate invocation without repeated batches;
the useful result is the phase split measured inside this single call.

The profile is saved in `results/eager_profile_30_cpu2.tsv`.

## Medians From Dimension 3 Through 25

Each algorithm has 100 hot samples per dimension. Small calls are batched so
each sample visits about 310,000 subsets; each batch time is divided by its
number of complete calls. This gives exactly 10,000 calls per sample at
dimension 5 and one call per sample from dimension 19 onward. Algorithm order
rotates between samples. Every batch matched the exact count and checksum.

Times are nanoseconds per complete traversal; changes are relative to eager.

| n | Subsets | Calls/sample | Eager (ns) | Gosper (ns) | Change | DFS (ns) | Change |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 7 | 44,286 | 72.959 | 7.011 | -90.39% | 20.426 | -72.00% |
| 4 | 15 | 20,667 | 93.593 | 18.230 | -80.52% | 49.433 | -47.18% |
| 5 | 31 | 10,000 | 124.475 | 49.694 | -60.08% | 107.417 | -13.70% |
| 6 | 63 | 4,921 | 173.500 | 131.837 | -24.01% | 222.194 | +28.07% |
| 7 | 127 | 2,441 | 252.373 | 353.509 | +40.07% | 457.847 | +81.42% |
| 8 | 255 | 1,216 | 401.710 | 989.395 | +146.30% | 932.805 | +132.21% |
| 9 | 511 | 607 | 704.215 | 2,211.764 | +214.08% | 1,882.928 | +167.38% |
| 10 | 1,023 | 304 | 1,355.127 | 4,447.709 | +228.21% | 3,796.531 | +180.16% |
| 11 | 2,047 | 152 | 2,633.773 | 9,053.727 | +243.76% | 7,577.023 | +187.69% |
| 12 | 4,095 | 76 | 5,247.526 | 18,143.092 | +245.75% | 15,134.309 | +188.41% |
| 13 | 8,191 | 38 | 10,674.355 | 36,454.487 | +241.51% | 30,316.882 | +184.02% |
| 14 | 16,383 | 19 | 21,848.684 | 72,794.947 | +233.18% | 60,546.053 | +177.12% |
| 15 | 32,767 | 10 | 44,337.500 | 145,410.400 | +227.96% | 121,027.100 | +172.97% |
| 16 | 65,535 | 5 | 90,000.000 | 291,191.600 | +223.55% | 241,666.600 | +168.52% |
| 17 | 131,071 | 3 | 181,229.000 | 582,659.667 | +221.50% | 484,444.333 | +167.31% |
| 18 | 262,143 | 2 | 396,468.750 | 1,166,895.500 | +194.32% | 971,020.250 | +144.92% |
| 19 | 524,287 | 1 | 874,437.000 | 2,334,374.500 | +166.96% | 1,943,708.000 | +122.28% |
| 20 | 1,048,575 | 1 | 1,862,646.000 | 4,742,332.000 | +154.60% | 3,918,978.500 | +110.40% |
| 21 | 2,097,151 | 1 | 3,862,895.000 | 9,891,560.000 | +156.07% | 7,934,248.000 | +105.40% |
| 22 | 4,194,303 | 1 | 8,135,498.000 | 19,478,370.500 | +139.42% | 16,276,516.500 | +100.07% |
| 23 | 8,388,607 | 1 | 16,334,162.500 | 39,001,823.500 | +138.77% | 32,469,721.000 | +98.78% |
| 24 | 16,777,215 | 1 | 33,031,554.500 | 77,817,501.500 | +135.59% | 64,899,047.000 | +96.48% |
| 25 | 33,554,431 | 1 | 68,567,837.000 | 153,740,024.500 | +124.22% | 128,442,822.500 | +87.32% |

Gosper is fastest through dimension 6. Eager buckets become fastest at
dimension 7 for a complete traversal. DFS overtakes Gosper at dimension 8 but
does not overtake eager in this range.

The exact output is in
`results/all_generators_3_to_25_median_100_cpu2.tsv`. The earlier isolated
dimension-5 run remains in `results/all_generators_5_median_100_cpu2.tsv`.

## Reproduce

```bash
./experiments/copositivity_subset_generation_2026-07-31/run_once.sh
```

The saved results are `results/one_hot_run_cpu2.tsv`,
`results/eager_profile_30_cpu2.tsv`,
`results/all_generators_5_median_100_cpu2.tsv`, and
`results/all_generators_3_to_25_median_100_cpu2.tsv`.
