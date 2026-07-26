# FracESSA computer timing comparison

> Frozen benchmark snapshot from 2026-07-26. It records that experiment, not
> the current worktree.

## Bottom line

The current native build on the Apple M1 Max is the fastest recorded run found.
On five hard matrices present in every result set (IDs 30, 31, 32, 33, and 35),
the observed historical run combinations are:

| Saved run | Relative time | Interpretation |
|---|---:|---|
| Current M1 Max | 1.00x | Current optimized worktree |
| Fastest prior JSON (`20260103_031737`) | 1.75x slower | Fastest of 142 saved historical runs |
| Latest prior JSON (`20260212_003121`) | 1.92x slower | Latest saved historical run |
| Earliest prior JSON (`20251211_193114`) | 3.55x slower | Earliest saved historical run |
| Old REF baseline | 6.05x slower | Old modified REF executable |

These are observed computer-plus-program comparisons, not pure CPU benchmarks.
The old JSON files do not record hostname, CPU, compiler, commit, or build flags,
so the effect of hardware cannot be separated from code changes.

## Hard matrices

Seconds reported by the executables:

| Matrix ID | M1 Max now | Latest prior | Latest / now | Old REF baseline | Baseline / now |
|---:|---:|---:|---:|---:|---:|
| 30 | 0.328 | 0.530 | 1.61x | 2.165 | 6.59x |
| 31 | 0.657 | 1.188 | 1.81x | 4.628 | 7.05x |
| 32 | 3.325 | 7.793 | 2.34x | 15.513 | 4.67x |
| 33 | 0.274 | 0.490 | 1.79x | 1.735 | 6.34x |
| 34 | 21.039 | 53.415 | 2.54x | 96.841 | 4.60x |
| 35 | 0.105 | 0.222 | 2.11x | 0.620 | 5.89x |

Across all 35 matrix timings, the current sum is 25.855 seconds, the latest
prior sum is 63.822 seconds (2.47x slower), and the old baseline sum is
122.205 seconds (4.73x slower).

## Current run

- Machine: Apple M1 Max, ARM64 Linux, performance core 2.
- Compiler: GCC 16.1.1.
- Optimization: `-O3`, LTO, `-funroll-loops`, `-march=native`, and
  `-mcpu=native`; no `-ffast-math` or `-Ofast` correctness relaxation.
- Protocol: safe parser, candidates enabled, internal C++ timing, repeated
  small cases, and all 35 expected ESS counts matched.
- Raw result: `experiments/speed_comparison_2026-07-26/results/current_worktree.json`.

## Historical coverage

- FracESSA: all 142 timestamped JSON result files from December 2025 through
  February 2026 were considered. Their hard-matrix profiles range from 1.75x
  to 4.40x the current timing.
- GitHub: tags and release builds do not contain additional machine timing
  records; expired Actions artifacts cannot establish runner timings.
- C# EFR: projects and `Stopwatch` timing code exist back to 2016, but no saved
  numeric timing output was found in readable text files. They establish the
  project history, not a defensible numeric comparison.

## `binfmt` error

The archived REF executable is Linux x86-64, while this M1 Max system is ARM64.
Linux handed the foreign executable to its `binfmt_misc` dispatcher, but no
x86-64 emulator such as FEX or QEMU is registered. This is an architecture
dispatch error, not a FracESSA error. Historical binaries were not rerun.
