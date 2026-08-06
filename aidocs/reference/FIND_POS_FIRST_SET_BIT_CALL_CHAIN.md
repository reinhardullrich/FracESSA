# `find_pos_first_set_bit` Production Call Chain

Last verified: 2026-08-06

Scope: production C++ under `cpp/include/` and `cpp/src/`; tests are excluded.

## Contract

`bs64::find_pos_first_set_bit(bits)` returns the zero-based position of the
lowest set bit. It delegates to the branch-free `ctz64()` primitive and requires
`bits != 0`. Every production caller below establishes that precondition.
`lowest_set_bit_as_bit()` is separate: it directly isolates the lowest bit and
is defined for zero.

## Direct Calls

| Location | Caller | Purpose |
|---|---|---|
| `cpp/src/checkstab.cpp:22` | `fracessa::check_stability` | Report the deterministic reference index `m` when logging is enabled. |
| `cpp/src/find_candidate_fast.cpp:443` | `find_candidate_fast::find` | Visit the next outside-support strategy in the heuristic payoff check. |
| `cpp/src/find_candidate_test.cpp:443` | `find_candidate_test::find` | Visit the next outside-support strategy in the experimental payoff check. |

`find_pos_next_set_bit()` now has no production caller. The adaptive-cone
copositivity checker does not scan support bitsets.

## Analyzer Call Graph

```text
main / fracessa_core.compute_matrix
  -> fracessa::fracessa
     -> fracessa::analyze_support
        -> fracessa::check_stability
           -> find_pos_first_set_bit
           -> lowest_set_bit_as_bit
           -> CopositivityChecker::is_strictly_copositive

fracessa::analyze_support
  -> find_candidate_fast::find / find_candidate_test::find
     -> find_pos_first_set_bit
```

## Separate Index-Extraction Path

`bs64::extract_set_indices()` does not call this helper. It masks to the active
dimension, then repeatedly executes `ctz64(bits)` and `bits &= bits - 1`.
Production callers are:

- `cpp/src/find_candidate_fast.cpp:378`
- `cpp/src/find_candidate_safe.cpp:226`
- `cpp/src/find_candidate_safe.cpp:227`
- `cpp/src/find_candidate_safe.cpp:299`
- `cpp/src/find_candidate_safe.cpp:301`
- `cpp/src/find_candidate_test.cpp:378`

Regenerate the occurrence list with:

```bash
rg -n 'find_pos_first_set_bit|find_pos_next_set_bit|lowest_set_bit_as_bit|extract_set_indices' \
  cpp/include cpp/src
```
