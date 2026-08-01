# `find_pos_first_set_bit` Production Call Chain

Last verified: 2026-08-01

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
| `cpp/src/checkstab.cpp:31` | `fracessa::check_stability` | Select the Bee pivot index `m`. |
| `cpp/src/checkstab.cpp:185` | `fracessa::check_stability` | Select the next reduction coordinate. |
| `cpp/include/linalg/copositive_fraction.hpp:94` | `CopositivityCheckerV3::is_copositive_hadeler` | Handle a one-dimensional subset. |

`find_pos_next_set_bit()` now has no production caller. Copositivity extracts
each larger subset once and indexes its fixed array directly.

## Analyzer Call Graph

```text
main / fracessa_core.compute_matrix
  -> fracessa::fracessa
     -> fracessa::analyze_support
        -> fracessa::check_stability
           -> find_pos_first_set_bit
           -> lowest_set_bit_as_bit
           -> linalg::is_strictly_copositive
              -> CopositivityCheckerV3::is_strictly_copositive
                 -> CopositivityCheckerV3::is_copositive_hadeler
                    -> find_pos_first_set_bit
```

## Separate Index-Extraction Path

`bs64::extract_set_indices()` does not call this helper. It masks to the active
dimension, then repeatedly executes `ctz64(bits)` and `bits &= bits - 1`.
Production callers are:

- `cpp/include/linalg/copositive_fraction.hpp:100`
- `cpp/src/find_candidate_fast.cpp:92`
- `cpp/src/find_candidate_fast.cpp:94`
- `cpp/src/find_candidate_safe.cpp:172`
- `cpp/src/find_candidate_safe.cpp:174`
- `cpp/src/checkstab.cpp:104`

Regenerate the occurrence list with:

```bash
rg -n 'find_pos_first_set_bit|find_pos_next_set_bit|lowest_set_bit_as_bit|extract_set_indices' \
  cpp/include cpp/src
```
