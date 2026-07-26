# `find_pos_first_set_bit` Production Call Chain

Last verified: 2026-07-27

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
| `cpp/src/checkstab.cpp:24` | `fracessa::check_stability` | Select the Bee pivot index `m`. |
| `cpp/src/checkstab.cpp:103` | `fracessa::check_stability` | Select the next reduction coordinate. |
| `cpp/include/linalg/copositive_fraction.hpp:103` | `CopositivityCheckerV3::is_copositive_hadeler` | Handle a one-dimensional subset. |
| `cpp/include/linalg/copositive_fraction.hpp:111` | same | Start principal-submatrix row iteration. |
| `cpp/include/linalg/copositive_fraction.hpp:113` | same | Start principal-submatrix column iteration. |

## Analyzer Call Graph

```text
main / fracessa_core.compute_matrix
  -> fracessa::fracessa
     -> fracessa::search_one_support
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

- `cpp/src/findeq.cpp:19`
- `cpp/include/fracessa/matrix_server.hpp:42`
- `cpp/include/fracessa/matrix_server.hpp:71`
- `cpp/include/fracessa/matrix_server.hpp:97`
- `cpp/include/fracessa/matrix_server.hpp:131`

Regenerate the occurrence list with:

```bash
rg -n 'find_pos_first_set_bit|find_pos_next_set_bit|lowest_set_bit_as_bit|extract_set_indices' \
  cpp/include cpp/src
```
