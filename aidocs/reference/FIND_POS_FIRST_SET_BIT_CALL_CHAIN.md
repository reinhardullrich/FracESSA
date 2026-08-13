# `find_pos_first_set_bit` Production Call Chain

Status: historical production snapshot from 2026-08-09. Tests were excluded. The 2026-08-12 coposit integration removed the
FracESSA-owned Hadeler and connected-component callers described below; current coposit internals live in `external/coposit/`.

## Contract

- `bs64::find_pos_first_set_bit(bits)` returns the zero-based position of the lowest set bit. It delegates to branch-free `ctz64()`
  and requires `bits != 0`.
- `bitset_multiword::find_pos_first_set_bit()` scans words from low to high, applies `ctz64()` to the first nonzero word, and requires
  a nonempty set.
- `bs64::lowest_set_bit_as_bit()` is separate: it isolates the lowest bit directly and is defined for zero.
- The former `find_pos_next_set_bit()` helper no longer exists. Iteration clears the current lowest bit or extracts all indices once.

## Direct Production Calls At The Snapshot

| Location | Caller | Purpose |
|---|---|---|
| `cpp/src/fracessa.cpp` | `basic_fracessa::log_candidate` | Find the deterministic reference strategy for optional diagnostic logging. |
| `cpp/src/fast_candidate_filter.cpp` | `fast_candidate_filter::passes_from_indices` | Visit outside-support strategies in the one-word binary64 payoff check. |
| `cpp/include/linalg/copositive_integer.hpp` | `CopositivityChecker::is_strictly_copositive_hadeler` | Read the only index of a one-coordinate principal subset. |
| `cpp/include/linalg/copositive_integer.hpp` | `CopositivityChecker::are_negative_components_strictly_copositive` | Visit vertices while constructing negative-entry graph components. |
| `cpp/include/fracessa/support_generator_non_circular.hpp` | `NonCircularSupportGenerator<bitset_multiword>::activate_pending` | Bucket a newly forbidden multiword support by its lowest strategy. |
| `cpp/include/fracessa/support_generator_circular.hpp` | `CircularSupportGenerator<bitset_multiword>::activate_pending` | Bucket a newly forbidden multiword orbit member by its lowest strategy. |

Every caller establishes nonemptiness through its loop, candidate, principal-subset, component, or forbidden-support invariant.

## Historical Call Graph

```text
basic_fracessa::run
  -> basic_fracessa::analyze_support
     -> basic_fracessa::log_candidate (logging only)
        -> one-word or multiword find_pos_first_set_bit
     -> basic_fracessa::check_stability
        -> CopositivityChecker::are_negative_components_strictly_copositive
           -> one-word or multiword find_pos_first_set_bit
           -> CopositivityChecker::is_strictly_copositive
              -> CopositivityChecker::is_strictly_copositive_hadeler
                 -> bs64::find_pos_first_set_bit (one-coordinate one-word subset only)

basic_fracessa::analyze_support
  -> fast_candidate_filter::passes
     -> bs64::find_pos_first_set_bit (one-word outside-strategy iteration)

multiword support generator::generate
  -> activate_pending
     -> bitset_multiword::find_pos_first_set_bit
```

## Separate Index-Extraction Path At The Snapshot

Neither `extract_set_indices()` implementation calls `find_pos_first_set_bit()`. Both repeatedly apply `ctz64()` and clear the bit
with `word &= word - 1`; the multiword version repeats that operation for every word.

Production used index extraction in:

- `cpp/include/linalg/copositive_integer.hpp` for Hadeler subsets and disconnected components;
- `cpp/src/exact_candidate_solver.cpp` for candidate systems, scaled reduced $B$, and outside-payoff checks;
- `cpp/src/fast_candidate_filter.cpp` for binary64 candidate systems; and
- `cpp/include/fracessa/circular_affine_symmetry.hpp` for multiword support permutations.

Regenerate the current occurrence list with:

```bash
rg -n 'find_pos_first_set_bit|find_pos_next_set_bit|lowest_set_bit_as_bit|extract_set_indices' \
  cpp/include cpp/src
```
