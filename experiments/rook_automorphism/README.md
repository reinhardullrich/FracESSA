# Rook-Automorphism Integer Search Plan

## Goal

Search the three-variable compact rook matrices from
[`research/ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md`](../../research/ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md) over expanding integer
triples $(A,B,C)$ and retain the largest safely verified ESS count found for every dimension and rook group.

The previous numerical seed was

$$
(A,B,C)=(7,7,15).
$$

For the published $3\times8$, $n=24$ game this gives 15,120 ESS. That is the only established rook-family incumbent in the source
notes; the other $7/7/15$ templates are proposed inputs whose ESS counts still need to be measured.

## Initial Scope

Use only rows with a literal three-variable compact string. This gives 15 grids:

$$
2\times3,
2\times5,
3\times4,
2\times7,
3\times5,
2\times9,
4\times5,
3\times7,
2\times11,
3\times8,
2\times13,
4\times7,
2\times15,
3\times10,
5\times6.
$$

Their dimensions are $6,10,12,14,15,18,20,21,22,24,26,28,30,30,30$. Prime dimensions and non-coprime grids have no compact
three-variable string and are outside the first experiment. Full symmetric rook matrices can be a separate experiment later.

Do not require $A$, $B$, and $C$ to be distinct. Equalities enlarge the automorphism group but may still produce the best ESS
count and are part of the requested integer search.

## Deterministic Value Stages

Use the scalar order

```text
0, +1, -1, +2, -2, +3, -3, ...
```

Each newly introduced scalar is one **value stage**. At stage $v$, form all triples from the values introduced through $v$, but
emit only triples containing $v$ at least once. Earlier triples were completed in earlier stages and are not generated again.

For example:

| Stage | Available values | New raw triples before exact redundancy removal |
|---|---|---:|
| `0` | `0` | 1 |
| `+1` | `0,+1` | $2^3-1^3=7$ |
| `-1` | `0,+1,-1` | $3^3-2^3=19$ |
| `+2` | `0,+1,-1,+2` | $4^3-3^3=37$ |
| `-2` | `0,+1,-1,+2,-2` | $5^3-4^3=61$ |

For each value stage, generate one deterministic stream ordered by dimension, grid rank, and triple order. Send that complete
stream directly to one multiprocessing run. Workers take whichever matrix is next; there is no batching, waiting, or separate run
by family or dimension. The `+2` stage is complete only after every enabled matrix through $n=30$ has returned. Then the `-2`
stage begins.

The value stage controls generation order only; it is not a save boundary.

The search is infinite, so the runner takes an explicit final value stage. Increasing it continues the same result store rather
than starting again.

## Exact Redundancies

The requested combinations can be covered without solving or storing mathematically identical games repeatedly:

1. Emit a nonzero triple only when $\gcd(|A|,|B|,|C|)=1$. Positive common scaling preserves all candidates and ESS, so later
   multiples are skipped rather than divided and stored again. Do not identify a triple with its negative: negative scaling
   exchanges payoff maxima and minima.
2. The $m\times k$ game $(A,B,C)$ is permutation-similar to the transposed $k\times m$ game $(B,A,C)$. They have the same ESS
   count. The deterministic generator therefore emits a triple only when $A$ occurs no later than $B$ in
   `0, +1, -1, +2, -2, ...`. When it later reaches $(B,A,C)$, that job is skipped. No swapped alias is stored.
3. Handle $(0,0,0)$ analytically as zero ESS instead of sending the maximally degenerate zero game through the support search.

This covers every distinct ESS outcome requested by the integer search. Only the primitive, first-generated representative is
stored, using the table's `A = same row`, `B = same column`, `C = neither` convention. There is no alias table.

## Multiprocessing Pipeline

Use the maintained PyFracESSA multiprocessing API:

```text
generator of Matrix jobs
    -> run_multiprocessing("fast", jobs, include_candidates=False)
    -> one parent result consumer
    -> in-memory progress and record update
    -> atomic state.json replacement every checkpoint interval
```

Each job carries `n`, `m`, `k`, template rank, value-stage index/value, and $(A,B,C)$ in `Matrix.metadata`. Workers only compute. The
parent is the sole writer.

Use one experiment file, `state.json`; do not create an append log or experiment database. The Python main thread keeps current
state in memory and atomically replaces this file through a temporary file every 10 seconds during initial testing. Make the
interval configurable as `--checkpoint-seconds`; use 30 seconds for long searches once the runner is established. Also write on
normal completion and on a handled interrupt.

The root of `state.json` is grouped first by decimal dimension and then by the `m x k` rook group. Each group has exactly two
objects, `work` and `result`:

```json
{
  "6": {
    "2x3": {
      "work": {
        "active_stage": null,
        "last_received": null,
        "completed_through": null,
        "received_in_stage": 0,
        "expected_in_stage": 0,
        "best_unverified_fast": null
      },
      "result": null
    }
  }
}
```

For every dimension plus rook group, `work` stores:

- `active_stage`: the index and value of the stage currently being searched;
- `last_received`: the most recently returned job index, value stage, and $(A,B,C)$, for human-visible progress;
- `completed_through`: the contiguous completion frontier—the largest deterministic job index for which that job and every
  smaller job index in the group have returned successfully;
- the value stage and $(A,B,C)$ belonging to `completed_through`;
- received and expected counts for the active value stage.

Multiprocessing completes jobs out of order, so `last_received` alone cannot be a resume point: an earlier job may still be
running. In memory, each group keeps a set of returned job indices above `completed_through`. Whenever the next consecutive index
arrives, advance the frontier and continue consuming consecutive indices already present in that set.

The later in-memory set is deliberately not written to disk. On restart, regenerate the deterministic stream from
`completed_through + 1`; jobs that had finished beyond a missing older job are recomputed. This may repeat more than one checkpoint
interval when one old job runs much longer than the following jobs, but it can never skip an unfinished matrix. That is the chosen
trade-off for keeping `state.json` small and simple.

`result` is `null` until the group has a safely verified maximum. It then stores result details only for that specific dimension
and rook group. Every maximum contains:

- `matrix`: the complete substituted FracESSA matrix string, including its `dimension#` prefix;
- `A`, `B`, and `C`: the integer triple that produced it;
- `candidate_count` and `ess_count`: the complete totals;
- `candidate_structure`, mapping each support size to its complete candidate count;
- `ess_structure`, mapping each support size to its complete ESS count;
- `safe_elapsed_ns`: the native safe timing;
- `safe_fallback`: the safe-fallback value.

The structures use JSON objects with decimal support sizes as keys, matching the canonical matrix database convention. The file
never stores candidate vectors, individual supports, payoffs, or candidate rows.

After a value stage is complete, report:

- the best safely verified ESS count and its first record-setting triple for each dimension and rook group;
- $\gamma=\operatorname{ESS}^{1/n}$ for comparisons across dimensions;
- elapsed work, completed triples, failures, and the largest completely finished value stage.

The runner accepts `--max-dimension`, initially 30. If $n=30$ or $n=28$ proves impractical, restart with a lower value. Existing
higher-dimensional results and progress remain stored and can be resumed later by raising the limit again.

The maintained multiprocessing API deliberately has no per-matrix timeout. Keep the first runner on that simple API. If a
high-dimensional job actually blocks useful progress, add the
calibration runner's disposable-child timeout pattern in a second step; do not build a custom replaceable-worker pool before it is
needed. Stopping and restarting an incomplete value stage may repeat returned jobs above a stalled contiguous frontier, but it
never loses a verified maximum or skips an unfinished job.

## Fast Search And Strict-Record Verification

`fast` returns exact results for the supports it keeps, but its floating-point rejection is not a proof and can undercount ESS.
The experiment deliberately uses this fixed workflow:

1. Run every generated triple with `fast`.
2. Give each returned ESS count and job key to the Python main thread, which updates the in-memory completed ranges,
   `last_received`, and current fast maximum for that dimension and rook group.
3. Run that triple with `safe` immediately only when its fast count is strictly larger than the current safely verified record for
   its dimension and rook group. The native result supplies the complete candidate and ESS counts and structures without returning
   individual candidate rows.
4. Replace the record only with the safe result. Keep the new record in memory and include it in the next periodic atomic
   `state.json` replacement; print it to standard output immediately.
5. An equal fast count is a tie, not a new record, and does not trigger `safe`.

Fast results arrive in multiprocessing completion order. That order can change how many intermediate strict records require a
safe check, but it cannot turn an equal count into a record. A checkpoint taken before a safe result returns contains the completed
fast job but not the unverified record. Therefore each group's `work` object also keeps its best unverified fast tuple; startup
safe-checks any such tuple that exceeds that group's verified record before resuming fast jobs.

The $3\times8$ seed $(7,7,15)$ initializes the verified `24 -> 3x8` record at 15,120 ESS. The same seed is run once for every
other grid to establish its fast result; it triggers safe verification because every other dimension-and-group result starts at
$-1$. Thus each group's first completed fast result establishes a real safe-verified baseline even when its ESS count is zero.

Every retained record is exact because it has passed `safe`, but this remains a heuristic search: a false rejection in `fast`
could hide a tuple whose true safe count would have set a record.

## Record Validation

Whenever safe verification establishes a new record:

1. store its exact ESS count and timing immediately;
2. rerun the final winner once in `safe` with candidates enabled before considering it for the canonical matrix database;
3. compare its complete candidate output on the repeated run, not only the ESS count.

## First Implementation

Keep the first implementation to one Python runner and `state.json`. It needs only:

- template substitution from the 15 literal compact strings;
- value-stage generation with primitive-triple and swapped-$A/B$ suppression;
- the existing `fracessa.run_multiprocessing()` API;
- one contiguous completion frontier plus an in-memory later-completion set per rook group, atomic `state.json`, and resumability;
- a compact best-results report.

Do not add a generic automorphism library, add an experiment database, or modify the canonical test database for this search.

## Running The Search

Run the source-only self-check without loading the native extension:

```bash
python3 experiments/rook_automorphism/search.py --self-check
```

Search through the `-2` value stage for every included group through dimension 30:

```bash
python3 experiments/rook_automorphism/search.py --through=-2
```

The default checkpoint interval is 10 seconds and the default state file is `experiments/rook_automorphism/state.json`. Use
`--checkpoint-seconds 30` for a long run and `--max-dimension` to postpone the largest groups. Restarting the same command resumes
from the contiguous completion frontier stored in that file.
