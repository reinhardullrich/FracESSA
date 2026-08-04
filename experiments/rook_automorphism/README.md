# Rook-Automorphism Integer Search Plan

## Goal

Search the three-variable compact rook matrices from
[`research/ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md`](../../research/ROOK_AUTOMORPHISM_GROUPS_N2_TO_N30.md) over expanding integer
triples $(A,B,C)$ and retain the largest ESS count found for every dimension and rook group.

The previously known numerical example was

$$
(A,B,C)=(7,7,15).
$$

For the published $3\times8$, $n=24$ game this gives 15,120 ESS. The current run retains already verified $(7,7,15)$ results as
incumbents, but they do not advance the search frontier: enumeration starts at $0,+1,-1,+2,-2,\ldots$.

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

Give every dimension and rook group the same deterministic job indices. Merge those streams by job index and then by ascending
dimension and grid rank, and send the merged stream to one multiprocessing run. A value stage controls queue priority only; it is
neither a barrier nor a save boundary. If one matrix from `-1` is still running, idle workers immediately continue with queued `+2`
jobs from any group.

The search is effectively infinite. The default final value stage is `+1000000`; stages are generated lazily, so this does not
precompute or store a million stages. The bounded multiprocessing window stays full while the saved per-group frontiers remain the
restart authority.

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

The root of `state.json` uses one `dimension:group` key for each search, for example `10:2x5`. Each group has exactly two objects,
`work` and `result`:

```json
{
  "10:2x5": {
    "work": {
      "last_received": null,
      "completed_through": null
    },
    "result": null
  }
}
```

For every dimension plus rook group, `work` stores:

- `last_received`: the most recently returned job index, value stage, and $(A,B,C)$, for human-visible progress;
- `completed_through`: the contiguous completion frontier—the largest deterministic job index for which that job and every
  smaller job index in the group have returned successfully;
- the value stage and $(A,B,C)$ belonging to `completed_through`.

Multiprocessing completes jobs out of order, so `last_received` alone cannot be a resume point: an earlier job may still be
running. In memory, each group keeps a set of returned job indices above `completed_through`. Whenever the next consecutive index
arrives, advance the frontier and continue consuming consecutive indices already present in that set.

The later in-memory set is deliberately not written to disk. On restart, regenerate the deterministic stream from
`completed_through + 1`; jobs that had finished beyond a missing older job are recomputed. This may repeat more than one checkpoint
interval when one old job runs much longer than the following jobs, but it can never skip an unfinished matrix. That is the chosen
trade-off for keeping `state.json` small and simple.

`result` is `null` until the group has a maximum. It then stores result details only for that specific dimension
and rook group. Every maximum contains:

- `matrix`: the complete substituted FracESSA matrix string, including its `dimension#` prefix;
- `variables`: the integer triple $[A,B,C]$ that produced it;
- `candidate_count` and `ess_count`: the complete totals;
- `candidate_structure`, mapping each support size to its complete candidate count;
- `ess_structure`, mapping each support size to its complete ESS count;
- `fast_elapsed_ns`: the native fast timing.

The structures use JSON objects with decimal support sizes as keys, matching the canonical matrix database convention. Variable
triples and both structure objects are kept on one line to make the state easier to scan. The file never stores candidate vectors,
individual supports, payoffs, or candidate rows.

After the requested search range is complete, report:

- the best ESS count and its first record-setting triple for each dimension and rook group;
- $\gamma=\operatorname{ESS}^{1/n}$ for comparisons across dimensions;
- the final requested value stage.

The runner accepts `--max-dimension`, initially 30. If $n=30$ or $n=28$ proves impractical, restart with a lower value. Existing
higher-dimensional results and progress remain stored and can be resumed later by raising the limit again.

The maintained multiprocessing API deliberately has no per-matrix timeout. Keep the first runner on that simple API. If a
high-dimensional job actually blocks useful progress, add the
calibration runner's disposable-child timeout pattern in a second step; do not build a custom replaceable-worker pool before it is
needed. Stopping and restarting an incomplete search may repeat returned jobs above a stalled contiguous frontier, but it
never loses a maximum or skips an unfinished job.

## Fast Search

`fast` returns exact results for the supports it keeps. Its floating-point rejection can miss candidates and therefore hide a
record, but it cannot invent candidates. The experiment deliberately uses this fixed workflow:

1. Run every generated triple with `fast`.
2. Give each result and job key to the Python main thread, which updates the completed ranges and `last_received`.
3. If its ESS count is strictly larger than the current record, retain its complete counts, structures, matrix, variables, and fast
   timing; print it immediately and include it in the next periodic atomic `state.json` replacement.
4. Treat an equal count as a tie, not a new record.

Fast results arrive in multiprocessing completion order. That order can change which equal maximum is retained first, but it cannot
change the maximum count. The zero game is handled analytically. Existing incumbents affect only the record to beat; they never mark
any search triple as completed. This remains a heuristic search because a false rejection in `fast` could hide a true record.

## Record Validation

Before adding a final winner to the canonical matrix database:

1. rerun it once in `safe` with candidates enabled;
2. compare its complete candidate output, not only the ESS count.

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

Start the effectively continuous search through the default `+1000000` limit:

```bash
python3 experiments/rook_automorphism/search.py
```

The runner uses CPUs 3 through 8 by default, leaving CPUs 0, 1, 2, and 9 free. Pass plain CPU IDs to override that set:

```bash
python3 experiments/rook_automorphism/search.py --cpu 3 4 5 6
```

Run only one or more selected dimensions while preserving every other saved frontier:

```bash
python3 experiments/rook_automorphism/search.py --dimension 14 18 --through 1000000000
```

One worker is created per selected CPU, and Linux CPU affinity restricts the parent and inherited workers to that set. Startup,
pipeline progress, new records, and the final summary are printed to the
terminal. The default checkpoint interval is 10 seconds and the default state file is `experiments/rook_automorphism/state.json`.
Use `--checkpoint-seconds 30` for a long run, `--dimension` for an explicit set, `--max-dimension` to postpone the largest groups,
or `--through` for a deliberately bounded test. `--dimension` and `--max-dimension` are mutually exclusive. Restarting the same
command resumes each selected group from its contiguous completion frontier stored in that file; unselected groups remain unchanged.
