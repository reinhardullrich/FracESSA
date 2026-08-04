# Dimension-24 Rook Symmetry Breaking

This experiment starts from the published $3\times8$ rook game with 15,120 ESS and reduces its column symmetry from arbitrary
permutations to rotations and reflections. Strategies are cells $(r,c)$ with three rows and eight circularly ordered columns.

The nine off-diagonal parameters are:

- $u_1,\ldots,u_4$: same row and circular column distance $1,\ldots,4$;
- $v_0,\ldots,v_4$: different rows and circular column distance $0,\ldots,4$.

Their compact FracESSA order is

$$
(v_1,v_2,u_3,v_4,v_3,u_2,v_1,v_0,u_1,v_2,v_3,u_4).
$$

The scaled baseline is

$$
(u_1,u_2,u_3,u_4,v_0,v_1,v_2,v_3,v_4)=(700,700,700,700,700,1500,1500,1500,1500),
$$

which is strategically equivalent to the published $(7,15)$ game. Scaling makes an integer perturbation of one a small rational
change. The runner first checks the baseline and all 18 individual $\pm1$ changes. It then performs larger coordinate sweeps, pair
sweeps, and reproducible random local perturbations around the best result found so far. It stores only the ESS count and timing;
candidate rows are deliberately not requested or retained.

By explicit decision this is a `fast`-only exploratory search with no separate safe verification and no candidate output. Fast
rejection can miss a true candidate, so the experiment may undercount. It never writes to the canonical matrix database.

Run the source self-check:

```bash
python3 experiments/rook_symmetry_breaking_2026-08-04/run.py --self-check
```

Run continuously on CPU 2 until interrupted:

```bash
python3 experiments/rook_symmetry_breaking_2026-08-04/run.py --cpu 2
```

`results.csv` is append-only and synchronized about every ten seconds. Restarting the same command reads the file and skips every
completed coefficient vector, so an interrupted experiment can continue without losing its recorded work.

## Run result

The search ran on CPU 2 from 2026-08-04 00:52:00 UTC through 09:55:21 UTC and stopped cleanly when the user returned. The saved CSV
contains 57,308 consecutive, unique matrices: the baseline plus 57,307 perturbations. Every evaluation completed successfully; there
were no timeouts, errors, or whole-matrix safe fallbacks.

No perturbation exceeded the 15,120-ESS baseline. Of the perturbations, 14,242 tied it; the best strictly lower result was 15,024 ESS.
As specified above, these are fast-mode exploratory counts and were not separately safe-verified.
