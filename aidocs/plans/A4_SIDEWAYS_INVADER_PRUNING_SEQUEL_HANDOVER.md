# A4 Sideways Invader Pruning: Sequel Handover

## Start Here

Continue in `/home/reinhard/projects/fracessa` on branch `main`. At handover time, `HEAD` is `41dfa8ff` and the worktree contains
many intentional uncommitted changes from concurrent work. Preserve them. Read `AGENTS.md` before acting, use Ponytail for code work,
and obtain Reinhard's approval before changing C++ or Python source.

The authoritative background is [A4_SIDEWAYS_INVADER_PRUNING.md](A4_SIDEWAYS_INVADER_PRUNING.md). Model A4 is isolated under
`models/a4`; production behavior is unchanged.

## Current A4 Rule

For an internally valid exact candidate on support `S`, with common payoff `u`, define the complete strict-invader set

\[
V(x)=\{j\notin S:(Ax)_j>u\}.
\]

Every ESS support containing `S` must intersect `V(x)`. A4 therefore rejects later supports `T` satisfying

\[
S\subsetneq T\subseteq N\setminus V(x).
\]

Current restrictions:

- sideways pruning is enabled only for non-circular matrices;
- rules are learned only from singleton and pair supports, `|S| <= 2`;
- `V(x)` is always complete and nonempty; it is never truncated;
- there is no current restriction on `|V|` and no usefulness threshold;
- there is no independent rule-count cap;
- the cardinality restriction permits at most
  \(n+\binom n2=n(n+1)/2\) rules, hence `O(n^2)` rules;
- collecting all invaders over all eligible pairs can nevertheless cost `O(n^3)` exact payoff comparisons.

Do not retain only part of `V`: that would strengthen the clause and can make pruning incorrect. Skipping an entire optional clause is
safe because it only gives up pruning.

The restrictions are hard-coded in two places:

- `models/a4/src/exact_candidate_solver.cpp`: after the first strict invader, supports larger than two return immediately;
- `models/a4/src/fracessa.cpp`: only singleton and pair `invader_interval` results are passed to the generator.

## Retained Implementation

Relevant files:

- `models/a4/include/fracessa/support_generator_non_circular.hpp`
- `models/a4/include/fracessa/exact_candidate_solver.hpp`
- `models/a4/include/fracessa/fracessa.hpp`
- `models/a4/src/exact_candidate_solver.cpp`
- `models/a4/src/fracessa.cpp`
- `models/a4/verify.py`
- `models/a4/README.md`

The ordinary generator indexes each clause by the lowest required bit of `S`. During descending-bit DFS, a clause becomes active when
the partial support first contains all of `S`. Active clauses are checked before a completed support is emitted. Pending clauses are
activated only at the next cardinality, so callbacks cannot invalidate active pointers.

The implementation has separate one-word and multiword paths. The lowest-bit buckets, pending staging, and active-rule pointers are
already lean. Do not replace them with SAT, a BDD, a generalized constraint engine, hashes, or watched literals without new evidence;
those broader designs were tested and were dramatically slower.

Circular matrices deliberately retain A2's bracelet/orbit generator. They do not apply sideways clauses. The current solver still
collects complete invader sets for rejected circular singleton and pair supports, but those sets are never consumed. Returning after
the first invader for circular matrices is therefore exactly equivalent; only the size of the speed gain is empirical.

## Verification and Benchmark State

The maintained A4 verification matched A2 on 29 representative ordinary, circular, affine-circular, multiword, and database cases:

```bash
CC=/usr/lib64/ccache/cc CXX=/usr/lib64/ccache/c++ cmake -S models/a4 -B models/a4/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLINT_INCLUDE_DIR="$PWD/.local/flint-3.6.0/include" \
  -DFLINT_LIB="$PWD/.local/flint-3.6.0/lib/libflint.so"
cmake --build models/a4/build --parallel
/home/reinhard/.local/bin/python3.12 models/a4/verify.py --target-seconds 0
```

The saved timing session is `a2_a4_invader_20260822` in the ignored database `experiments/model_timings.sqlite3`. It selected the 1,071
stored matrices with `0 < safe_calibration_ns < 100,000,000,000`, excluding dimension two. The target was 0.5 seconds, timeout 30
seconds, scheduler CPU 1, and workers CPUs 2 through 9. A2 rows came from the unchanged stored baseline; only A4 was rerun.

Results:

- A4: 1,070 correct, zero mismatches, one timeout (`2493`);
- A2: 1,069 correct, zero mismatches, two timeouts (`2453`, `2493`);
- matrix `2453`: A2 timed out at 30 seconds; A4 finished in 0.169 ms;
- 1,069 paired successful rows: median A4/A2 `1.025651634`, or 2.565% slower;
- paired outcomes: 308 faster, 685 slower, 76 exact timing ties;
- 176 paired matrices improved by at least 5%, 70 by at least 2x, and 17 by at least 10x;
- largest absolute regression: matrix `211`, +0.330 seconds;
- largest ratio regression: matrix `2169`, 77.98x, from 0.374 ms to 29.196 ms.

The earlier size-class table mixed circular and non-circular matrices. Correct interpretation:

| Group | Rows | Median slowdown |
|---|---:|---:|
| All non-circular | 714 | 7.1128% |
| Small, all symmetries | 172 | 7.1128% |
| Small, non-circular | 119 | 13.9548% |
| Medium, non-circular | 214 | 10.2165% |
| Large, non-circular | 198 | 7.2644% |
| Super-large, non-circular | 183 | 4.9333% |

The equality between `all non-circular` and `small, all symmetries` is real but coincidental: both medians use the same two central
measurements, matrix IDs `2049` and `2280`. Overall medians are not averages of subgroup medians.

## Correctness Concern Found During Review

The generator inherited A2's shortcut: if one cardinality emits no support, stop all larger cardinalities. That is valid when every
pruning rule is an ordinary upward cone. It is not valid for arbitrary sideways clauses, because adding an invader can make a larger
support valid again.

Direct generator counterexample:

```text
{0} requires 1
{1} requires 2
{2} requires 0
```

Every pair violates one clause, but `{0,1,2}` satisfies all clauses. The current A4 generator emits the singletons, emits no pair,
then stops before the valid triple. This reproduces directly against the generator API.

No actual matrix mismatch is currently known: the complete 1,071-matrix run, the maintained verification corpus, 200,000 random
six-dimensional matrices, and one million random seven-dimensional matrices all matched A2. Nevertheless, the generic generator
shortcut is unjustified once sideways clauses exist unless an additional theorem about payoff-derived clauses is supplied.

Smallest safe repair: preserve the empty-cardinality shortcut only while no invader clause has ever been installed. Apply the same
guard to the one-word and multiword generators, and add a regression check using the three-rule example above. Do this before further
performance tuning.

## Ranked Next Work

After the correctness repair, make one change at a time and compare against the saved A2 rows and current A4 session.

1. **Skip provably useless rules.** Let

   \[
   f=n-|S|-|V|.
   \]

   A rule can reject at most `2^f - 1` later supports. When `f=0`, it rejects none and must not be stored. This is exact, not a
   heuristic.

2. **Stop collecting invaders when they cannot be consumed.** Pass an explicit `collect_invaders` policy into the exact solver. It
   should be true only for non-circular singleton and pair supports. This centralizes the duplicated cardinality policy and restores
   immediate rejection for circular supports.

3. **Specialize singleton comparisons.** For `S={i}`, probabilities are exactly one and the invader test is simply
   `A(j,i) > A(i,i)`. Avoid general FLINT multiplication by one.

4. **Make activation cardinality-aware.** When a clause becomes active, let `r` be the remaining number of support elements to choose.
   If fewer than `r` future non-invaders remain, every completion must choose an invader. The clause is guaranteed satisfied and need
   not remain active. This is exact and costs one popcount on the one-word path.

5. **Only then benchmark usefulness thresholds.** `f=1` rejects one later support, `f=2` at most three, and so on. Skipping whole
   clauses below a threshold is correctness-safe but performance-empirical. Compare thresholds rather than guessing.

Do not optimize active-rule removal or the multiword representation first. The 15 paired benchmark matrices above dimension 64
already improved on average, while the clearest waste is complete invader collection and weak clauses on small non-circular games.

## Required Decision Boundary

No C++ or Python source was changed during the final review. The sequel session must present the exact intended source patch and obtain
Reinhard's approval before implementing it.
