# Sokoban Move-Optimal Solver

A parallel A* solver that computes provably minimum-move solutions for a
Sokoban variant in which every box must end on a goal square and the keeper
must also finish standing on a goal square.

Everything lives in one file: [`solver.cpp`](solver.cpp).

## Building

```bash
g++ -O3 -march=native -std=c++17 -pthread solver.cpp -o solver
```

`-O2` also works (~5% slower). Requires GCC or Clang on a 64-bit target
(`__int128` and GCC builtins are used). Transparent-hugepage support is used
when available (Linux) but is not required.

## Running

```bash
./solver <boardId> [timeoutSec] [--solution|-s] [--no-dominance] [dbg]
```

| Argument     | Meaning                                                        | Default |
|--------------|----------------------------------------------------------------|---------|
| `boardId`    | one of the built-in boards `p1` … `p15`                        | `p15`   |
| `timeoutSec` | wall-clock limit for the search                                 | `600`   |
| `--solution` | also reconstruct and print one optimal move sequence (`-s`, `sol`) | off |
| `--no-dominance` | disable keeper-dominance pruning (`nodom`), for comparison | on |
| `dbg`        | print per-box heuristic details and exit (debugging aid)       | off     |

Positional arguments are `boardId` then `timeoutSec`. Flags may appear anywhere.

```
$ ./solver p15 600 --solution
[p15] boxes=8 goals=9 live=32 h0=32 timeout=600s
  ...expanded 4120564, f=87 gen=17269394 pruned=3658143 best=-1 0s      <- progress every ~4M expansions
  ...
RESULT p15 optimal=184 expanded=21399195 pruned=27691734 generated=76487677 seconds=2.0 threads=24
SOLUTION p15 moves=184 RdRdRdRRUrDDDrdLL...
```

`RESULT` reports `optimal=<N>` on success, or `UNSOLVABLE`, `NONE` (no solution
found within the search's f-limit, see Limits), `TIMEOUT` (with the proven
lower bound reached), or `UNSUPPORTED` (board exceeds a limit, see below). Two
abnormal exits (code 2, message on stderr) are also possible: `out of memory`
if the up-front hash allocation fails, and `hash full` if a board generates
more unique states than the fixed-size table holds. The remedy for the latter
is raising the `1 << 27` in `run()` and rebuilding.

In a `SOLUTION` string, lowercase `u d l r` are keeper walking moves and
uppercase `U D L R` are pushes. Its length always equals the reported
optimum.

### Built-in boards

Boards are digit-string grids hardcoded in the `BOARDS` map:
`0` floor, `1` wall, `2` box, `3` keeper, `4` goal, `5` box-on-goal,
`6` keeper-on-goal. To add a board, add an entry to the map. But note that
`setup()` performs no validation, so a new board must satisfy: all rows the
same length (a short row silently parses as floor and yields a confidently
wrong answer), exactly one keeper (`3`/`6`. A keeperless board crashes), and
digits `0` to `6` only (anything else is silently floor). The grid border acts as
an implicit wall, so boards need not be wall-enclosed. Surplus goals are
allowed (they may stay empty). One goal must be left free for the keeper.

Verified optima: p1=6, p2=15, p3=13, p4=17, p5=12, p6=13, p7=47, p8=22, p9=34,
p10=59, p11=51, p12=41, p13=78, p14=26, p15=184. All solve in well under a
second except p15 (~2 s / ~1.5 GB on a 12-core desktop).

## Algorithm

### State space

A state is `(box bitmask, keeper cell)`. Boxes are a bitmask over the live
cells. Squares from which a box can still reach some goal (computed by
reverse BFS from every goal). Everything else is provably dead. Boards with
≤ 32 live cells use a `uint32` mask, and the whole state key packs into 40 bits
(`mask<<8 | keeper`), so one hash slot is a single 64-bit word holding key and
value together. Boards with 33–128 live cells switch to a 128-bit mask with a
wider hash. The entire search is templated on the mask type, so small boards
pay nothing for the large case.

Limits: ≤ 255 cells total (8-bit keeper index), ≤ 128 live cells, and
`f = g + h < 6000` (`MAXF`). A board whose initial heuristic already reaches
6000 is reported `UNSUPPORTED`, and a (pathological) board whose optimum is
≥ 6000 moves would report `NONE` rather than being found.

Memory: the best-g hash is allocated eagerly and sized by live-cell count.
32 MiB for boards with < 24 live cells, 1 GiB for 24–32 live cells (p7,
p15), and an estimate-based size (up to ~800 MiB) on the wide path. Machines
with less free RAM than that get a clean `out of memory` abort before the
search starts. The fast path also allocates the keeper-dominance table at a
quarter of the hash's slots (256 MiB for p15). It is optional: if that
allocation fails, the search runs without dominance pruning.

### Search: move-optimal A* over macro pushes

Successors are macro pushes: "walk to a push position, push once", with
edge cost (shortest keeper walk) + 1. `g` therefore counts real keeper
moves, so the answer is the true minimum-move solution, while the state space
stays the compact `(boxes, keeper)`. An explicit per-step walk graph is never
searched.

The open list is a Dial (bucket) queue indexed by integer
`f = g + h`. The heuristic `h` is the sum over boxes of the precomputed push
distance to the nearest goal (`DMIN`). It is admissible and consistent for
these macro edges, which makes the search safe to run in parallel: the first
time a state is claimed for expansion it already has its optimal `g`, so states
close on first pop with no reopening. At expansion `h(parent) = f − g` exactly,
so a successor's `f` is an O(1) delta (`− DMIN[from] + DMIN[to]`) rather than a
loop over boxes.

When a state with all boxes on goals is expanded, the finishing cost is the
shortest walk to a free goal. The best such total is the incumbent `best`.
Levels are processed in increasing `f`, so once `f ≥ best` every cheaper
possibility has been exhausted and `best` is the proven optimum.

### Pruning

* Dead squares. A push onto a cell that cannot reach any goal is never
generated (from the reverse-BFS table).
* Freeze deadlock. A push that leaves a box immovably wedged (against
walls / dead squares / other frozen boxes, recursively, on both axes) is
discarded unless the box lands on a goal.
* Incumbent gate. Successors with `f ≥ best` are never queued (they can
never be expanded), though their `g` is still recorded for reconstruction.
* Boards where a starting box sits on a dead square are reported `UNSOLVABLE`
immediately.

#### Keeper-dominance pruning

Many states share a box layout and differ only in where the keeper stands.
p15's 49 M states cover just 6 M box layouts, about 8 keeper cells each.
A popped state `(B, kb, gb)` is *dominated* by an already-expanded state
`(B, ka, ga)` with the same boxes when the keeper can walk from `ka` to `kb`
and `ga + walk(ka, kb) ≤ gb`. Every push the dominated state could make, the
dominating one can make at least as cheaply by walking to `kb` first (triangle
inequality on walks), and the same holds for the finishing walk onto a free
goal. So expanding it can never lower any `g`, and it is skipped. It stays
closed in the hash, so solution reconstruction is unaffected: a g-tight
predecessor always has an exact `g`.

Optimality holds under parallelism too. States with the same boxes have the
same `h`, so a dominator always has a strictly smaller `f` and was expanded at
an earlier level. Within one level, two distinct variants can never dominate
each other, since that would need a zero-length walk.

The implementation (`DomHash`, fast path only) is a lock-free side table keyed
by the box mask. It remembers the first 4 variants expanded per layout, i.e.
the lowest-`g` ones, because a layout's variants come out in increasing `g`.
Its slot is prefetched alongside the main hash slot. The keeper BFS that every
expansion runs anyway does the test (`keeperBFSDom`). It stops as soon as it
reaches a stored variant's cell within that variant's walk budget `g − ga`, so
pruned states also skip most of their BFS. The table is purely an
accelerator: when it is full or a probe runs long, the check is skipped.

| p15, 12-core / 24-thread desktop | expanded | pruned | 1 thread | 12 threads | 24 threads | peak RSS |
|---|---|---|---|---|---|---|
| `--no-dominance` | 49,090,929 | 0 | 20.7 s | 4.2 s | 3.1 s | 1.21 GB |
| default | 21,399,255 | 27,691,674 | 14.3 s | 2.7 s | 2.0 s | 1.47 GB |

(Times include the bitboard keeper BFS below. Before it: 26.9 / 5.8 / 4.2 s
without dominance pruning, 17.6 / 3.8 / 2.9 s with it.)

The generated-state count is identical with and without the pruning (76.49 M
single-threaded), which confirms that every skipped expansion was redundant.
Checking against *every* expanded variant (measured offline) would prune at
most 57.5% of states, vs. the 56% the 4-variant table achieves.

### Parallelism

Each `f`-level is expanded in waves by a persistent pool of worker threads
(up to 32, sized to the hardware). All shared search state lives in one
open-addressed concurrent hash (linear probing, fixed size, no rehash) that
maps a state to its best-known `g` plus a closed bit:

* `improve(state, g)`. CAS-insert or CAS-lower the stored `g`.
* `claim(state, g)`. Atomically verify the entry still has exactly this `g`
and is unclosed, and set the closed bit. The winner expands the state,
duplicates and stale queue entries fail cheaply.

On the fast path each slot is one atomic `uint64` (key and value in the same
word, pure relaxed CAS. Genuinely lock-free). On the wide path a slot is an
atomic control word (`EMPTY → WRITING → READY|value`) plus side arrays for the
128-bit key, published with release/acquire ordering. Readers briefly spin on
a slot mid-publication, so it is wait-free in practice but not formally
lock-free.

Workers accumulate successors in per-`(thread, f)` chunk buffers that are
bulk-moved into the global buckets after each wave, and keep private counters.
The hot loop performs no shared-cacheline writes besides the hash itself.

### Solution reconstruction (`--solution`)

The search stores only `g` values. No parent pointers. So the path is
recovered afterwards from the completed hash: find a finishing state whose
`g + finishing-walk = optimum`, then repeatedly invert the last push (the
keeper stands on the cell the box vacated) to a predecessor whose stored `g`
is tight: `g_pred + walk + 1 = g`. Consistency guarantees every stored `g`
on such a chain is the true distance from the start, so greedy tight-predecessor
steps provably reach the initial state with no backtracking. Reversing the
collected segments and appending the finishing walk yields a complete optimal
move string (verified to replay correctly on all 15 boards).

### Performance notes

The search is memory-latency bound, so the implementation leans on:
2 MB-aligned, `MADV_HUGEPAGE`-backed hash memory (the table is hit at random.
4 KB pages would thrash the TLB). Software prefetch of upcoming hash slots.
Bucket entries packed into single `uint64`s on the fast path. And the chunked
merge above. The 1 GiB table (2²⁷ slots) used for large fast-path boards was
measured fastest for the hardest built-in board. p15 pops ~49 M states (21 M
expanded, 28 M dominance-pruned) in ~2 s on a 12-core (24-thread) desktop.
Larger tables measurably hurt (TLB spread), smaller ones lengthen probe chains.

After dominance pruning, the per-state keeper BFS was the largest single cost
(47% of single-thread time on p15). A queue BFS spends most of that on
mispredicted branches: for every neighbor it asks "wall? box? already
visited?", and the answers are data-dependent. `keeperBFSbb` replaces it with
a bitboard flood fill. The board is laid out as a 256-bit grid with each row
padded by one never-free column, so a cell's neighbors are one-bit and one-row
shifts, and each distance layer is a few dozen branch-free 64-bit operations
(`next = neighbors(frontier) & free & ~visited`). Distances are written only for
the cells the search reads afterwards: push spots next to boxes, goals, and the
dominance targets. That keeps `kReach()`/`kdist` working unchanged and gives
exactly the same search (identical single-thread counts). p15: 17.6 s → 14.3 s
on one thread, 2.9 s → ~2 s on 24. Boards whose padded grid exceeds 255 bits
fall back to the queue BFS.
