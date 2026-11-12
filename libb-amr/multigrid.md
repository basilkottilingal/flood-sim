# Multigrid: a plan

*Planned 2026-11-11. A first version is done (2026-11-12), smaller than
this plan: the solver is `liblinear-algebra/multigrid.[ch]`.*

**Done, as decided then:**
- `foreach_level (l)` and `foreach_level_or_leaf (l)` (A), on the tree's
  lists by level (`tree_caches ()`), in `block-loops.h`.
- `ghosts_update_level (s, l)` (C, in part): the ghosts of the grid at depth
  l, the parents coarser than l restricted from that grid first. The
  restriction is needed after all: a coarse leaf of the grid reads its
  neighbours of its own level, which may be parents. Called explicitly by
  the solver: the level loops give `ghosts_by_caller`, and cwrap writes no
  update or mark around them.
- The residual goes down by the scalar's restriction (`tree_restrict`), the
  correction up by its prolongation (`block_prolong`): no `coarse ()` in a
  level loop, and none of B, D or E needed yet.
- The root (question 4): settled, as the tree always keeps level 1
  (`tree ()`, `adapt_wavelet ()`).

**Still open:** stale bits per level and lazy updates by cwrap (C), then
`coarse ()` in `foreach ()` (B), the `filling` marker (D) and the reach per
loop (E). The plan below is as it was written.

The aim: geometric multigrid on the block tree, written with cwrap's stencil
syntax, as in Basilisk: `coarse ()` read from loops over levels, and in
`foreach ()` itself.

## What a V-cycle needs

1. **Smoothing on one level:** a relaxation on every cell of level `l`,
   parents included, reading its neighbours on that level: ordinary offsets,
   `s[1]`, `s[-1]`, ...
2. **Restriction of the residual,** level `l` to `l - 1`: a loop over the
   parents, with `foreach_child ()` (as `foreach_restriction` does now).
3. **Prolongation of the correction,** level `l - 1` to `l`: a loop over
   level `l` reading `coarse (dp)`, `coarse (dp, ±1)` and `child.x`.

So the addition is **loops over levels**, in which both same-level offsets and
`coarse ()` are allowed. `coarse ()` in `foreach ()` (the leaves) is the case
where each leaf reads its parent.

## A. Loops over levels

Declared with `@foreach` as the others:

| Loop | Points | Use |
|---|---|---|
| `foreach_level (l)` | every node at level `l`, leaves and parents | smoothing; prolonging the correction |
| `foreach_level_or_leaf (l)` | level `l`, and the leaves coarser than `l` | the grid at depth `l` of an adaptive tree, as Basilisk's |
| `foreach_coarse_level (l)` | the parents at level `l` | restriction, with `foreach_child ()` |

Built on the `Levels` lists the ghost core already has (`levels_new ()` in
`block-ghosts-impl.h`). Their block loop gives `coarse_sblocks`, the parent
node's blocks, so that `coarse ()` works in them.

## B. `coarse ()` in `foreach ()`

- `foreach_leaf` also gives the parent's blocks, and `foreach_cell` the coarse
  index `cc` and the half `sgn` of each cell, from the leaf's place in its
  parent: one lookup per leaf.
- The root of a single-block tree has no parent: `coarse ()` there needs a
  guard (see the questions).

## C. Ghost state by level (the real design)

Today: one bit per scalar, "the leaves' ghosts may be stale"
(`ghosts_dirty`, `block-ghosts.h`). `coarse ()` needs more:

- `coarse (s)` with no offset reads the parent's interior: valid once the
  leaves are restricted, which `ghosts_update ()` does first.
- `coarse (s, ±1)` may reach the parent's ghosts, which only
  `ghosts_update_levels ()` fills.
- A loop over level `l` needs the ghosts of level `l` alone, made **without
  restriction**: a correction on level `l` is a field of its own, which the
  restriction of finer cells must not overwrite.

Proposed:

- **Stale bits per level:** 64 scalars by at most 21 levels, a few hundred
  words.
- **`ghosts_update_level (s, l)`:** the ghosts of level `l` only, by copy from
  that level, prolongation from `l - 1` and boundary conditions; no
  restriction.
- **Marking:**
  - a write in `foreach_level (l)` marks level `l` stale;
  - a write in `foreach ()` marks the leaves stale, and the parents too
    (their restriction is out of date);
  - a read at an offset in `foreach_level (l)` updates level `l` first, if
    stale;
  - `coarse ()` reads need level `l - 1` (or the parents) up to date.

## D. Telling filling loops from level loops

Now, a loop whose cell loop gives `cc` or `fc` is taken as a filling loop:
no offset on its own level, no ghost tracked. That is right for
`foreach_prolongation` and `foreach_restriction`, which run inside the ghost
update, and wrong for `foreach_level (l)`, which needs both.

Proposed, keeping the convention of names (no new keyword in `@foreach`):
filling loops' cell loops give a name, `filling`, and only loops that give it
get the restrictive rules.

## E. The reach of `coarse ()`, per loop

In a prolongation only the coarse block's first ghost layer is valid
(`BAMR_COARSE_REACH = 1`). In a level loop, after a level update, all
`GHOSTS` layers are. So the limit comes from the loop: its cell loop gives a
constant `coarse_reach` (an `enum` inside it), and the `_Static_assert` cwrap
already writes checks against that name.

## Order

1. **D and E**, small and keeping behaviour: the `filling` marker, the reach
   per loop. Every test as before, bit for bit.
2. **C:** `ghosts_update_level (s, l)` and stale bits per level, tested
   against `ghosts_update_levels ()` on one level.
3. **A:** the loops over levels, with `coarse ()` in them.
4. **B:** `coarse ()` in `foreach ()`, with the parents kept up to date.
5. **Check:** a Poisson multigrid solver in a libb-amr test, written with
   these loops on an adaptive tree, reducing the residual by the expected
   factor per V-cycle: the proof the semantics are right.

## Questions to settle first

1. Basilisk's three loops over levels, with its names and meanings?
2. Stale bits per level, and `ghosts_update_level (s, l)` with no restriction,
   as the basis?
3. Filling loops marked by a name their cell loop gives (`filling`)?
4. `coarse ()` at the root of a single-block tree: an error at run time, or a
   rule that it needs one level of refinement at least?
