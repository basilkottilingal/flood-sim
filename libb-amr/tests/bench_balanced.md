# memindex vs voxel-memindex on 2:1 balanced trees

`sh tests/bench_balanced.sh` (2026-11-09). Intel i7-8650U, 8 GB (about
3 GB free), Linux 6.8, gcc 11.4 -O2, THP `[madvise]`, memindex with
MADV_NOHUGEPAGE. Second of two runs; the first agreed within ~10%.

The trees are built as `tree.c` builds them: cells touching a front
(circle/sphere of radius 0.35, or a meandering river two finest cells wide)
are refined, a family at a time, after their 3^d - 1 neighbours exist
(`tree_refine ()`, `covering ()`). Every level has its own array; nodes are
inserted in pre-order (Morton order, parents first).

- B/node: (RSS + page tables) after building, per node of every level.
- dfs: depth-first traversal from the root in Morton order, each child
  looked up (`children_of ()`), its record read; ns per node.
- list: the nodes in pre-order, each looked up (`index_node ()` over a
  list, as `ghosts_update ()` does); ns per node.
- nbr: for every leaf, its 3^d - 1 neighbours at its level, each falling
  back to the covering coarser node (`block_ghosts ()`); ns per leaf, and
  per lookup. About 20% of the lookups fall back to a coarser level.
- +cur: voxel with one `VMemCursor` per level.

| case | impl | nodes (leaves) | B/node | build ns | dfs ns | dfs+cur ns | list ns | list+cur ns | nbr ns/leaf | nbr+cur ns/leaf | ns/lookup (+cur) | free ns |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2D dense depth 10 | mem | 1398101 (1048576) | 9 | 14.5 | 5.5 | - | 3.7 | - | 38 | - | 4.8 (-) | 12.6 |
| 2D dense depth 10 | voxel | 1398101 (1048576) | 8 | 12.5 | 8.1 | 4.7 | 7.7 | 3.0 | 61 | 51.8 | 7.7 (6.5) | 11.9 |
| 2D shell depth 15 | mem | 638149 (478612) | 615 | 568.6 | 20.4 | - | 14.2 | - | 106 | - | 10.5 (-) | 420.0 |
| 2D shell depth 15 | voxel | 638149 (478612) | 25 | 22.6 | 13.8 | 7.1 | 11.1 | 4.7 | 101 | 76.8 | 10.0 (7.6) | 28.9 |
| 2D river depth 15 | mem | 1286501 (964876) | 234 | 261.8 | 11.7 | - | 8.2 | - | 71 | - | 7.2 (-) | 238.3 |
| 2D river depth 15 | voxel | 1286501 (964876) | 23 | 21.6 | 12.5 | 6.6 | 8.9 | 4.0 | 92 | 72.6 | 9.3 (7.3) | 21.0 |
| 3D dense depth 6 | mem | 299593 (262144) | 77 | 97.8 | 32.9 | - | 8.4 | - | 272 | - | 10.8 (-) | 37.0 |
| 3D dense depth 6 | voxel | 299593 (262144) | 8 | 13.8 | 25.2 | 7.3 | 8.7 | 4.4 | 310 | 274.2 | 12.3 (10.9) | 12.3 |
| 3D shell depth 8 | mem | 576009 (504008) | 263 | 362.0 | 45.9 | - | 12.8 | - | 371 | - | 11.6 (-) | 177.1 |
| 3D shell depth 8 | voxel | 576009 (504008) | 15 | 19.2 | 27.9 | 8.7 | 10.9 | 5.5 | 431 | 380.1 | 13.5 (11.9) | 17.9 |
| 3D river depth 10 | mem | 374513 (327699) | 446 | 548.2 | 38.9 | - | 9.6 | - | 321 | - | 9.5 (-) | 327.7 |
| 3D river depth 10 | voxel | 374513 (327699) | 18 | 21.9 | 28.7 | 8.6 | 12.1 | 5.0 | 472 | 393.2 | 14.0 (11.7) | 20.6 |
| 3D shell depth 9 | mem | 2287433 (2001504) | 252 | 317.5 | 44.5 | - | 11.9 | - | 310 | - | 9.7 (-) | 209.3 |
| 3D shell depth 9 | voxel | 2287433 (2001504) | 14 | 19.3 | 28.3 | 8.5 | 12.5 | 5.0 | 462 | 392.1 | 14.4 (12.3) | 20.4 |

## Reading it

- **Morton traversal (dfs).** Without a cursor, voxel is slower than
  memindex only on the uniform 2D tree (8.1 against 5.5 ns per node). It is
  about equal on the 2D river, and faster on the 2D circle (13.8 against
  20.4) and in 3D (25-29 against 33-46). With a cursor per level, voxel is
  faster in every case:
  - 1.2x on the uniform 2D tree;
  - 1.8-2.9x on the 2D fronts;
  - 4.5-5.3x in 3D (8.5 against 44.5 ns per node on the sphere at depth 9).

  memindex reads a family's children from 2 rows (2D) or 4 rows (3D), each
  on a different page. voxel finds all 2^d children in one bottom node.
- **Over a level list (list).** Without a cursor, memindex is 2x faster on
  the uniform 2D tree and within about 25% either way elsewhere. With a
  cursor, voxel is 1.2-3x faster in every case.
- **Neighbours (nbr).** This is where memindex keeps an edge. With cursors,
  voxel compares as follows:
  - 28% faster on the 2D circle;
  - equal (within 2%) on the 2D river, the uniform 3D tree and the 3D
    sphere at depth 8;
  - 22-26% slower on the 3D river and the 3D sphere at depth 9;
  - 36% slower on the uniform 2D tree.

  Without cursors, voxel is up to 1.6x slower (equal on the 2D circle). The
  cursor helps less here for two reasons. The box around a leaf on the edge
  of its 4x4x4 (3D) or 8x8 (2D) group spills into the next bottom node. And
  a fifth of the lookups fall back to a coarser level, each through that
  level's own cursor.
- **What the difference costs.** The largest gap is about 80 ns per leaf,
  for its whole box of 26 neighbours (3D). `block_ghosts ()` then copies
  the block's ghost layers from those neighbours: thousands of cells, an
  estimated few microseconds per scalar for 16^3 blocks. So the gap would
  be around a percent of a ghost update or less. This is an estimate, not
  a measurement; timing `ghosts_update ()` with voxel wired in would
  confirm it.
- **Memory, build, free.** On every non-uniform tree, voxel takes 10-25x
  less memory (14-25 against 234-615 B per node). It builds 12-25x faster
  and frees 10-16x faster. On the uniform trees both are close, except in
  3D, where voxel uses 10x less memory, builds 7x faster and frees 3x
  faster.
