# memindex vs voxel-memindex vs nano-vdb on 2:1 balanced trees

`sh tests/bench_balanced.sh`. Intel i7-8650U laptop, 8 GB, Linux 6.8, gcc 11.4 `-O2 -mpopcnt`, THP
`[madvise]`, memindex with `MADV_NOHUGEPAGE`. Run on 2026-11-10 with a
browser busy in the background (load about 1.5), so absolute times are
slower and noisier than the earlier runs; compare implementations within a
case, where they ran back to back. Cache and dTLB miss counters were not
available (`perf_event_paranoid` = 4).

Implementations: `mem` (memindex), `voxel` (voxel-memindex, 64-way nodes),
`nvdb` (nano-vdb, dense VDB nodes), `nvdb-compact` (nano-vdb with masks and
popcount rank, `-DNVDB_COMPACT`). Cursors: voxel's bottom node, nano-vdb's
whole-path accessor. Every implementation returned the same values as
memindex (checksums).

Second of two runs; the first agreed within about 15%.

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
- +cur: one cursor (voxel) or accessor (nano-vdb) per level.

| case | impl | nodes (leaves) | B/node | build ns | dfs ns | dfs+cur ns | list ns | list+cur ns | nbr ns/leaf | nbr+cur ns/leaf | ns/lookup (+cur) | free ns |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2D dense depth 10 | mem | 1398101 (1048576) | 9 | 17.2 | 18.0 | - | 7.4 | - | 109 | - | 13.6 (-) | 12.8 |
| 2D dense depth 10 | voxel | 1398101 (1048576) | 8 | 15.7 | 17.2 | 6.7 | 8.7 | 3.8 | 140 | 67.2 | 17.5 (8.4) | 12.9 |
| 2D dense depth 10 | nvdb | 1398101 (1048576) | 9 | 15.0 | 18.8 | 10.5 | 6.8 | 5.8 | 151 | 76.6 | 18.9 (9.6) | 8.5 |
| 2D dense depth 10 | nvdb-compact | 1398101 (1048576) | 9 | 38.9 | 34.5 | 9.1 | 16.5 | 6.9 | 214 | 76.3 | 26.8 (9.6) | 21.4 |
| 2D shell depth 15 | mem | 638149 (478612) | 615 | 736.3 | 54.3 | - | 19.3 | - | 228 | - | 22.6 (-) | 427.3 |
| 2D shell depth 15 | voxel | 638149 (478612) | 25 | 26.3 | 23.7 | 11.6 | 13.1 | 5.7 | 221 | 135.1 | 21.9 (13.4) | 40.2 |
| 2D shell depth 15 | nvdb | 638149 (478612) | 30 | 30.3 | 27.4 | 9.8 | 11.2 | 8.2 | 248 | 102.8 | 24.6 (10.2) | 14.1 |
| 2D shell depth 15 | nvdb-compact | 638149 (478612) | 13 | 55.8 | 33.5 | 11.4 | 18.3 | 8.4 | 268 | 107.2 | 26.6 (10.6) | 32.2 |
| 2D river depth 15 | mem | 1286501 (964876) | 234 | 310.7 | 46.1 | - | 21.3 | - | 207 | - | 20.8 (-) | 315.6 |
| 2D river depth 15 | voxel | 1286501 (964876) | 23 | 33.3 | 24.3 | 10.3 | 12.4 | 5.7 | 213 | 95.2 | 21.4 (9.6) | 25.4 |
| 2D river depth 15 | nvdb | 1286501 (964876) | 28 | 48.6 | 30.7 | 12.5 | 16.9 | 6.3 | 226 | 99.8 | 22.7 (10.0) | 16.1 |
| 2D river depth 15 | nvdb-compact | 1286501 (964876) | 13 | 40.7 | 32.3 | 10.5 | 20.4 | 8.0 | 264 | 97.5 | 26.6 (9.8) | 39.3 |
| 3D dense depth 6 | mem | 299593 (262144) | 77 | 118.7 | 37.9 | - | 9.6 | - | 531 | - | 21.1 (-) | 45.0 |
| 3D dense depth 6 | voxel | 299593 (262144) | 8 | 15.4 | 25.7 | 19.6 | 9.1 | 5.7 | 537 | 460.8 | 21.3 (18.3) | 12.5 |
| 3D dense depth 6 | nvdb | 299593 (262144) | 9 | 13.1 | 26.1 | 21.5 | 7.3 | 5.9 | 535 | 428.4 | 21.2 (17.0) | 9.7 |
| 3D dense depth 6 | nvdb-compact | 299593 (262144) | 9 | 37.9 | 31.1 | 22.8 | 11.1 | 6.3 | 661 | 483.9 | 26.3 (19.2) | 33.4 |
| 3D shell depth 8 | mem | 576009 (504008) | 263 | 343.4 | 50.4 | - | 15.6 | - | 653 | - | 20.5 (-) | 168.4 |
| 3D shell depth 8 | voxel | 576009 (504008) | 15 | 22.3 | 28.0 | 20.8 | 11.5 | 6.8 | 712 | 600.9 | 22.3 (18.8) | 18.4 |
| 3D shell depth 8 | nvdb | 576009 (504008) | 27 | 22.2 | 33.3 | 25.2 | 12.0 | 10.4 | 679 | 560.5 | 21.3 (17.6) | 17.5 |
| 3D shell depth 8 | nvdb-compact | 576009 (504008) | 12 | 41.0 | 33.2 | 23.7 | 14.7 | 8.4 | 806 | 597.3 | 25.2 (18.7) | 34.1 |
| 3D river depth 10 | mem | 374513 (327699) | 446 | 560.1 | 45.3 | - | 12.8 | - | 657 | - | 19.5 (-) | 291.7 |
| 3D river depth 10 | voxel | 374513 (327699) | 18 | 21.4 | 29.1 | 21.0 | 14.4 | 6.5 | 788 | 623.5 | 23.4 (18.5) | 22.2 |
| 3D river depth 10 | nvdb | 374513 (327699) | 44 | 29.0 | 31.9 | 26.7 | 12.2 | 8.9 | 716 | 587.4 | 21.3 (17.4) | 11.5 |
| 3D river depth 10 | nvdb-compact | 374513 (327699) | 14 | 41.7 | 32.6 | 24.5 | 11.5 | 7.1 | 837 | 611.9 | 24.9 (18.2) | 33.3 |
| 3D shell depth 9 | mem | 2287433 (2001504) | 252 | 329.1 | 52.8 | - | 13.3 | - | 629 | - | 19.7 (-) | 215.8 |
| 3D shell depth 9 | voxel | 2287433 (2001504) | 14 | 20.1 | 29.3 | 20.7 | 13.6 | 6.4 | 747 | 614.4 | 23.4 (19.2) | 20.0 |
| 3D shell depth 9 | nvdb | 2287433 (2001504) | 26 | 22.0 | 32.8 | 27.2 | 11.3 | 10.5 | 672 | 561.1 | 21.0 (17.6) | 10.4 |
| 3D shell depth 9 | nvdb-compact | 2287433 (2001504) | 12 | 40.6 | 32.3 | 25.0 | 12.3 | 8.0 | 796 | 603.6 | 24.9 (18.9) | 34.6 |

## Reading it

- **Memory.** On every non-uniform tree: nano-vdb compact 12-14 B per node,
  voxel 14-25, nano-vdb dense 26-44, memindex 234-615. On the uniform trees
  every implementation is at 8-9 B, but memindex at 77 B in 3D.
- **Build and free.** On the non-uniform trees, the trees build 6-28 times
  faster than memindex; compact nodes cost 1.4-2.9 times dense ones (except
  the 2D river, where they were faster). nano-vdb dense frees fastest in
  every case.
- **Morton traversal (dfs).** With cursors, the trees take 6.7-27 ns per
  node, voxel usually the fastest, against 18-54 ns for memindex: 1.7-5.5
  times faster.
- **Over a level list.** With cursors, 4-11 ns per node for every tree,
  against 7-21 for memindex.
- **Neighbour box.** With cursors, every tree beats memindex. In 2D by
  30-55%, the three trees alike. In 3D by 11-19% for nano-vdb's dense
  accessor, 4-9% for its compact one, and 2-13% for voxel's bottom-only
  cursor. The accessor helps where a box spills into the next leaf: it
  restarts from the cached lower node, not from the root.
- In these runs memindex's neighbour box was slower, relative to the trees,
  than in the earlier runs of 2026-11-09 on a quiet machine, where it led
  voxel by up to 36%. Rerun on a quiet machine before weighing differences
  of a few tens of percent.
