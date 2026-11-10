# memindex vs voxel-memindex vs nano-vdb: isolated lookups

`sh tests/bench_memindex.sh`. Intel i7-8650U laptop, 8 GB, Linux 6.8, gcc 11.4 `-O2 -mpopcnt`, THP
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

nano-vdb's dense nodes skip the scattered 3D cases (`random` L12, L20):
about 300 KB per isolated point. memindex finished the 3D shell at level 10
this time (2 GB); with less memory free it is killed.

- B/point: (RSS + page tables) after building, per point.
- ns columns: average per operation. `nbr`: the 2*dim face neighbours of
  every point, Morton order. `uniform`: indices drawn over the whole level
  (the hit rate in brackets).
- kB left after free: RSS + page tables still held once every point is
  freed (pools keep freed nodes hot).

| case | impl | points | B/point | build ns | build faults | ordered ns | shuffled ns | nbr ns | nbr+cursor ns | uniform ns (hits) | free ns | kB left after free |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2D dense L11 | mem | 4194304 | 8 | 10.2 | 8200 | 5.5 | 22.8 | 9.6 | - | 21.4 (100%) | 99.5 | 140 |
| 2D dense L11 | voxel | 4194304 | 8 | 20.7 | 8454 | 22.1 | 61.1 | 13.7 | 4.8 | 50.0 (100%) | 124.2 | 33420 |
| 2D dense L11 | nvdb | 4194304 | 9 | 14.1 | 9482 | 5.9 | 68.3 | 14.1 | 6.0 | 66.3 (100%) | 119.5 | 37552 |
| 2D dense L11 | nvdb-compact | 4194304 | 9 | 27.5 | 9647 | 12.9 | 81.4 | 19.3 | 5.1 | 78.4 (100%) | 151.7 | 3100 |
| 2D shell L16 | mem | 164492 | 2944 | 2245.2 | 105176 | 41.3 | 84.5 | 62.2 | - | 837.5 (0%) | 2410.5 | 53012 |
| 2D shell L16 | voxel | 164492 | 92 | 95.0 | 4153 | 62.6 | 113.1 | 28.1 | 12.7 | 12.1 (0%) | 237.0 | 14720 |
| 2D shell L16 | nvdb | 164492 | 113 | 72.1 | 5200 | 24.1 | 132.6 | 30.2 | 9.9 | 19.2 (0%) | 116.4 | 18204 |
| 2D shell L16 | nvdb-compact | 164492 | 22 | 43.4 | 1075 | 14.9 | 52.5 | 25.7 | 7.2 | 15.1 (0%) | 81.5 | 3452 |
| 2D river L16 | mem | 2097152 | 168 | 193.5 | 69809 | 9.1 | 64.0 | 14.3 | - | 957.3 (0%) | 963.0 | 66236 |
| 2D river L16 | voxel | 2097152 | 19 | 27.4 | 10288 | 15.7 | 145.4 | 21.1 | 6.7 | 16.4 (0%) | 274.2 | 38048 |
| 2D river L16 | nvdb | 2097152 | 22 | 22.7 | 12207 | 8.5 | 104.6 | 22.2 | 7.3 | 24.7 (0%) | 218.8 | 45308 |
| 2D river L16 | nvdb-compact | 2097152 | 12 | 34.4 | 6888 | 11.7 | 119.4 | 24.4 | 5.5 | 17.7 (0%) | 200.2 | 24380 |
| 2D random L14 | mem | 795296 | 2304 | 1500.7 | 446388 | 47.7 | 50.6 | 107.6 | - | 47.5 (38%) | 1867.0 | 4296 |
| 2D random L14 | voxel | 795296 | 331 | 166.2 | 72400 | 52.6 | 144.3 | 35.1 | 13.3 | 134.1 (38%) | 590.3 | 257072 |
| 2D random L14 | nvdb | 795296 | 368 | 181.6 | 79398 | 42.2 | 181.2 | 34.0 | 10.9 | 163.9 (38%) | 326.9 | 285768 |
| 2D random L14 | nvdb-compact | 795296 | 52 | 81.3 | 13847 | 15.2 | 186.5 | 23.1 | 6.4 | 154.7 (38%) | 369.6 | 40552 |
| 2D random L20 | mem | 49989 | 8370 | 5753.4 | 54085 | 66.0 | 111.1 | 99.3 | - | 102.0 (0%) | 2797.5 | 570648 |
| 2D random L20 | voxel | 49989 | 1907 | 1280.6 | 40031 | 192.6 | 271.4 | 121.0 | 75.5 | 50.8 (0%) | 1024.2 | 93088 |
| 2D random L20 | nvdb | 49989 | 7426 | 4836.3 | 136701 | 248.3 | 364.6 | 161.7 | 77.4 | 54.4 (0%) | 623.0 | 362500 |
| 2D random L20 | nvdb-compact | 49989 | 462 | 465.2 | 9873 | 34.8 | 255.2 | 38.5 | 16.6 | 32.0 (0%) | 644.9 | 22548 |
| 3D dense L7 | mem | 2097152 | 32 | 45.9 | 16642 | 6.8 | 34.4 | 17.1 | - | 37.7 (100%) | 108.8 | 704 |
| 3D dense L7 | voxel | 2097152 | 8 | 14.7 | 4228 | 8.4 | 63.5 | 13.4 | 6.6 | 64.8 (100%) | 98.7 | 16748 |
| 3D dense L7 | nvdb | 2097152 | 8 | 13.4 | 4248 | 6.4 | 48.1 | 12.7 | 5.5 | 47.8 (100%) | 41.8 | 17056 |
| 3D dense L7 | nvdb-compact | 2097152 | 8 | 33.8 | 4389 | 10.5 | 65.2 | 17.2 | 7.1 | 63.4 (100%) | 207.4 | 4816 |
| 3D shell L9 | mem | 530552 | 1027 | 1820.8 | 133218 | 21.2 | 43.7 | 26.1 | - | 37.8 (0%) | 1086.2 | 2748 |
| 3D shell L9 | voxel | 530552 | 47 | 41.5 | 6519 | 15.8 | 72.9 | 19.8 | 12.0 | 34.6 (0%) | 210.3 | 24444 |
| 3D shell L9 | nvdb | 530552 | 98 | 52.9 | 12977 | 16.3 | 67.7 | 18.6 | 9.3 | 15.8 (0%) | 46.3 | 51016 |
| 3D shell L9 | nvdb-compact | 530552 | 16 | 43.6 | 2501 | 9.7 | 37.4 | 16.7 | 10.0 | 9.2 (0%) | 71.9 | 8540 |
| 3D shell L10 | mem | 2105960 | 2062 | 2761.8 | 1059956 | 98.4 | 110.5 | 27.0 | - | 65.7 (0%) | 1577.3 | 14908 |
| 3D shell L10 | voxel | 2105960 | 47 | 33.8 | 25935 | 11.8 | 104.0 | 16.4 | 9.6 | 15.4 (0%) | 250.0 | 97488 |
| 3D shell L10 | nvdb | 2105960 | 98 | 50.3 | 51420 | 14.5 | 96.3 | 16.7 | 9.2 | 18.7 (0%) | 127.7 | 201880 |
| 3D shell L10 | nvdb-compact | 2105960 | 16 | 38.5 | 9766 | 8.2 | 103.7 | 16.8 | 8.8 | 10.6 (0%) | 165.6 | 33100 |
| 3D river L12 | mem | 851968 | 344 | 571.2 | 74780 | 9.8 | 60.1 | 14.9 | - | 22.6 (0%) | 630.5 | 21688 |
| 3D river L12 | voxel | 851968 | 26 | 25.6 | 5871 | 12.1 | 72.4 | 17.6 | 9.2 | 8.2 (0%) | 142.4 | 21464 |
| 3D river L12 | nvdb | 851968 | 57 | 35.0 | 12310 | 9.2 | 53.4 | 14.7 | 6.7 | 4.7 (0%) | 44.4 | 47192 |
| 3D river L12 | nvdb-compact | 851968 | 14 | 40.9 | 3219 | 9.1 | 55.6 | 17.2 | 8.8 | 4.5 (0%) | 97.7 | 11448 |
| 3D random L12 | mem | 99914 | 5393 | 8153.8 | 161325 | 72.7 | 116.6 | 59.8 | - | 159.4 (0%) | 2620.8 | 131088 |
| 3D random L12 | voxel | 99914 | 1468 | 989.2 | 58993 | 173.8 | 253.3 | 77.2 | 45.6 | 72.9 (0%) | 699.2 | 143200 |
| 3D random L12 | nvdb-compact | 99914 | 456 | 655.5 | 19342 | 39.8 | 239.4 | 29.8 | 14.3 | 47.6 (0%) | 958.2 | 44528 |
| 3D random L20 | mem | 30000 | 16716 | 14702.9 | 94088 | 169.1 | 232.0 | 87.9 | - | 66.9 (0%) | 2960.1 | 591632 |
| 3D random L20 | voxel | 30000 | 3636 | 2550.0 | 49403 | 475.1 | 568.6 | 176.1 | 152.2 | 46.4 (0%) | 1600.6 | 106528 |
| 3D random L20 | nvdb-compact | 30000 | 7367 | 6963.8 | 97902 | 279.0 | 456.2 | 96.7 | 69.6 | 29.8 (0%) | 1113.7 | 215832 |

## Reading it

- **Memory.** nano-vdb compact is the smallest on every level that is not
  full: 12-22 B per point on the shells and rivers (voxel 19-92, nano-vdb
  dense 22-113, memindex 168-2944). On scattered points it takes 52-462 B in
  2D, against 331-1907 (voxel) and 2304-8370 (memindex), and 456 B on the 3D
  points at level 12 (voxel 1468, memindex 5393). Its weak case is scattered
  3D points at level 20: 7.4 KB a point, because the masks of a 3D upper
  node are 4 KiB of bits whatever it holds; voxel needs 3.6 KB, memindex
  16.7 KB. Dense nano-vdb is worse still there, hence skipped; in 2D, whose
  upper nodes have 32^2 slots, it ran, at 7.4 KB a point at level 20.
- **Build.** On the shells and rivers the three trees build 6-52 times
  faster than memindex (far fewer page faults); on scattered points 1.2-12
  times. Compact nodes cost 2-2.5 times dense ones on full levels
  (inserting into packed arrays).
- **Present entries, random order.** memindex is fastest on the full
  levels, the 2D river and the scattered points. Compact nano-vdb beats it
  on the 2D shell (53 against 85 ns) and the 3D shell at level 9 (37
  against 44); dense and compact nano-vdb both beat it on the 3D river and
  the 3D shell at level 10.
- **Neighbours with a cursor.** One of the nano-vdb variants is the fastest
  in every case but the 2D full level (voxel), and every tree's cursor
  beats memindex's plain lookups except voxel on the scattered 3D points at
  level 20.
- **Absent entries** (0% hits): 4.5-73 ns for the trees, which stop at the
  first empty slot, against 23-957 ns for memindex, which reads untouched
  pages.
- **After freeing.** Among the trees, compact nano-vdb keeps the least (its
  nodes are small); memindex keeps its page tables (570-590 MB on the
  level-20 cases).
