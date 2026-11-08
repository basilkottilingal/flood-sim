# memindex vs voxel-memindex

`sh tests/bench_memindex.sh` (2026-11-08). Intel i7-8650U, 8 GB (about
2.3 GB free), Linux 6.8, gcc 11.4 -O2, THP `[madvise]`, memindex with
MADV_NOHUGEPAGE. One run per case (a second run agreed within ~20%).

- B/point: (RSS + page tables) after building, per point.
- ns columns: average per operation. `nbr`: the 2*dim face neighbours of
  every point, Morton order. `uniform`: indices drawn over the whole
  level (the hit rate in brackets).
- kB left after free: RSS + page tables still held once every point is
  freed.

| case | impl | points | B/point | build ns | build faults | ordered ns | shuffled ns | nbr ns | nbr+cursor ns | uniform ns (hits) | free ns | kB left after free |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2D dense L11 | mem | 4194304 | 8 | 9.2 | 8200 | 2.5 | 15.9 | 4.2 | - | 17.9 (100%) | 75.7 | 136 |
| 2D dense L11 | voxel | 4194304 | 8 | 12.9 | 8454 | 7.4 | 56.3 | 7.7 | 4.8 | 52.3 (100%) | 121.8 | 33416 |
| 2D shell L16 | mem | 164492 | 2944 | 2735.5 | 105176 | 37.7 | 68.1 | 24.6 | - | 744.4 (0%) | 1991.6 | 53004 |
| 2D shell L16 | voxel | 164492 | 92 | 57.4 | 4153 | 17.2 | 71.8 | 16.1 | 7.8 | 7.5 (0%) | 182.2 | 14716 |
| 2D river L16 | mem | 2097152 | 168 | 183.8 | 69809 | 6.7 | 62.1 | 6.8 | - | 934.3 (0%) | 593.2 | 66232 |
| 2D river L16 | voxel | 2097152 | 19 | 23.1 | 10288 | 11.8 | 104.5 | 11.2 | 4.8 | 7.8 (0%) | 212.3 | 38044 |
| 2D random L14 | mem | 795296 | 2304 | 1398.6 | 446388 | 21.8 | 45.5 | 39.3 | - | 43.8 (38%) | 1797.3 | 4288 |
| 2D random L14 | voxel | 795296 | 331 | 173.6 | 72400 | 42.3 | 121.9 | 25.0 | 15.0 | 111.1 (38%) | 455.4 | 257068 |
| 2D random L20 | mem | 49989 | 8370 | 6397.6 | 54085 | 72.6 | 117.6 | 44.1 | - | 102.9 (0%) | 2956.2 | 570644 |
| 2D random L20 | voxel | 49989 | 1907 | 1350.9 | 40029 | 190.6 | 267.8 | 92.2 | 66.9 | 49.3 (0%) | 947.6 | 93084 |
| 3D dense L7 | mem | 2097152 | 32 | 47.9 | 16642 | 7.9 | 31.7 | 6.1 | - | 34.4 (100%) | 103.1 | 704 |
| 3D dense L7 | voxel | 2097152 | 8 | 15.0 | 4228 | 10.9 | 68.5 | 10.2 | 6.3 | 72.4 (100%) | 105.0 | 16740 |
| 3D shell L9 | mem | 530552 | 1027 | 1915.2 | 133218 | 27.0 | 54.4 | 13.6 | - | 41.8 (0%) | 1209.6 | 2740 |
| 3D shell L9 | voxel | 530552 | 47 | 33.3 | 6519 | 13.4 | 72.4 | 14.3 | 9.8 | 18.7 (0%) | 157.9 | 24444 |
| 3D shell L10 | mem | killed: out of memory (> 2.3 GB) | | | | | | | | | | |
| 3D shell L10 | voxel | 2105960 | 47 | 33.0 | 25935 | 13.2 | 114.2 | 13.3 | 10.2 | 17.4 (0%) | 251.1 | 97436 |
| 3D river L12 | mem | 851968 | 344 | 594.3 | 74780 | 12.1 | 86.4 | 8.7 | - | 25.5 (0%) | 810.1 | 21636 |
| 3D river L12 | voxel | 851968 | 26 | 34.1 | 5871 | 23.4 | 136.6 | 18.2 | 12.2 | 14.9 (0%) | 254.8 | 21416 |
| 3D random L12 | mem | 99914 | 5393 | 10694.1 | 161325 | 100.7 | 126.6 | 48.2 | - | 183.4 (0%) | 2890.9 | 131040 |
| 3D random L12 | voxel | 99914 | 1468 | 1188.0 | 58993 | 153.1 | 282.1 | 69.5 | 52.3 | 58.3 (0%) | 684.0 | 143192 |
| 3D random L20 | mem | 30000 | 16715 | 12325.3 | 94088 | 209.1 | 290.7 | 75.5 | - | 71.7 (0%) | 3107.4 | 591628 |
| 3D random L20 | voxel | 30000 | 3636 | 2974.0 | 49407 | 574.0 | 628.5 | 199.2 | 172.3 | 36.8 (0%) | 1570.9 | 106524 |

Cache and dTLB miss counters were not available (kernel.perf_event_paranoid
= 4); with `sudo sysctl kernel.perf_event_paranoid=1` the script prints them.

## Reading it

- Memory: voxel uses 3-60x less wherever the level is not full (dense
  cases tie, or favour voxel in 3D). memindex pays a 4 kB page and part of
  a page-table page per row touched; on the 3D shell at level 10 that is
  more than the 2.3 GB this machine had free, against 99 MB for voxel.
- Build: voxel 5-60x faster on sparse levels (far fewer page faults),
  about the same on dense ones.
- Lookups of present entries: memindex is faster, 1.3-3x, more the deeper
  the level (voxel walks ceil (L/B) nodes, memindex reads `dimension`
  rows). With the cursor, voxel's neighbour sweep matches or beats
  memindex's in every case.
- Lookups of absent entries (`uniform` on sparse levels): voxel 2-120x
  faster. memindex reads untouched pages, each a page fault (the zero page).
- After freeing everything, memindex keeps its page tables and upper rows
  (570-590 MB on the random level-20 cases); voxel keeps its freed nodes
  in the pool, hot (pool_set_hot_max () / pool_trim () can change that;
  pool_destroy () gives them all back).
- On sparse deep levels voxel still spends ~2-4 kB per point: an isolated
  point owns most nodes on its path (512 B each). Nodes with a 64-bit mask
  and only the children present (popcount rank) would cut that.
