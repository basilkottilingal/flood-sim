# libhmesh

A simplicial manifold of dimension `m` immersed in `R^D`, `0 <= m <= D <= 3`
(`D` cwrap's `-dimension=`): points, a curve (a river network), a surface, a
volume. A half-edge mesh, generalised to every level of the complex, with
fields on every level stored as a structure of arrays.

```c
manifold (2);                          /* hmesh: a surface in R^D */
uint32_t a = hm_vertex_add (hmesh, xa), b = ..., c = ...;
uint32_t t = hm_add (hmesh, 2, (uint32_t []) { a, b, c });
int depth = hm_scalar_new (hmesh, 2, "h");     /* on triangles */
int flux  = hm_scalar_new (hmesh, 1, "q");     /* on edges */
hm_foreach (hmesh, 2, s)
  hm_val (hmesh, depth, 0)[s] = ...;
```

## The model

**Levels.** Level `k`, `0 <= k <= m`, holds the k-simplices. Every face of
every simplex is stored once, at its own level, with its own id: a
triangle's three edges are level 1 simplices, so that each level carries
its own scalars and vectors (depths on triangles, fluxes on edges, heights
on points).

**Halves.** A k-simplex `s = (v0 ... vk)`, `k >= 1`, has `k + 1` halves,
`h = (k + 1) s + i`, one for each facet: the facet opposite `vi`. Each half
has

| | |
|---|---|
| `hm_half_vertex ()` | `vi`, the vertex opposite its facet |
| `hm_half_facet ()` | the facet's id, at level `k - 1` |
| `hm_half_sign ()` | the incidence number: +1 if the orientation `s` induces on the facet, `(-1)^i (v0 ... vi ... vk without vi)`, is the facet's own |
| `hm_half_next ()` | `i + 1 mod k + 1`, computed, not stored |
| `hm_half_twin ()` | the next half, of another k-simplex, on the same facet |

So every simplex but a point has a twin and a next:

| k | halves | twin cycles about |
|---|---|---|
| 1 | 2 half-points (edge ends) | a point: the edges at it |
| 2 | 3 half-edges | an edge: the triangles on it |
| 3 | 4 half-faces | a triangle: the tetrahedra on it |

**Twins are cycles.** At level `m`, inside the manifold, a facet has two
halves and twin is the half-edge mesh's twin: `twin (twin (h)) == h`, the
neighbour across, of the opposite sign when the orientation is consistent.
On a boundary the cycle is `h` alone (`hm_half_boundary ()`). Where more than
two meet, as at a confluence of three river branches, the cycle holds them
all, so a junction needs no special case. Below level `m` a cycle holds
every coface, in no particular order: the edges at a point of a surface.

For a triangle the halves are the half-edges of a half-edge mesh: half `i`
goes from `v(i+1)` to `v(i+2)` (`hm_half_origin ()`), and
`twin (next (next (h)))` turns about the origin of `h`.

**Uses.** A (k-1)-simplex keeps one half on it at level `k`, `hm_use ()`,
the entry into its cycle: `hm_use (M, 0, v)` starts the walk through the
edges at the point `v`.

**Signs.** The signs are the boundary operator: `d (v0, v1) = v1 - v0`, and
`d d = 0` at every level (tested). A flux on an edge, added to a triangle
with `hm_half_sign ()`, is a conservative divergence on any orientation of
the edges.

## Storage: fast push and remove anywhere

Per level, a structure of arrays: vertices, twins and facets at
`[(k + 1) s + i]` (the facet's sign in its top bit), one use per simplex
below the top level, the points' coordinates `x[0], x[1], x[2]`, and each
component of each field an array of its own. Ids are 32 bit.

- **Push.** An id from the stack of freed ids, the last freed first (its
  memory is still in cache), else the end of the arrays, which double as
  they fill. A new simplex's field values are zero.
- **Remove, anywhere.** Its bit is cleared and its id pushed on the stack,
  in O(1). Its halves leave their facets' cycles, which costs the length of
  each cycle.
- **Ids are stable** until `hm_compact ()`: twins, facets, uses and any id
  held outside the mesh stay valid through any sequence of additions and
  removals (a flip, a split, a collapse).
- **Iteration.** `hm_foreach (M, k, s)` visits the ids in use, in order,
  skipping holes 64 at a time (a bit per id). `hm_dense ()` says there are
  no holes, so that a plain loop over `0 ... hm_count () - 1` will do (for
  SIMD).
- **Compaction.** `hm_compact (M, remap)` moves every level's simplices down
  to `0 ... n - 1` in the order of their ids, and every reference follows.
  `remap[k]`, if asked for, maps old ids to new, for links held outside the
  mesh (to blocks of the tree, to ranks).

`hm_find (M, k, v)` finds a simplex from its vertices, in any order, with no
hash table: the k-simplices on the facet `v[1] ... v[k]` are its use's
cycle, and the one sought is the one whose half there is opposite `v[0]`.
`hm_add ()` finds or makes each face this way, so a mesh is built from its
top simplices alone.

## API

| | |
|---|---|
| `manifold (dim)`, `manifold_free ()` | the program's manifold, `hmesh`; -1 unless `0 <= dim <= D` |
| `hmesh_new (dim)`, `hmesh_free (M)` | the same, for a mesh of one's own |
| `hm_vertex_add (M, x)` | a point at `x` (D coordinates) |
| `hm_add (M, k, v)` | the k-simplex `v[0] ... v[k]`, its faces found or made, each new face in the orientation its first coface induces; the one already there if any, with its own orientation; `HM_NONE` for a vertex repeated or not in use, or `k > m` |
| `hm_find (M, k, v)` | the k-simplex of vertices `v`, any order, or `HM_NONE` |
| `hm_remove (M, k, s)` | `s` alone; -1 while a (k+1)-simplex has it as a face |
| `hm_remove_closure (M, k, s)` | `s`, and every face of it that nothing else uses, points included |
| `hm_compact (M, remap)` | no holes; `remap` may be `NULL` |
| `hm_field_new (M, k, name, ncomp, esize)` | a field on level `k`: `ncomp` arrays of `esize` byte elements |
| `hm_scalar_new`, `hm_vector_new` | `HmReal` fields of 1 and D components (`HmReal` is double, or float with `-DBAMR_SINGLE`, as libb-amr's `Real`) |
| `hm_field (M, f, c)`, `hm_val (M, f, c)` | component `c`'s array, `void *` or `HmReal *`; it moves when its level grows |
| `hm_x (M, d)` | the points' coordinate `d` |
| `hm_measure`, `hm_volume_signed`, `hm_centroid` | length, area or volume (a Gram determinant, any `k <= D`), signed when `k = m = D`, centroid |
| `hm_check (M)` | every invariant, O(size), for tests |
| `hm_set_domain (M, origin, size)` | the box `[origin, origin + size]^D` (closed) the points must stay in; -1 if a point is outside already |
| `hm_vertex_move (M, v, x)`, `hm_touch (M)` | a checked move; after writing `hm_x ()` directly |
| `hm_compact_order (M, order, remap)` | the compaction, each level in a given order |
| `hm_clip (M, origin, size)`, `hm_tree_clip (M)` | the mesh cut to the box, which then holds |
| `hm_domain_faces (M, k, s)`, `hm_half_on_domain (M, k, h)` | the faces of the box a simplex, or a boundary half's facet, lies on |

The library and the program must agree on `D`: a program built with
another dimension does not link (`hmesh_dimension_is_D`).

## Crossing the box (hmesh-clip.c)

A manifold built larger than the domain is cut to it:
`hm_clip (M, origin, size)`, or `hm_tree_clip (M)` for the tree's domain.
Build the mesh with its points anywhere, clip, and the box then holds
(`hm_set_domain ()`). For points that moved out, write them through `hm_x ()`
and clip again.

- **One plane at a time** (Sutherland and Hodgman, on the complex). A point
  within `HM_CLIP_SNAP` (1e-10) of the box's size from the plane is moved
  onto it, so no sliver is made. The others are in or out.
- **One new point an edge.** Each edge that crosses gets one point where it
  crosses. The edge is stored once, so every simplex on it shares the point
  and the cut mesh stays conforming.
- **The part inside of each simplex across the plane**, maximal ones only
  (the top level, and lower ones that are the face of nothing):
  - a segment: its end out moved to the cut;
  - a triangle or quadrangle: the fan from its lowest point;
  - a polyhedron: the cone from its lowest point to the fans of the faces it
    is not on.

  Both sides of a shared face choose by the same global ids, so they agree
  (Dompierre et al.). No point is added inside.
- **What lies out goes**, with the faces nothing uses then.
- **Orientation kept.** Each new simplex has the sign of the determinant of
  its points' barycentric coordinates in the simplex it comes from.
- **Fields.** A new point's HmReal values are interpolated along its edge
  (exact for linear fields). A new maximal simplex takes the values of the
  one it comes from. Other new faces are zero.

**Twin and next on the box.** A half whose facet is on the box is a
boundary half: `twin (h) == h`, and `next` is as ever. `hm_half_on_domain
(M, k, h)` tells it from the manifold's own boundary by the faces of the box
its facet lies on, as a mask: bit `2 d` the lower face along `d`, `2 d + 1`
the upper. `hm_domain_faces (M, k, s)` gives the same for any simplex. The
cut points lie on the faces exactly (the coordinate is the plane's), so
these tests are exact.

Tested (`tests/test_clip.c`), each against what the cut must give exactly:

| case | result |
|---|---|
| a straight curve | its length inside, its two ends on x = 0 and x = 1 |
| a circle crossing every side of the square | four arcs, eight ends on the box, Euler 4 |
| a rotated triangulation covering the square | area 1, boundary all on the box, of length 4 |
| one half out | area 1/4, boundary 1 on the box and 1 off |
| a tilted plane in the cube, cut by its sides and by z = 0, 1 | area 2/3 sqrt (13/4) |
| a rotated cube of tetrahedra covering the cube | volume 1, boundary 6, on all six faces |

Every case also checks the invariants, orientations, Euler characteristic,
`d d = 0`, the twins' signs, and the fields. The suite also covers snapping,
isolated points and a dangling edge, and clipping twice. Two bugs made on
purpose are caught: fans not from the lowest point (the cube's cut faces
stop matching) and no orientation fix.

## In the tree (hmesh-tree.h)

The manifold inside libb-amr's tree, after `hm_tree_clip ()` if it
crossed the domain's boundary.

- **The box.** `hm_tree_bind (M)` takes the tree's domain from
  `tree_get_origin ()` and `tree_get_size ()` (libb-amr unchanged), and from
  then on a point outside it is refused, by `hm_vertex_add ()` and
  `hm_vertex_move ()`; `hm_check ()` finds one written outside through
  `hm_x ()`. If the domain moves (`tree_shift_resize ()`), the next cache
  binds again, and fails if that leaves a point outside.
- **The keys.** 64 bit Morton keys, x in the lowest bit as in
  `index_child ()`, so that sorted keys run along `tree_caches ()->leaves`.
  52 bits an axis in 1D, 32 in 2D, 21 in 3D (`hmesh-morton.h`). The tree
  must not be deeper than the keys can tell its cells apart:
  `tree_get_depth_max () + log2 (block_m) <= bits`, so 17 with blocks of 16
  in 3D, checked by `hm_tree_bind ()` and the cache.
- **The cache.** `hm_tree_cache (M)` gives the points sorted along the curve
  (`order[]`, `key[]`), each point's leaf (`leaf_of[]`, an index into
  `tree_caches ()->leaves`, so its block and its `pid`), each leaf's range of
  points (`start[]`, `stop[]`, `hm_foreach_leaf_point ()`), the points in no
  leaf, and the duplicates (two points closer than `size / 2^bits` on every
  axis). One walk along the leaves and the sorted points together; no point
  located by descending the tree. It is made again only when
  `tree_generation` or the mesh's generation (points added, removed, moved,
  renumbered) changed.
- **The sort.** To the tree's finest cells only (`D (depth_max + log2 m)`
  bits, 11 a radix pass), as that is all the leaves and cells need. It
  starts from the last order: nothing to do if no point changed its cell,
  insertion if a few did, else radix.
- **The imaginary tree** with at most one point a leaf is the sorted keys:
  the smallest cell about a point holding no other is set by the longest
  common prefix of its key with its neighbours'. Nothing more is stored.
- **Memory order.** `hm_compact_morton (M, remap)` renumbers every level
  along the curve (points by their whole keys, a k-simplex by its
  centroid's), so that a leaf's points are the ids `start[q] ... stop[q] - 1`
  and their fields are read in memory order.

The cache's cost, 2^20 points in about as many cells (a uniform tree,
blocks of 16), -O2, a loaded laptop, ns a point:

| | 1D | 2D | 3D |
|---|---|---|---|
| from nothing | 52 | 64 | 70 |
| after every point moved 1e-5, ids in no order | 46 | 59 | 52 |
| after `hm_compact_morton ()`, nothing moved | 27 | 42 | 47 |
| after `hm_compact_morton ()`, every point moved | 43 | 54 | 67 |

For comparison, a step of libvof's advection is 30-60 ns a cell (2D), and
a curve or surface in the domain has far fewer points than the grid has
cells: rebuilding the cache at every step is cheap. Most of its time is
memory traffic, not the sort itself: computing the keys (about 10 ns a
point), the radix passes (about 4 ns a point each, 2 to 6 of them), and
the walk.

## Tests

```sh
sh tests/run.sh            # 1-3D, double and float, ASan and UBSan, then -O2
QUICK=1 sh tests/run.sh    # without the -O2 runs
```

- **Pool.** Ids reused, the last freed first. Holes skipped in order.
  Compaction's map, with the coordinates moved along.
- **Meshes.** A curve with a junction (three edges at a point, the sum of
  the signs at each point), a triangulated square, flat in 2D and tilted in
  3D, an octahedron, and a cube of 6 n^3 Kuhn tetrahedra. For each: counts,
  Euler characteristic, the twin cycles (two of opposite signs inside, one on
  the boundary: 4n edges, 12 n^2 faces), `d d = 0`, length, area and
  volume, and the one-ring of a point by `twin (next (next (h)))`.
- **Churn.** Half the top simplices out at random, with their closure or
  with their unused faces down to the edges, then back in. Invariants after
  each stage, counts and Euler characteristic as before, fields on every
  level still the measures of their simplices, then the same after
  compaction.
- **Refusals.** And `hm_check ()` finding eight kinds of corruption.
- **In the tree** (`test_tree.c`, linked with libb-amr). A tree refined
  over half its domain, the domain moved and resized: points outside, NaN,
  on the upper faces and the corners; a tree too deep for the keys (3D);
  each point in the one leaf whose box holds it, by brute force (a face
  shared by two leaves to the upper one), the ranges, the keys, the
  duplicate found; again after the points moved, were written through
  `hm_x ()`, were removed and added, after the tree was refined, after the
  domain moved from under a point and back; `hm_compact_morton ()`. A
  cache made wrong on purpose (one leaf off) is caught.

Measured at -O2 on a loaded laptop (taking out a random simplex and adding
it back; the time is dominated by memory misses at random ids, not by the
pool):

| | per out and in |
|---|---|
| a point, 2^20 of them, with a scalar and a vector | 70-170 ns |
| an edge of a curve of 2^17 | 90 ns |
| a triangle of 131072 (256^2 squares) | 0.5-0.9 us |
| a tetrahedron of 82944 (24^3 cubes) | 1.2-1.4 us |

A triangle or a tetrahedron costs more because `hm_add ()` finds each face
through the cycles at its vertices, and each step of a cycle lands somewhere
else in memory.

## Next

- **Periodic boxes.** The boundary halves on a face linked to those on the
  opposite one. Their facets are distinct simplices, a translation apart,
  so the twin cycle, which joins the halves of one facet, would need an
  identification of the two.
- **Clipping as points move.** Today the whole mesh is clipped again. Only
  the simplices near the box need to be.
- **Ranks.** Each point's rank under MPI, from its leaf's `pid`.
- **Order.** The halves' vertices, twins and facets could be interleaved
  into one record, if cycle walks are where the time goes.
- **Topological operations.** Edge flip, split and collapse, built on
  `hm_add ()` and `hm_remove ()` (ids stay stable through them).
- **Ordered cycles.** The cofaces about an edge of a volume, or about a
  point of a surface, sorted by angle.
- **cwrap grammar.** Scalars and vectors on levels, and loops over
  simplices and halves.
