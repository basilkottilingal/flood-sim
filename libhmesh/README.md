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

The library and the program must agree on `D`: a program built with
another dimension does not link (`hmesh_dimension_is_D`).

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

- **Coupling.** Each point linked to its leaf block of libb-amr, and its
  rank under MPI. The map from `hm_compact ()` is there for that.
- **Order.** Compaction in Morton order of the centroids, so that the
  mesh's memory follows the tree's. The halves' vertices, twins and facets
  could also be interleaved into one record, if cycle walks are where the
  time goes.
- **Topological operations.** Edge flip, split and collapse, built on
  `hm_add ()` and `hm_remove ()` (ids stay stable through them).
- **Ordered cycles.** The cofaces about an edge of a volume, or about a
  point of a surface, sorted by angle.
- **cwrap grammar.** Scalars and vectors on levels, and loops over
  simplices and halves.
