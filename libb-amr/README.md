# block based amr

## Ghost cells

Every block stores its m^d interior cells and `GHOSTS` (default 2) ghost
layers around them, as PARAMESH/FLASH, Athena++ and Parthenon do. A stencil
reads its neighbours in place:

```c
ghosts_update (u);                     /* block-ghosts.h */
BlockShape shape = tree_block_shape ();
foreach_leaf_parallel () {
  const Real * a = block_data (sblocks, u);
  Real * b = block_data (sblocks, v);
  block_cells_simd (shape, 0)
    b [c] = cell_at (a, 1, 0) + cell_at (a, -1, 0) + cell_at (a, 0, 1)
          + cell_at (a, 0, -1) - 4 * cell_at (a, 0, 0);
}
```

`ghosts_update (s, ...)`:
1. **restricts**, finest level first, every parent from its children
   that a leaf may read: those of the coarsest leaf's level and finer, as
   a parent is read only by a leaf of its own level beside it (a uniform
   tree: none). `ghosts_update_levels ()` and `tree_restrict ()` restrict
   every parent. The loops run over the tree's list of parents
   (`tree_caches ()->parents_by_level`), not every node;
2. **fills the ghosts**, coarsest level first, each block's from
   - the node of its own level beside it (a leaf, or a restricted parent):
     rows copied;
   - otherwise the leaf one level coarser (2:1 balance), by the scalar's
     prolongation. That reads the coarse leaf's interior and first ghost
     layer, which are already valid because that level was done first;
   - outside the domain, the scalar's boundary condition, one axis after
     the other, so that corners come out right.

Within a level the blocks are independent and run on all threads.
`ghosts_update_levels ()` fills parents too, for coarse-level stencils.

**Per-scalar operations**, anywhere, at any time: a prolongation and a
restriction (`ScalarOps`, `scalar_set_ops ()`), and a boundary condition per
face of the domain, each an integer id. A boundary condition is bound with
cwrap's `ghost_ops ()`, face by face, `boundary.left (op)`, `.right`,
`.bottom`, `.top`, `.back`, `.front`, or `boundary.all (op)` for every face:
`ghost_ops (u.x, boundary.all (bc_neumann), boundary.left (bc_odd))`. A
plain `boundary (op)`, or a face the dimension lacks, is an error;
`boundary.all` after a face, in one call, a warning, as it binds that face
again. `scalar_get_boundary (s, face)` reads one (`GHOST_LEFT`, ...). 0 takes the default: `prolong_linear`, `restriction_mean`,
`bc_neumann`. Also given: `prolong_minmod` (no new extrema, for a water
depth), `prolong_injection` and `bc_odd`. One prolongation fills both the
children of a refined block and the ghosts at a level jump, so the two
agree. The operators are cwrap's `@prolongation`, `@restriction` and
`@boundary` definitions (`block-ops.h`). Prolongations and restrictions are
written with stencils over two levels, `foreach_prolongation (pr)` with
`coarse (s, i, j)` and `child.x`, `foreach_restriction (p)` with
`foreach_child ()` (`block-ghosts.h`); boundary conditions with
`boundary_cells ()`. They may read other scalars of the same nodes, for
coupled operators such as `h` with `zb`. `ghosts_update ()` calls them
through a `switch` on the id, each case inlined, not through a function
pointer. The file that dispatches must see every definition:
`block-ghosts.c` for the given ones; a program with operators of its own
builds its own in its place (`block-ghosts.h`), as `tests/test_ghosts.c`
does.

**Coupled operators**, which read other scalars, as the water depth `h`
prolonged from the surface `h + zb`: the operator names them in a second
parameter list, `@prolongation prolong_lake (const Prolongation * pr) (scalar
zb)`, and `ghost_ops (h, prolongation (prolong_lake, zb))` binds them.
`ghost_ops_set ()` refuses a wrong count or a cycle. Wherever several
scalars' ghosts are made, `ghosts_update ()`, a refinement or a coarsening
(`foreach_scalar_ordered`), a scalar comes after those it depends on, a stale
one of those updated with it; and marking one stale marks those that depend
on it. What is depended on may be scalars, vertex scalars, vectors, face
vectors or tensors, `(vector n, scalar alpha)`, flattened to their scalars
by the binding, and of the types the operator declares: cwrap checks where
it sees the declaration, `ghost_ops_bind ()` as the binding runs. The scalars
depended on may be local: freeing one while a scalar in
use depends on it stops the run (`ghost_ops_freeing ()`). To bind a coupled
operator for a while, keep `scalar_get_ops (s)` and give it back with
`scalar_restore_ops ()`.

**Operators that fill several scalars.** A coupled operator may also write
the scalars it depends on, at the same cells: libvof's volume fraction `c`
fills its interface's normal `n` and position `alpha` with it, at level jumps
(`prolong_plane`) and on parents (`restriction_plane`). Theirs then do
nothing: `prolong_none` and `restriction_none` (given). In one ghost update
every restriction runs first, finest level first, then the ghosts level by
level, each scalar after those it depends on, so `n`'s parents and jump
ghosts are `c`'s; the walls of a scalar whose prolongation is `prolong_none`
are made again after the prolongations that depend on it, which may have
written what they mirror. Their ghosts at level jumps are made only with
`c`'s: update `c` with them (`tests/test_ghosts.c`, 11). The null operators
are never called: no parent is visited for a `restriction_none` scalar, and
no family or level-jump region for a `prolong_none` one; `adapt_wavelet ()`
refuses a `prolong_none` scalar, which has no wavelet of its own.

**Derived scalars** are not prolonged. Fill the ghosts of the scalars they
come from, then compute them on the halo, `block_cells (shape, GHOSTS)`, as
Parthenon's `FillDerived` does. With cwrap, compute them on the interior,
`foreach () d[] = f (u[]);`, and let the stale ghosts be filled as for any
scalar (below): one more update, but right at every kind of boundary.

**Refining and coarsening keep ghosts valid.**
- `tree_refine ()` prolongs the children from the parent (interior and
  ghosts), then fills the children's ghosts.
- `tree_coarsen ()` fills the ghosts of the new leaf, and of the finer
  leaves that touch it, which read its children before.

Every other leaf reads the same values as before. So a balance cascade
never reads a stale halo, as long as the ghosts were valid before. Call
`ghosts_update ()` after changing values; `adapt_wavelet ()` does it for
every scalar first.

**Layout:** compact by default; `-DBAMR_ALIGN=V` starts every interior row
on a multiple of V Reals, for aligned vector loads. Code that needs the
interior coordinates of a cell takes `ci - shape.ox`, `cj - shape.g`, and
`block_at ()` goes the other way.

**Cost** (3D, m = 16, g = 2, 568 leaves with one level jump; serial, -O2,
on a loaded machine, so rough):

| | stored halos | gather (before) |
|---|---|---|
| fill, restriction included | 23–31 µs/leaf | 14–23 µs/leaf (restriction only) |
| 7-point stencil sweep | 12–16 µs/leaf | 16–25 µs/leaf, gather included |
| one fill + one sweep | 35–46 µs/leaf | 30–47 µs/leaf |

The fill costs about what the gather did, and buys linear (not piecewise
constant) ghosts at level jumps. Every further sweep that reads the same
ghosts then runs at stencil speed. The price is memory: 1.95× the interior
in 3D for m = 16 (1.42× for m = 32).

Not done yet: periodic boundaries, and flux correction (refluxing) at
coarse-fine faces, which conservation needs whichever way ghosts are
stored. See Future plans, "Ghost operators and adaptation".

**The same stencil with cwrap's `foreach ()`** (see `../libc-pp/README.md`,
"Stencils"):

```c
ghosts_update (u);
foreach ()
  v[] = u[1,0] + u[-1,0] + u[0,1] + u[0,-1] - 4 * u[];
```

cwrap expands it into `foreach_leaf ()` and `foreach_cell (shape, 0)`,
with one `block_data ()` per leaf for each scalar used (`const` for `u`, only
read) and `u[i,j]` as `cell_at`. It is serial over the leaves. In the body,
`x`, `y`, `z` and `Delta` are the cell's centre and size, from
`tree_get_origin (axis)` and `tree_get_size ()`, the domain that
`tree_shift_resize ()` sets.

**Stale ghosts** (`block-ghosts.h`): a bit per scalar, set while its ghosts
may be out of date. `foreach ()` marks the scalars it writes, and before a
`foreach ()` that reads one at an offset updates those marked, in one
`ghosts_update_dirty ()`; so the `ghosts_update (u)` above is no longer
needed after a `foreach ()` writes `u`. Every update clears the bits of the
scalars it updates, and a new or freed scalar starts clean. Writes through
`block_data ()` by hand mark nothing: call `ghosts_update ()`, or
`ghosts_mark_dirty ()`.

## Faces and fluxes

A `face_vector F` (cwrap) holds at each cell its lower faces: `F.x` at
i - 1/2. A block's faces along x are at i = 0 .. m, face m in the ghost
layer, so a face vector takes the same blocks as any scalar.

    foreach_face ()                               /* every face, once */
      F.x[] = (s[] - s[-1]) / Delta;              /* rotated for F.y, F.z */
    foreach () {
      double div = 0;
      foreach_dimension ()
        div += (F.x[1] - F.x[]) / Delta;          /* what leaves the cell */
      ...
    }

- **Each face is made once,** by its owner. A block makes its faces 0 to
  m - 1, and face m too where no node of its level is above it (a wall,
  a coarser leaf). So at a level jump the fine side makes the faces.
- **Fluxes agree at level jumps (refluxing), automatically.** A loop that
  writes F marks it stale, as for any scalar. The next loop that reads it,
  at an offset or at its cell, has its ghosts made first, by cwrap. For a
  face scalar, the ghosts are made as for a cell scalar, except for the
  faces where two blocks meet:
  - the parents' faces are the mean of their children's
    (`restriction_face`), including their upper face m;
  - a block keeps the face m it made;
  - a coarse leaf takes the face it shares with a finer region from that
    region's parent.

  Every coarse face is then the mean of the fine faces on it, so what
  leaves one cell enters the next. The sum of the net fluxes over the
  leaves is zero to round-off (`tests/test_faces.c`).
- **Refining** gives new faces by `prolong_face`. A fine face lying on a
  coarse one takes its value, one across a coarse cell's middle the mean
  of its two faces. **Coarsening** restricts by `restriction_face`. These
  are the defaults of a face scalar (by `scalar_location ()`).
- In `foreach_face ()`, `x`, `y`, `z` are the face's centre.
  `foreach_face_x ()` (`_y`, `_z`) loops over one axis only.

## Vertices

A `vertex_scalar v` (cwrap) holds at each cell its lower corner. A block's
vertices are at 0 .. m along each axis, m in the ghost layer.

    foreach_vertex ()
      v[] = zb_at (x, y);                       /* x, y: the vertex */
    foreach ()
      c[] = (v[] + v[1] + v[0,1] + v[1,1]) / 4; /* the cell's corners */

- **Each vertex has one owner,** the block where its indices are all below
  m, and the block on an upper wall for its vertices there.
  `foreach_vertex ()` visits those.
- **Every other copy gets the owner's value** when the ghosts are made,
  lazily, as for any scalar. Within a level that is a copy. Where a coarse
  leaf is above fine blocks, their vertices on it are the coarse ones:
  copied where they coincide, interpolated between (hanging vertices,
  `prolong_vertex`). Where fine blocks are above, the coarse leaf takes
  their parent's, the fine ones (`restriction_vertex`, an injection). So a
  vertex field is continuous at level jumps, and a linear one exact
  everywhere, refined and coarsened too (`tests/test_vertices.c`).
- The ghosts beyond a wall are made as a cell scalar's (a mirror of the
  cells), not a vertex's: read them with care.

## Tree caches

The loops over the tree (`foreach_leaf ()`, `foreach_leaf_parallel ()`,
`foreach_node ()`, the ghost update, restriction, block allocation) read
lists of the tree's nodes, `tree_caches ()` (`tree.h`): every node, the
leaves, and each by level, in Morton order. They are made again, by one
walk, only when read after the tree changed (`tree ()`, `tree_refine ()`,
`tree_coarsen ()`, `tree_destroy ()`): a `foreach_leaf ()` no longer walks
the tree (0.37 µs over 687 leaves in 3D, from 7.1 µs). The tree must not
change during such a loop: that stops the run. See `traversal-caches.md`.

## MPI: partition, plans, pruning (tree-mpi.h)

Every process builds the same whole tree. Then rank 0 partitions it and
tells the others; each process works out what it exchanges for the ghosts,
and frees the subtrees it does not need. The partition, the plans and the
pruning are computed from the tree with no MPI call, so they run and are
tested in one process for any number of processes. With `-D_MPI`,
`tree_mpi_partition ()` applies them on `MPI_COMM_WORLD`, and
`ghosts_update ()` exchanges blocks.

- **Ranks** (`tree.h`). `tree_npe ()` and `tree_rank ()` read two statics
  of `tree.c`. With `-D_MPI` they are `MPI_COMM_WORLD`'s, read by `tree ()`,
  which starts MPI if the program has not (`MPI_Init_thread`,
  `MPI_THREAD_FUNNELED`) and finalises it at exit; without, constants: 1
  and 0. Ranks are 0-based, as MPI's.
- **The owner, `node->pid`.** The root gets `tree_rank ()`, children their
  parent's pid, so a tree made on a process is that process's until it is
  partitioned. `ghosts_update_list ()` makes the ghosts of, and restricts,
  the nodes with `pid == tree_rank ()` only, with or without MPI (on one
  process, every node).
- **Loops over this process's leaves.** `tree_caches ()->owned` lists the
  leaves with `pid == tree_rank ()`, in Morton order, made again when the
  tree or the owners change. `foreach_leaf ()`, `foreach_leaf_parallel ()`,
  the loops over faces and vertices, and so `foreach ()`, walk it;
  `foreach_level (l)` and `foreach_level_or_leaf (l)` skip the nodes owned
  elsewhere. On one process nothing changes.
- **Reductions** (cwrap's `reduction (op:v)`, libc-pp README): after the
  loop, `tree_reduce (&v, type, op)` combines v over the processes
  (`MPI_Allreduce`, `+ * min max` of int, long, long long, their unsigned,
  float, double), rank 0's starting value counted once. `tree ()`
  distributes at once, so a sum counts the tree once; with
  `tree_replicated ()`, before the partition every process owns the whole
  tree, and a sum counts it npe times.
- **Partition.** The leaves are cut along the Morton curve into npe
  contiguous pieces by weight. Leaf q goes to `floor (npe (W_q + w_q / 2) /
  W)`, the piece that holds its middle; with unit weights this is computed
  in integers, and the pieces differ by one leaf at most. A parent goes to
  the process of its first leaf. A refinement, a coarsening or a pruning
  makes the partition stale (`tree_partition_valid ()`); a pruned tree is
  refused.
- **The pieces by their first keys.** `tree_morton_key (p)` (`tree.h`) is
  p's place on the curve, its first cell's bits interleaved at level
  `BAMR_MAX_DEPTH`. The partition keeps each process's first key
  (`key_first`, npe + 1, an empty process taking the next one's), and
  `tree_partition_owner (p)` finds the owner of any index by a binary
  search of them, without the node: what a process will ask about nodes it
  does not hold (milestone 2, the halo).
- **Neighbours as bit maps.** `tree_neighbours (p)` gives three masks of
  3^D - 1 bits, numbered as `block_ghosts ()` numbers its regions: `same`
  (a leaf of p's level), `fine` (a parent of p's level), `coarse` (no node
  there: the leaf one level coarser covers it). A wall is no bit. The 2:1
  balance makes these the only cases. `pid [dir]` is the rank of the node
  read there (that node, or the coarse leaf), -1 outside.
- **The plan.** A process restricts its own parents and makes its own
  nodes' ghosts, in the order of one process. Before each step it receives
  whole blocks (interior and ghosts, every scalar of the update) from the
  owners, in three kinds of exchange per level:

  | exchange | before | what is received |
  |---|---|---|
  | `TREE_MPI_RESTRICT` | restricting level l | the children of its parents owned elsewhere (restricted already if parents) |
  | `TREE_MPI_GHOSTS` | the ghosts of level l | the nodes of its `same` and `fine` bits (cells), and of its `coarse` bits (level l - 1, their ghosts just made, for the prolongation) |
  | `TREE_MPI_FACES` | the same, face scalars only | the upper children of the parent below a leaf along each axis: the faces m that the reflux reads |

  Each exchange is one list per peer, sorted by level, then k, j, i, on
  both sides. So one message per peer of known length is enough, with no
  size sent first. There is one plan for leaves (`ghosts_update ()`) and
  one for every node (`ghosts_update_levels ()`). The plan keeps the
  whole tree's depth and coarsest leaf level: the levels the update walks
  and exchanges at, the same on every process after pruning.
- **Pruning, for real.** A process needs its nodes, what it receives, and
  their ancestors. `tree_mpi_prune (plan)` finds the largest subtrees with
  nothing needed below their root, and `tree_prune (p)` (`tree.c`) frees
  each: the nodes below p, their blocks, their entries in the level index,
  the counts (`tree_get_nblocks ()`, `tree_get_nleaves ()`, the levels of
  `tree_get_depth ()` and `tree_get_depth_min ()`). The root p stays, with
  its blocks and its pid: `BAMR_NODE_IS_PRUNED`, no children, not a leaf;
  in the lists of every node and by level, in none of leaves and parents.
  A parent that is received but whose children are not keeps its block
  only. A family whose one child is needed keeps all its children (they
  are allocated together). It cannot be undone; refining or coarsening a
  pruned tree stops the run.
- **`tree ()` distributes at once** (milestone 2). On several processes it
  refines to the first level with 4 npe leaves at least (below depth max)
  and calls `tree_mpi_partition (NULL)`; on one process it is as before.
  `tree_replicated ()` builds the tree without distributing: the whole
  tree on every process, to refine as wanted before `tree_mpi_partition
  ()` (the tests, an initial mesh), since a distributed tree cannot be
  refined yet.
- **`tree_mpi_partition (weight)`.** Rank 0 partitions, with its weights,
  and sends every node's pid (in Morton order) and the loads to ranks 1 to
  npe - 1 by `MPI_Bcast`; the others set the pids and record the same
  partition. A tree of another size than rank 0's stops the run. From
  there no process uses the whole tree: it prunes to its skeleton and
  makes its plan from it (below). `tree_mpi_end ()` frees the plan and the
  partition; the tree stays pruned.
- **The skeleton, active and inactive nodes.** A process holds its own
  nodes, their ancestors, and N (p) of each own node p: p's same-level
  neighbours, the coarser leaf covering an absent one, for a leaf the
  children of a neighbour with children that touch it (the reflux), p's
  parent and children (the restriction). That is everything the plans
  read, and its structure is the whole tree's: the owner of a neighbour
  is its `pid`, here. N () goes both ways (N (x) holds p when N (p) holds
  x, by 2:1 balance), so the ranks of N () of a process's nodes, its
  peers, name it too, and every exchange is point to point between peers,
  sizes first: no `MPI_Alltoall`, nothing of size npe. A node a process
  neither owns nor receives is inactive (`BAMR_NODE_IS_INACTIVE`, beside
  `BAMR_NODE_IS_PRUNED`): its flags and pid, but no data. Blocks go by
  family: a family with an active member has them all, one with none has
  none (`grid_scalar_family_blocks ()`), and a scalar made later gets
  blocks in active families only. A node a family brings, out of the
  skeleton, is of a structure unknown here: pruned.
- **The plan from the skeleton:** what a process receives, from its own
  nodes and the tree here (`plan_needs ()`); what it sends, as its peers
  ask it (each sends each peer its receive lists, 7 ints a node). The
  levels are the whole tree's (`MPI_Allreduce`).
- **`tree_mpi_repartition (weight)`** (milestone 3), on the distributed
  tree, collective: weights for the leaves of this process
  (`tree_caches ()->owned`, NULL: 1 each).
  1. **The new cut by prefix sums:** the leaves and the weight before this
     process's (`MPI_Exscan`), the totals (`MPI_Allreduce`); each leaf to
     the process `tree_partition ()` would give it on the whole tree (in
     integers for unit weights). Each process's first key, load and leaves
     by one `MPI_Allreduce` each (min of the keys of its leaves).
  2. **Migration:** every node of this process that changes hands (by the
     new keys), sent to its new owner, whole: its index, whether a leaf,
     the block of every scalar in use, and N () of it as structure only,
     so that the new owner holds the skeleton around it. The new owner
     makes it (`tree_insert ()`, inactive until it has data) and copies the
     blocks in. Senders and receivers from the two partitions' first keys:
     to the new owners of a process's nodes, from the old owners of its
     new piece of the curve, both contiguous ranges of ranks.
  3. **Then the skeleton again:** every node's owner from the new keys
     (`key_first`, kept for this), pruning to the skeleton, the plan,
     active and inactive nodes.

  It returns the nodes this process sent. Integer weights give the whole
  tree's cut bit for bit; other weights may differ at a tie, the sums
  being made in another order.
- **Families of leaves never split** (milestone 4, on one process too).
  Each leaf goes by its middle as above, then every family whose children
  are all leaves goes whole to the process of its child 0. The pids still
  never decrease along the curve, and a process's load is within one
  family's weight of the mean. Merging makes new families of leaves that
  may span processes; `tree_mpi_repartition ()` gathers them: the owner of
  child 0, knowing from its skeleton that the family is all leaves (the
  siblings are in N () of each other), tells its peers the family's
  process, one round. So a family is always one process's at adapting.
- **Adapting the distributed tree** (milestone 4): `adapt_wavelet ()`
  with a plan applied.
  - The parents' plan receives every node's parent too, with its ghosts:
    what the wavelet reads.
  - Each process flags its own leaves, as on one process.
  - **Refining, in rounds** (`tree_mpi_refine ()`): a leaf is divided
    (`tree_divide ()`: children prolonged, no balance, no ghosts) once
    every index of its level beside it holds a node. A coarser leaf
    covering one is divided first: here if this process's, else asked of
    its owner (a peer), and the leaf waits a round. Each round, what was
    divided goes to every peer, which divides its copy (`tree_graft ()`:
    leaves, inactive, no blocks). The rounds end when no process has
    anything waiting or sent (`MPI_Allreduce`), a few a level at most.
  - **Coarsening, level by level, finest first** (`tree_mpi_coarsen ()`):
    each process merges its own families (`tree_merge ()`) where 2:1
    balance allows, from its skeleton (the touching children of the
    neighbours are in N () of the family's children), and sends the
    merges to its peers (`tree_ungraft ()`). A level's merges can let a
    coarser family merge, as on one process; merges of one level do not
    interact.
  - Then the skeleton and plan again (`tree_mpi_redistribute ()`), and
    `tree_mpi_repartition (NULL)`. A leaf out of the skeleton (a sibling
    a family brings) is forgotten at each skeleton (`tree_forget ()`:
    pruned, unknown): no one would tell this process when it changes, so
    it must not claim to know it.
  - **What comes out:** the tree and every cell-centred scalar as on one
    process, bit for bit (the order of the rounds does not matter: a leaf
    is prolonged from its own cells and ghosts as they were at the start,
    a family merged from its children's interiors). A face or vertex
    scalar's restriction reads the children's ghost layer, which one
    process updates as it coarsens family after family, so its merged
    values depend on that order: on several processes they may differ in
    the last bits.
- **`tree_coarsen ()`'s 2:1 check, exact** (milestone 4, on one process
  too): coarsening p is refused only if a child of a neighbour that
  touches p has children (`tree_child_touches ()`). It refused before if
  any child of a neighbour had children: more than 2:1 needs, and more
  than a process's skeleton knows. One process now coarsens a little more
  often; `tests/bit/compare.sh` differs for adapting runs.
- **In `ghosts_update_list ()`** (`#ifdef _MPI`, with a plan applied): the
  exchange before each restriction level and each ghost level, at the
  plan's levels. Without `_MPI`, every node has pid 0 = rank: nothing
  changes but one comparison per block.

### Tests

`tests/test_partition.c`, in one process, 1-3D, double and float:

| what | against |
|---|---|
| ranks, owners | `tree_npe () == 1`, `tree_rank () == 0`; after `tree ()` every pid 0; a refined leaf's children get its pid |
| unit-weight cut | the closed form `floor (npe (2 q + 1) / (2 n))`; pieces within one leaf; none empty when n >= npe |
| weighted cut | loads summed again from the pids; max <= W / npe + w_max, and <= the optimal contiguous cut (bisection with greedy packing) + w_max |
| parents | the smallest pid of the leaves inside the parent's box, by search |
| keys, owners, this process's leaves (every cut above) | `tree_morton_key ()` increasing along the leaves, a node's its first leaf's (found by descending); `tree_partition_owner ()` of every node its pid; `owned` the leaves of pid 0, in order |
| bit maps and ranks, every node and direction | the leaves that overlap the neighbour's box on the finest lattice; the pid of the node read, by search |
| plans for npe = 1, 4, 7 | what each process must receive, from that geometry, from the owners; every send list equal, item for item, to the matching receive list |
| pruning, each process of npe = 1, 4, 7 (sphere), 2^D + 3 (uniform: some own nothing) | the tree made again for each; a node kept iff it is the root or its parent has something needed strictly below, by a search of boxes; every kept node at its place with its flags and pid, the pruned ones `IS_PRUNED` without children; every other node gone from the level index; `tree_pruned ()`, the node and leaf counts, the depth and coarsest level (0 with no leaf left), the lists; the partition stale and refused; a pruned node not pruned again; a scalar made and freed, `grid_scalar_check ()` |
| bad arguments, staleness | npe 0, weight 0 or NaN refused, partition kept; stale after refining |

`tests/test_mpi.c`, with `-D_MPI`, MPI started by `tree ()` and
finalised at exit: a cell scalar, a face vector and a
vertex scalar on a tree refined about a sphere. The references are
`ghosts_update ()` on the whole tree, before the partition (each process
owns every node then). Each process also makes the partition itself, to
compare with what rank 0 sends, and then gives every node back its own
rank. After `tree_mpi_partition ()`, every node needed is still there, and
the pids and the partition are those made locally, its first keys too;
the owner of every node of the whole tree, freed here or not, comes from
the keys alone; each process's `owned` list is its kept leaves of its
rank, in order, and over the processes every leaf once. Reductions in
`foreach ()` (`+` of cells and levels, `*` with one factor 3 at the
corner cell, `min` and `max` of the level), `foreach_face ()` (`+`, by
axis) and `foreach_level (2)` give the whole tree once, against the same
loops before the partition, which gave it npe times. Then every kept block
a process does not own is poisoned (NaN), and so are the ghosts of its own
leaves and its own parents. After `ghosts_update ()`, every cell that the
update writes in an owned node must equal the reference bit for bit, for
leaves and for every node; a block read but pruned would be a crash under
ASan. Last, `tree_mpi_end ()` leaves the tree pruned, and an update with no
plan runs. Runs in 1-3D, double and float, ASan and UBSan, on 1, 2, 3, 4,
7, 10, 13 and 16 processes of `tests/mpi-fork`: all pass (3D in double on
16 processes needs more than this machine's 11 GB free, each rank holding
four whole-tree snapshots: one rank is killed; in float it passes).

Deliberate bugs, one line each (`test_partition` failures, then
`test_mpi` on 3 processes, in 2D and 3D; "crash": a rank stops, by ASan
or a missing message):

| bug | partition 2D, 3D | mpi 2D, 3D |
|---|---|---|
| pruning ignores what is received | 66, 66 | crash, crash |
| nested roots pruned (not the largest only) | 11, 11 | 0, 0 |
| freed leaves left in the level counts | 3, 3 | 0, 0 |
| a pruned node keeps `HAS_CHILDREN` | crash, crash | 0, 0 |
| children not given their parent's pid | 1, 1 | 8, 8 |
| ranks > 0 do not set the pids received | - | crash, crash |
| a received parent keeps its children | 26, 24 | 0, 0 (more kept than needed) |
| pruning leaves the partition valid | 14, 14 | 0, 0 |
| `pid` of a coarse neighbour -1 | 1, 1 | 0, 0 (the plan reads the node) |
| every node processed (`BLOCK_MINE` 1) | - | 0, 0 (wasted work, overwritten by what is received) |
| the local coarsest level for the restriction, not the plan's | - | 0, 0, also on up to 16 processes |

Before this rewrite, the plans (reflux children, coarse leaf sent early,
restriction one level short or first child only, `fine` and `same`
swapped, exchanges after the level or dropped) were caught as in the
previous commit's table; those checks are unchanged. The last bug changes
the exchanges only on a process whose kept leaves are all finer than the
tree's coarsest: none of the test's trees has one.

`tests/run.sh` runs it on `tests/mpi-fork/`, a stand-in `mpi.h`, so it
needs no MPI installed. The stand-in forks the ranks in `MPI_Init ()` and
joins each two by a socket pair. Every message carries its tag and length,
and a receive that does not match them stops the run. It has
`MPI_Init_thread`, `MPI_Initialized`, `MPI_Finalized`, `MPI_Bcast`, and
`MPI_Allreduce` of `+ * min max` over int, long, long long, their
unsigned, float and double. Rank 0's `MPI_Finalize ()` fails the run if
another rank failed, also when called at exit.

Milestone 1 (owned leaves, reductions, MPI in `tree ()`, first keys):
`test_mpi` passes in 1-3D on 1, 2, 3, 4, 7, 10, 13, 16 ranks in float and
1, 2, 4, 7, 13 in double. Deliberate bugs, one line each (failures in
`test_partition` / `test_mpi` on 3 ranks, 2D and 3D alike):

| bug | partition | mpi |
|---|---|---|
| `owned`: every leaf | 27 | 6 |
| `foreach_leaf ()` over every leaf | - | 3 |
| faces over every leaf | - | 3 |
| a level loop over every pid | - | 3 |
| the owners changed, the lists kept | 29 | 0 (pruning rebuilds them) |
| an empty process's first key 0 | 10, 6 | 0 |
| the owner's search off by one | 27 | 3 |
| a key's bits interleaved the wrong way | 31 | 3 |
| `tree_reduce ()`: min and max swapped | - | 3 |
| `tree_reduce ()`: no `MPI_Allreduce` | - | 3 |
| `tree ()` does not start MPI | - | 1 |

The min/max swap passed until the test reduced the rank itself: every
piece touches the coarsest and finest levels and both ends of x.

Milestone 2 (`tree ()` distributing, the halo by queries, the plan from
the local tree). `test_mpi` builds its tree with `tree_replicated ()` and,
before distributing, takes an oracle from the whole tree: the plan of
`tree_mpi_plan_new ()`, what is needed, and so what a process keeps.
After `tree_mpi_partition ()`, the plan made from the queries must equal
it entry for entry (both sides, every exchange, peer and index), and its
levels too; the nodes kept, by index, the oracle's; each needed node a
leaf, a parent or pruned as there; the pids; `tree_pruned ()` and the leaf
count against the flags. The nodes freed and made again are other
objects: the test finds them again by index, and poisons the new ones too
before the bit-for-bit ghost comparison. Passes in 1-3D, float, on 1, 2,
3, 4, 7, 10, 13, 16 ranks. `tests/test_tree_mpi.c`: `tree (5, 4)` itself,
1-3D, double and float, 1 to 16 ranks: every leaf of a process at the
first level with 4 npe leaves, the level's leaves once over the processes
and within one of the mean, owned as the keys say; an integer of the
global cell index exact in every ghost inside the domain after
`ghosts_update ()`; a reduction over the cells the level once.

| bug (one line) | mpi np 3, 2D, 3D | tree_mpi np 3, 5, 2D, 3D |
|---|---|---|
| queries for the own leaves only (not the parents) | crash | crash |
| nothing there: the index made, not the coarser leaf | crash | 0 |
| no reflux children made | 1, 1 | 0 (uniform: no reflux) |
| a leaf made pruned | 3, 3 | 0 |
| the owners of the nodes made not set | crash | crash |
| siblings made leaves | crash | crash at np 5 (2D) |
| the pruned count kept when a stub has children again | 3, 3 (np 7: 7, np 16: 16) | 3, 5 |
| sends not taken from the requests | crash | crash |
| the plan's levels local, not the whole tree's | 0; crash at np 16 (3D) | 0 |
| `tree ()` a level short (2 npe) | - | 15 (np 5, 2D), 9 (np 3, 3D) |
| no pruning to the own nodes first | 9, 9 | 0 |
| a parent answered as a leaf | crash (2D), 3 (3D) | 0 |
| the owner from the asking node, not the index | crash | crash |

"crash": a rank stops (an assertion, ASan, a node asked of a process that
does not own it, or a message the stand-in did not expect).

Milestone 3 (`tree_mpi_repartition ()`). `test_mpi` takes a second oracle
from the whole tree, with integer weights from each leaf's key (7 or 8 on
the second half of the curve, 1 or 2 on the first, so that the cut moves
far). After the first distribution it fills a field `d` on every node of
each process (`foreach_level ()`, parents too) and poisons every other
block of it. A repartition with unit weights must move nothing and leave
the first oracle's plan; one with the weights must give the second
oracle's cut (loads, leaf counts, first keys, imbalance), pids, nodes and
plan entry for entry; `d` exact on every node of its new owner; every
other scalar on the leaves as the last update left it (made again by
`fill ()` and the same update); and the ghosts bit for bit. 1-3D, float
on 1, 3, 7, 13, 16 ranks, double on 2, 5, 13 (3D: 5): all pass. Moved,
of the whole tree's nodes: 2D 85 of 197 on 3 ranks, 174 on 16; 3D 222 of
521 on 3, 465 on 16, with 44 million cells compared.

| bug (one line) | mpi np 3, 7 (2D, 3D) | 1D np 13 |
|---|---|---|
| `MPI_Exscan` ignored: every rank counts from 0 | 35, 83; 19, 83 | |
| the weights ignored | 21, 48; 21, 49 | |
| the cut by the running sum, not the leaf's middle | 16, 38; 0, 34 | |
| an empty process keeps its own key | 0 | 131 |
| the parents not sent | 3, 7 | 13 |
| the blocks not copied in | 3, 7 | |
| one scalar short | 3, 7 | 13 |
| a leaf received made pruned | crash | |
| the pids by the old partition | crash | |
| the halo not made again | crash | |
| the imbalance from the local load | 5, 13; 2, 12 | |

The empty process, the parents and the scalar left out were caught by
nothing until the test had more processes than leaves (1D, 13 ranks, in
`tests/run.sh` now) and checked every node and every scalar.

The skeleton (no queries, no `MPI_Alltoall`, inactive nodes). The oracle
marks the skeleton on the whole tree: the own nodes and N () of them, by
a search of its own. After each distribution and repartition: the nodes
kept those of the skeleton and their families; a skeleton node a leaf, a
parent or pruned as on the whole tree (another may be pruned, its
structure unknown); active iff the process owns or receives it; blocks
iff its family has an active member, also for a scalar made after
(`grid_scalar_check ()`, which now counts blocks in wholly inactive
families as an error); every node the plan reads in the skeleton; the
peers symmetric (a 0/1 row per rank through `MPI_Alltoall`, in the test
only) and every rank of the plan a peer. An inactive node has no block,
so a read the plan forgot is a NULL read, which ASan reports. 1-3D on 1,
2, 3, 4, 7, 13, 16 ranks: all pass. A process keeps 0.72, 0.51, 0.51,
0.37 of the 197 nodes (2D, 2, 4, 7, 16 ranks), and 0.71, 0.51, 0.49, 0.33
have blocks; 3D 0.82, 0.66, 0.65, 0.45 kept, 0.81, 0.66, 0.65, 0.43 with
blocks: little saved, halo families being mostly mixed (see "Further
options").

| bug (one line) | np 3, 7 (2D) | np 3, 7 (3D) | 1D np 13 |
|---|---|---|---|
| N (): no coarse leaf for an absent neighbour | 0, 0 | 0, 0 | crash |
| N (): no touching children (reflux) | 2, 14 | 0, 1 | crash |
| N (): the far half taken as touching | 1, 3 | 0, 1 | crash |
| N (): no parent | 0, 0 | 0, 0 | 0 |
| N (): no children | 0, 0 | 0, 0 | 0 |
| received nodes left inactive | crash | crash | crash |
| a new scalar in inactive families too | 2, 15 | 0, 3 | 31 |
| inactive families keep their blocks | 2, 28 | 0, 6 | 58 |
| no structure sent with a moving node | crash | 16, 40 | crash |
| empty old pieces taken as senders | 0, 0 | 0, 0 | 0 |
| no pruning to the skeleton | 18, 42 | 18, 42 | 78 |

Not caught: the parent and children rules of N (), and the coarse leaf
but in 1D. On these trees those nodes and ranks are in the skeleton and
among the peers by other relations too (a coarse leaf is also a neighbour
of the own node's parent; a parent's owner is a neighbour). The rules stay
as what the symmetry argument needs. "Empty old pieces" needs a new piece
of the curve spanning an empty old one, which no test tree makes. The pruned
count was caught by nothing until the test counted the flags. The local
levels differ only where a process lacks the finest or coarsest leaves:
none of the test's trees but 3D on 16 ranks has one.

Milestone 4 (adapting). `tests/test_adapt_mpi.c`: each process builds the
whole tree (`tree_replicated ()`), refined about a sphere, and adapts it
five times as one process does (`adapt_wavelet ()`, maxlevel 5, 4 in 3D):
three steps to a sharp front (`tanh`, width 0.01) whose centre moves,
which refine and coarsen, then two on a constant, which merge every family
down to minlevel (and make families of leaves that span processes). The
reference: the counts of each step, and every leaf with every cell of two
cell scalars (the adapted field and a linear one carried along) after a
ghost update, after the front and at the end. Then the same tree
distributed (`tree_mpi_partition ()`) and adapted the same way: the
counts of each step the reference's; after each step every process within
one family of the mean; at both checkpoints every leaf owned once and
every cell of both scalars (interior and ghosts) the reference's, bit for
bit. 1-3D on 1, 2, 3, 4, 7, 13, 16 ranks: all pass (2D: 304 leaves and
243 200 cells after the front; 3D: 1212 leaves, 19.4 million cells).
`test_adapt.c` adds the exact 2:1 check: a family beside a node whose far
child has children is coarsened, the tree stays balanced; once the near
child has children, not.

The test found two faults before it passed: the cut could leave a family
of leaves split after merges (the families were whole only when the cut
made them; now the repartition gathers them), and a process could hold a
sibling out of its skeleton as a leaf after its owner had divided it,
then trip on its children arriving with a move (now forgotten at every
skeleton).

| bug (one line) | adapt 2D, 3D | partition 2D | adapt_mpi np 3 2D, 3D; np 13 2D | mpi np 3 |
|---|---|---|---|---|
| coarsen: the far half taken as touching | 4, crash | 0 | 9, 9; 39 | 0 |
| coarsen: no 2:1 check | 2, crash | 0 | crash | 0 |
| the cut: families split | 0 | 35 | 3, 0; 0 | 31 |
| the repartition: split families not gathered | 0 | 0 | 0, 6; 26 | 0 |
| the plan: parents not received | 0 | 2 | 6, 6; crash | 0 |
| refine: a remote coarser leaf not asked | 0 | 0 | crash (no end) | 0 |
| refine: peers not told (no graft) | 0 | 0 | crash | 0 |
| coarsen: peers not told (no ungraft) | 0 | 0 | crash | 0 |
| coarsen: no 2:1 check under MPI | 0 | 0 | crash | 0 |
| coarsen: all levels in one round | 0 | 0 | 9, 3; 39 | 0 |
| a leaf out of the skeleton not forgotten | 0 | 0 | 0, 0; crash | 0 |
| no repartition after adapting | 0 | 0 | 3, 3; 13 | 0 |

The missing repartition was caught by nothing until the test checked the
balance after each step. A remote family taken as a candidate changes
nothing (its leaves are never flagged here), so its check is not listed.

On real ranks, with OpenMPI (`sudo apt install openmpi-bin
libopenmpi-dev`):

    for d in 1 2 3; do
      ../libc-pp/build/cwrap -dimension=$d -D_MPI -D_GNU_SOURCE -std=gnu11 \
        -O1 -g -Wall -Wextra -Wno-unused-parameter \
        -fsanitize=address,undefined $(mpicc.openmpi -showme:compile) \
        -I. -I.. tests/test_mpi.c tree.c tree-mpi.c grid-scalar.c \
        memindex.c voxel-memindex.c nano-vdb.c pool.c block-ghosts.c \
        adapt.c ../libutils/bits.c $(mpicc.openmpi -showme:link) -lm \
        -o /tmp/test_mpi_$d
      for np in 2 3 4 7; do
        ASAN_OPTIONS=detect_leaks=0 mpirun.openmpi --oversubscribe \
          -n $np /tmp/test_mpi_$d | tail -1
      done
    done

Not run on OpenMPI since the rewrite (not installed on this machine).
Ubuntu's mpich (4.2) does not run more than one rank: it is built for
PMIx, but its launcher, Hydra, is PMI-1 only, so every rank starts alone
("0 of 1"). Use OpenMPI, which installs beside it, through the `.openmpi`
names.

Measured in 2D, 197 nodes (3 rounds of refinement about the sphere): with
2 processes, 40 blocks received for leaves (54 for every node); with 3
processes, 92 (137). The pieces' surface in `test_partition`, as the
fraction of face pairs between leaves that are cut: 3D, 1576 leaves, 0.061
with 4 processes, 0.188 with 16. The nodes a pruned process keeps, of the
whole tree's: `test_partition`, 0.51 (2D) and 0.51, 0.47 (3D) with 4 and
7 processes; `test_mpi`, 2D 0.72, 0.69, 0.51, 0.51 with 2, 3, 4, 7
processes, 3D (521 nodes) 0.82, 0.82, 0.66, 0.65.

### Not done yet

What the distributed tree cannot do yet, most important first.

**Wrong or missing under MPI**

1. **Multigrid.** `ghosts_update_level ()` and the level loops make no
   exchange, so `mg_helmholtz ()` (liblinear-algebra) is wrong on a
   distributed tree. The level loops already skip other processes' nodes,
   and the residual norm is a reduction; the per-level exchange is the
   missing piece (milestone 5).
2. **`block_ghosts ()` called alone** makes no exchange: only
   `ghosts_update_list ()` and what calls it do.
3. **Refining or coarsening by hand:** `tree_refine ()` and
   `tree_coarsen ()` stop the run on a distributed tree; only
   `adapt_wavelet ()` adapts it.
4. **Coarsened face and vertex scalars** may differ from one process in
   the last bits: their restriction reads the children's ghost layer,
   which one process updates family after family (see "Adapting the
   distributed tree"). Cell scalars come out bit for bit.
5. **Reductions** count the tree npe times before a partition with
   `tree_replicated ()` (every process owns all of it then), and writes
   through a pointer, a member or a subscript in a stencil loop
   (`*e = ...`, `same [k] += ...`) are not reductions: serial only.

**Not tested where it will run**

6. **Real MPI:** since the rewrite, only on the stand-in `tests/mpi-fork`,
   up to 16 processes on one machine; not on OpenMPI ranks, nor on more
   than one node.
7. **The other libraries:** libvof not run under MPI; libhmesh (points
   that move between processes) has no MPI design yet; no output (libview)
   or checkpoint and restart of a distributed tree.
8. **Rules of the skeleton no test needs:** N ()'s parent and children
   rules, and its coarse leaf but in 1D (reached by other relations on
   the test trees); a new piece of the curve spanning an empty old one.

**Scale and speed**

9. **The start needs the whole tree on every process:** `tree ()`'s
   coarse uniform level, or a `tree_replicated ()` mesh refined first. An
   initial mesh too large for one process must be grown after
   distributing, by `adapt_wavelet ()`.
10. **O(npe) on each process:** the first keys (`key_first`), and the
    repartition's `MPI_Allreduce` of npe-long arrays (an `MPI_Allgather` of
    one key and one load per process would do).
11. **Communication rounds:** one exchange per level in a ghost update,
    and no computation overlapping it. See "MPI" in "Future plans"
    (dummy children for the halo, coarse buffers).
12. **Halo memory:** blocks go by family, so a mixed halo family carries
    unused blocks; blocks per node wait for the edge-cut partition (see
    "Block allocation", "Further options").
13. **Load weights:** given by hand, no cost measured per block;
    non-integer weights may not give the whole tree's cut at a tie.

**Later, as planned**

14. Periodic boundaries; the partition by edge cut; MPI with GPU offload.
    MPI is started `MPI_THREAD_FUNNELED`: no OpenMP thread may call MPI.

Suggested order: 1 and 2 (milestone 5), then 6 on OpenMPI ranks, then 7
for libvof.

## Tests

    QUICK=1 sh tests/run.sh      # after a change: about 40 s
    sh tests/bit/compare.sh      # every cell bit for bit against HEAD: 20 s
    sh tests/run.sh              # before a commit: the whole

`QUICK=1` runs the pool, then tree, blocks, ghosts, adapt and stencil in 1,
2 and 3 dimensions, double and float; the whole adds OpenMP, the aligned
layout, the other level indices and `MADV_FREE`. Each library source is
built once per build and shared by its tests, which run `JOBS` at once (by
default, the processors). `tests/bit/compare.sh [REF]` builds
`tests/bit/checksum.c` against the working tree and against `REF` (`HEAD`
by default): it refines, fills ghosts with given operators, adapts, and
hashes every cell of every block; `TIMING=1` also times
`ghosts_update_levels ()` (`tests/bit/timing.c`). A `REF` whose API
`checksum.c` does not build against fails, and says so.

`tests/run.sh` builds each test in `tests/` with cwrap (`CWRAP`, by default
`../libc-pp/build/cwrap`), in 1, 2 and 3 dimensions, in double and float
(`-DBAMR_SINGLE`), under AddressSanitizer and UBSan, and fails on any warning:
- `test_pool.c`: the hot and cold free lists, with page faults counted;
- `test_tree.c`: counts, the scalar hooks, Morton order, refinement, and the
  loops of `tree.h`; the tree's lists against a walk after refining and
  coarsening, and a change of the tree in `foreach_leaf ()` stopping it;
- `test_blocks.c`: every loop of `block-loops.h` against `block_cells`, in
  both layouts and every halo, and `block_at ()`;
- `test_ghosts.c`:
  - a linear field, extrapolated at the walls: every ghost is exact, at
    level jumps, walls and corners, and on parents with
    `ghosts_update_levels ()`;
  - `bc_neumann` and `bc_odd` mirror over every wall;
  - injection gives the covering coarse cell;
  - minmod is exact for a linear field and stays positive at a kink;
  - a custom restriction runs once per parent;
  - a derived scalar computed on the halo matches filling its ghosts;
  - refining and coarsening keep every leaf's ghosts equal to a fresh
    `ghosts_update ()`;
  - a freed id gets the default operations back.
- `test_adapt.c`: families allocated at once, prolongation (exact for a
  linear field, and keeping the mean), coarsening, 2:1 balance through
  random refining and coarsening, and `adapt_wavelet ()` on a front: it
  settles, keeps the integral, and coarsens a constant field to minlevel;
  with minlevel 0, never into the root, and a lone root is refined again.
  The wavelet of a node is the difference between its cells and its
  parent's prolongation of them, by the scalar's own operator
  (`block_prolong_into ()`): `test_ghosts.c` checks it is zero for a lake
  at rest over a wavy bed.

- `test_stencil.c`: cwrap's `foreach ()` and `s[i,j,k]`: every interior cell
  once and no ghost written, stencils in every dimension equal to `cell_at`
  bit for bit, offsets that are expressions, `+=` and `++` on a cell, a
  vector's components, an array of scalars and a pointer to one, a scalar
  parameter, and names in the body (`c`, `level`, `i`, `node`, `sblocks`)
  that the loops' own names do not hide; `x`, `y`, `z` and `Delta` on a
  moved, resized domain with two levels of leaves, against the centre of
  each cell computed apart; ghosts updated lazily, against `ghosts_update ()`
  by hand; `foreach_neighbor (1)` against the same box summed by hand, its
  ghosts updated first, and `x` of the neighbours. Also with OpenMP and
  `-DBAMR_ALIGN=4`.
- `test_ghosts.c` also: a coupled prolongation (lake at rest, `h + zb` the
  surface at level jumps, walls and refined leaves, the bed made first
  against the order of ids), a stale bed making the depth stale, and a
  cycle, a wrong count, and `scalar_set_ops ()` of it refused; coupled
  operators that write what they depend on, local scalars with
  `prolong_none` and `restriction_none` (every ghost of the written one
  made, its children too), a binding given back by `scalar_restore_ops ()`,
  adapting on a `prolong_none` scalar and freeing a scalar depended on,
  refused; a vector dependency bound whole, each of its components depended
  on, and a binding of other types refused as it runs.

- `test_partition.c` and `test_mpi.c` (`-D_MPI`, on 1 and 3 processes of
  `tests/mpi-fork`): see "MPI: partition, plans, pruning".
- `test_morton.c`: `tree_morton_key ()`, interleaved by masks (log2 of
  the bits in steps, as libhmesh's `spread ()`), against the reference
  one bit at a time: every value of each axis at the finest level (the key
  is the OR of one part an axis), every index of the levels with at most
  2^18 cells, 10^6 random indices, and a child's key from its parent's
  (`key (p) | c << dimension (BAMR_MAX_DEPTH - level - 1)`). Deliberate
  bugs caught: a wrong mask in 2D and in 3D, a wrong shift in 1D.

  | key at -O2, ns | 1D | 2D | 3D |
  |---|---|---|---|
  | one bit at a time | 17.5 | 41.2 | 68.2 |
  | masks | 1.5 | 9.0 | 9.9 |
- `test_index.c`: the three level indices side by side (memindex,
  voxel-memindex, nano-vdb, dense and compact) against a reference hash
  map: random assigns, overwrites and frees, clustered and scattered; every
  entry and the indices around it, by get and by cursor, out of range too;
  emptying (nano-vdb gives back every node), destroying and reusing; levels
  on both sides of nano-vdb's hashed root (above 4096).

It also runs `blocks`, `ghosts` and `adapt` with OpenMP (2D, 3D) and with
`-DBAMR_ALIGN=4` (1D to 3D), and `tree`, `blocks`, `ghosts` and `adapt` on
the other level indices: `-DBAMR_INDEX=2` (voxel-memindex), `3` (nano-vdb)
and `3 -DNVDB_COMPACT` (2D, 3D).

### On the GPU

    GPU=1 sh tests/run.sh        # the offload builds alone, on the device

`GPU=1` builds tree, blocks, ghosts, adapt, stencil, faces and vertices in
1, 2 and 3 dimensions, double and float, with OpenMP offload to an NVIDIA
GPU (gcc's `nvptx-none` offload compiler, `gcc-13-offload-nvptx`), and runs
them with `OMP_TARGET_OFFLOAD=MANDATORY`: no device is a failure, not a run
on the host. `GPU_ARCH` sets the PTX target (default `sm_80`, the newest gcc
13 knows; CUDA 13 dropped gcc's default, `sm_50`). There are no sanitizers,
which offload does not take. `test_blocks.c` is the test with offloaded
code: `gpu_cells` against `block_cells`, with `memcmp`.

The device must give the host's bits:

| Setting | Why | Without it (RTX A3000, gcc 13.3, CUDA 13.5, driver 620.70) |
|---|---|---|
| `-ffp-contract=off` | gcc fuses `a * b + c` into one `fma.rn` on the device, by default; read from the host's compile, not from `-foffload-options` | blocks 3D: 32 of 261 checks fail, double and float |
| `GOMP_NVPTX_JIT=-O0` | gcc then writes `mul` and `add` without a rounding modifier, which the driver's PTX compiler may fuse (PTX ISA, `mul`, `add`) and does at -O1 and up | blocks 3D: 32 of 261 checks fail, double and float |

In 2D the stencil is `... - 4 u`, and a product by 4 is exact: fused or
not, the bits agree, which is why only 3D (`- 6 u`) shows it. The 32 are
every `gpu_cells` check in 3D: four m, four layouts, two halos; the cells
differ in their last bit (2 to 35% of them), and match a host computing
`fma (-6, u, sum)` bit for bit. Both settings cost speed, so a run wants
neither; `GOMP_NVPTX_JIT` from the environment overrides `-O0`. With both,
42 of 42 runs pass on the device; with no device (in a sandbox), the six
blocks runs fail on `MANDATORY`, and the others, with no target region,
pass on the host.

The OpenMP builds of the whole suite take `-foffload=disable`: a gcc built
with `--enable-offload-defaulted` (Ubuntu's) offloads to every offload
compiler installed, and with `gcc-13-offload-nvptx` installed, 11 blocks
runs failed to link on `-fcf-protection` before it.

Measured on an RTX A3000 Laptop GPU (Ampere, compute 8.6). The device reads
host `malloc` memory here (the open kernel module's HMM) as `test_blocks.c`
maps nothing; without it, `block_table ()` and `target data map` are needed
(`block-loops.h`, 3. Which one).

### Profiling: the CPU against the GPU

    sh tests/run2.sh             # 2D and 3D, double and float
    GPU=0 sh tests/run2.sh       # the CPU only

Not a test: `tests/bench_gpu.c` times a ghost update piece by piece, on
the CPU and on the GPU, each as fast as it goes. No bit identity: `-O3
-march=native`, fma allowed, the PTX JIT at its default.

- **CPU:** the library's functions on every thread: `block_restrict ()`,
  `block_ghosts ()`, `block_prolong ()`, level by level, as
  `ghosts_update ()` runs them. The fill is timed with the walls off (a
  boundary condition that does nothing) and the jumps off
  (`prolong_none`), then with the jumps, then with everything. The jumps
  and walls are the differences, so they are noisy where small (even
  negative).
- **GPU:** every node's block in one device array, put there once, with
  kernels of the benchmark's own.
  - Restriction: the mean of the children.
  - Prolongation: every parent's children, as a refinement fills them.
  - The stencil: `v = u + k lap (u)`.
  - The ghosts: one table entry per ghost cell of every leaf, built on
    the host from the tree as `block_ghosts ()` would fill it: a copy, a
    linear prolongation at a level jump, or a wall, followed back to the
    copy or prolongation it ends in. That makes two launches a level,
    with no pass for the walls.
  - Upload and download of every block: what a ghost update on the host
    would cost the GPU.
- **Checked before any timing:** device against host, to `64 REAL_EPSILON`:
  - the leaves' ghosts, after the device's are spoilt first;
  - the restricted parents;
  - the stencil;
  - one prolongation.

  Each of seven one-line breaks fails it, in 2D and in 3D:
  - a copy's source shifted;
  - a jump's sign inverted;
  - a wall's mirror off by one;
  - the wrong child restricted;
  - a prolongation's half inverted;
  - the stencil without its Laplacian;
  - the jumps left out.

Measured on 2026-11-15: Intel i9-11950H (8 cores, 16 threads), RTX A3000
Laptop GPU (sm_86, 6 GiB), gcc 13.3, CUDA 13.5. Times are in ms, the
median of 7 or more calls.

The trees, uniform at a base level and refined on a shell of radius 0.3:
- **2D:** levels 4 to 8, m = 32: 3016 leaves, 42 MB a scalar in double.
- **3D:** levels 3 to 5, m = 16: 4880 leaves, 357 MB a scalar in double.

The columns:
- **CPU:** 16 threads.
- **GPU:** the benchmark's kernels on `gpu_cells ()` as it is now: one
  flat simd loop (below).
- **old `gpu_cells`:** the same before that change, its four loops
  collapsed and no simd. This changes only restriction, prolongation and
  the stencil.

| 3D | double CPU | double GPU | CPU / GPU | old `gpu_cells` | float CPU | float GPU | CPU / GPU | old `gpu_cells` |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| restriction | 11.2 | 1.41 | 7.9 | 7.55 | 5.60 | 1.06 | 5.3 | 2.57 |
| ghosts: copies (15.5 M) | 27.7 | 3.22 | 8.6 | | 19.0 | 2.66 | 7.1 | |
| ghosts: jumps (3.2 M) | 4.13 | 1.19 | 3.5 | | 1.99 | 1.05 | 1.9 | |
| ghosts: walls (0.3 M) | ~0 | 0.29 | | | ~0 | 0.28 | | |
| ghosts, the whole fill | 31.8 | 4.45 | 7.1 | | 19.5 | 3.75 | 5.2 | |
| `ghosts_update` (restriction + fill) | 43.9 | 5.88 | 7.5 | 12.0 | 28.0 | 4.81 | 5.8 | 6.27 |
| prolongation (20.4 M cells) | 24.6 | 2.34 | 10.5 | 47.7 | 15.4 | 1.48 | 10.4 | 14.0 |
| stencil (20.0 M cells) | 36.3 | 3.52 | 10.3 | 47.3 | 19.2 | 1.80 | 10.7 | 12.4 |
| step (`ghosts_update` + stencil) | 80.4 | 9.40 | 8.6 | 59.3 | 48.4 | 6.61 | 7.3 | 18.3 |
| upload, every block | | 62.8 | | | | 32.0 | | |
| download, every block | | 63.4 | | | | 33.5 | | |

| 2D | double CPU | double GPU | CPU / GPU | float CPU | float GPU | CPU / GPU |
|---|--:|--:|--:|--:|--:|--:|
| restriction | 0.63 | 1.10 | 0.6 | 0.12 | 1.03 | 0.1 |
| ghosts, the whole fill | 0.88 | 2.25 | 0.4 | 0.34 | 2.22 | 0.2 |
| prolongation | 2.07 | 1.23 | 1.7 | 1.07 | 1.11 | 1.0 |
| stencil | 3.94 | 0.52 | 7.6 | 0.73 | 0.41 | 1.8 |
| step | 7.23 | 3.88 | 1.9 | 2.10 | 3.66 | 0.6 |

After the library's `gpu_cells ()` changed, run2.sh gave the GPU column
again to within 1% (3D double step 9.40 ms, float 6.60).

The CPU times of that second run differed by up to 40% (3D double step
112 ms, 2D double 14.6 ms). The CPU step is noisy, and in several runs it
came out above the sum of its parts (`ghosts_update` + stencil). Read the
CPU ratios as good to about 1.5 times either way.

Reading it:

- **3D:** the GPU takes a step 8.6 times faster than 16 threads in double,
  7.3 in float. Every phase is faster on it.
  - The ghosts are the largest share on both, about half the step. The
    copies make most of the fill: 15.5 M copied ghost cells against 3.2 M
    from jumps.
  - The stencil and prolongation reach about 130 GB/s on the device.
- **The blocks must live on the device.** Moving every block once (63 ms
  in double, PCIe, pageable memory) costs 7 GPU steps. A ghost update on
  the host, the blocks moved there and back, would cost 126 ms a step in
  moves alone: more than the CPU's whole step (80 ms).
- **2D, on these small trees (42 MB), the GPU is launch-bound.**
  - Restriction and the ghosts launch once per level and kind: about
    0.2 ms a launch whatever its size, so they lose to the CPU.
  - Only the stencil, one launch, wins.
  - In float the CPU takes the step.
  - Larger trees or larger blocks change that: see the block sizes below.
- **`gpu_cells ()` was wasting the GPU.** gcc's nvptx back end makes an
  OpenMP thread a warp, and its 32 lanes work only in a simd loop. Without
  simd, one lane in 32 did the work.
  - With `simd` added to the four collapsed loops, it was 2.8 times
    slower still.
  - `gpu_cells ()` is now one flat loop over block x cell with `simd`.
    Each axis of the region visited (m + 2 halo cells) is padded to a
    power of 2, so that the flat index splits by shifts. The lanes that
    fall in the padding do nothing, and at halo 0 there is none, as m is
    a power of 2.
  - The signature is the same, and so are the cells visited and, on the
    host, their order (that of `block_cells`). Its locals are prefixed
    `bamr_`, so that none hides a name the loop's body uses (the
    benchmark's restriction has an `h` of its own).
  - It is 13 times the old loop on the 3D stencil and 20 times on
    prolongation, in double:

    | 3D, m 16, ms, double (float) | old `gpu_cells` | new | 16 CPU threads |
    |---|--:|--:|--:|
    | stencil | 47.3 (12.4) | 3.52 (1.79) | 36-39 (19-21) |
    | prolongation | 47.7 (14.0) | 2.38 (1.39) | 25-27 (15-16) |
    | restriction | 7.55 (2.57) | 1.42 (1.06) | 11-13 (5.6-5.8) |
    | step | 59.3 (18.3) | 9.40 (6.60) | 80-112 (45-48) |

    The CPU column is the range of two runs.
  - Checked:

    | test | result |
    |---|---|
    | `GPU=1 sh tests/run.sh`, bit for bit, with every halo and both layouts in `test_blocks.c` | 42 of 42 |
    | `sh tests/run.sh`, the whole host suite | 109 of 109 |
    | break: the padding test off by one (`>=` to `>`) | `blocks` fails, 6 of 6 runs |
    | break: cj off by one | `blocks` fails, 4 of 4 (2D, 3D; 1D has no cj) |
    | break: the last block left out | `blocks` fails, 6 of 6 |
    | `tests/run2.sh`, device against host to a tolerance | every check, 2D and 3D, double and float |
    | `tests/run2.sh`, GPU times against the benchmark's own flat loop | the same to 1% (3D step 9.40 ms double, 6.60 float) |

#### 2D, block sizes

    DIMS=2 ARGS2="5 9 4;4 8 5;3 7 6;2 6 7;1 5 8;6 10 4;5 9 5;4 8 6;3 7 7;2 6 8" sh tests/run2.sh

The same shell at two finest cells, 1/8192 and 1/16384, with blocks of
m = 16 to 256.
- **The trees:** base and depth shift with m, so that m 2^depth stays
  fixed. The tree column is base-depth-blockdepth. With m = 128 and 256,
  2:1 balance from the shell lifts the coarsest leaves to level 3.
- **Mcells:** the cells a step writes.
- **Times** in ms. GPU ghosts and GPU stencil are those phases alone;
  upload moves every block.
- **The GPU columns** ran on the benchmark's own flat loop, the same as
  the new `gpu_cells ()`.

| finest cell | m | tree | Real | Mcells | CPU step | GPU step | CPU / GPU | GPU, old `gpu_cells` | GPU ghosts | GPU stencil | upload |
|--:|--:|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| 1/8192 | 16 | 5-9-4 | double | 3.02 | 3.96 | 3.73 | 1.1 | 6.91 | 2.27 | 0.44 | 4.55 |
| 1/8192 | 16 | 5-9-4 | float | 3.02 | 1.03 | 3.55 | 0.3 | 4.59 | 2.22 | 0.35 | 2.26 |
| 1/8192 | 32 | 4-8-5 | double | 4.85 | 6.97 | 3.88 | 1.8 | 9.99 | 2.25 | 0.52 | 7.20 |
| 1/8192 | 32 | 4-8-5 | float | 4.85 | 1.90 | 3.66 | 0.5 | 5.56 | 2.22 | 0.41 | 3.45 |
| 1/8192 | 64 | 3-7-6 | double | 8.56 | 13.87 | 4.20 | 3.3 | 16.00 | 2.22 | 0.72 | 12.46 |
| 1/8192 | 64 | 3-7-6 | float | 8.56 | 4.93 | 3.85 | 1.3 | 7.41 | 2.21 | 0.52 | 5.85 |
| 1/8192 | 128 | 2-6-7 | double | 15.40 | 35.42 | 4.55 | 7.8 | 26.35 | 1.78 | 1.47 | 22.28 |
| 1/8192 | 128 | 2-6-7 | float | 15.40 | 14.04 | 3.53 | 4.0 | 10.06 | 1.75 | 0.74 | 11.62 |
| 1/8192 | 256 | 1-5-8 | double | 22.57 | 34.34 | 4.98 | 6.9 | 37.59 | 1.29 | 2.47 | 34.57 |
| 1/8192 | 256 | 1-5-8 | float | 22.57 | 21.42 | 3.43 | 6.2 | 13.55 | 1.28 | 1.28 | 18.50 |
| 1/16384 | 16 | 6-10-4 | double | 6.96 | 18.89 | 4.30 | 4.4 | 11.55 | 2.52 | 0.66 | 11.29 |
| 1/16384 | 16 | 6-10-4 | float | 6.96 | 7.22 | 3.97 | 1.8 | 6.19 | 2.45 | 0.46 | 5.59 |
| 1/16384 | 32 | 5-9-5 | double | 10.15 | 18.47 | 4.56 | 4.1 | 17.25 | 2.47 | 0.80 | 15.47 |
| 1/16384 | 32 | 5-9-5 | float | 10.15 | 16.23 | 4.13 | 3.9 | 7.92 | 2.42 | 0.57 | 9.35 |
| 1/16384 | 64 | 4-8-6 | double | 17.71 | 27.58 | 5.23 | 5.3 | 29.52 | 2.42 | 1.20 | 25.69 |
| 1/16384 | 64 | 4-8-6 | float | 17.71 | 89.01 | 4.54 | 19.6 | 11.71 | 2.42 | 0.80 | 13.91 |
| 1/16384 | 128 | 3-7-7 | double | 32.70 | 64.76 | 7.46 | 8.7 | 53.68 | 2.43 | 2.79 | 51.76 |
| 1/16384 | 128 | 3-7-7 | float | 32.70 | 30.33 | 5.25 | 5.8 | 19.52 | 2.37 | 1.26 | 25.05 |
| 1/16384 | 256 | 2-6-8 | double | 60.14 | 110.97 | 10.87 | 10.2 | 97.28 | 1.98 | 5.95 | 96.58 |
| 1/16384 | 256 | 2-6-8 | float | 60.14 | 53.59 | 6.79 | 7.9 | 34.52 | 1.94 | 2.97 | 47.29 |

In short, the GPU's lead on a step (CPU / GPU; below 1 the CPU wins):

| m | 1/8192 double | 1/8192 float | 1/16384 double | 1/16384 float |
|--:|--:|--:|--:|--:|
| 16 | 1.1 | 0.3 | 4.4 | 1.8 |
| 32 | 1.8 | 0.5 | 4.1 | 3.9 |
| 64 | 3.3 | 1.3 | 5.3 | 3.7 * |
| 128 | 7.8 | 4.0 | 8.7 | 5.8 |
| 256 | 6.9 | 6.2 | 10.2 | 7.9 |

\* The CPU step of that run is off (89 ms, see below): 3.7 is from its
phases instead, restriction + ghosts + stencil (2.1 + 2.7 + 12.0 = 16.9
ms) against the GPU's step (4.54 ms).

- **The GPU has a floor of about 3.5 ms a step,** whatever the size.
  - Restriction (about 1 ms) and the ghosts (about 2 ms) launch once per
    level and kind, at about 0.2 ms a launch.
  - Up to about 5 M cells a step, that floor is the step.
- **Per cell, the CPU stays about flat; the GPU gets cheaper with size.**
  - CPU: 1.2-1.6 ns a cell on the double stencil for every m.
  - GPU: from 0.27 ns at m = 16 (1/8192) to 0.05-0.14 ns from m = 128,
    as the launches are spread over more cells.
- **So the GPU's lead grows with m and with resolution:**
  - at 1/8192, from 1.1 (double) and 0.3 (float) at m = 16 to 7-8 and
    4-6 at m = 128-256;
  - at 1/16384, from 4.4 and 1.8 at m = 16 to 10 and 8 at m = 256;
  - the CPU wins only in float at m = 16 and 32 on the coarser shell.
- **The ghosts on the GPU cost about 2 ms whatever m is.** The CPU's
  shrink as m grows (from 7 ms at m = 16 to 2 at m = 256 at 1/16384,
  double), as the ghost cells' share falls as 4 g / m.
- **One CPU point is off:** 4-8-6 float took 89 ms a step, against 28 in
  double, while restriction, ghosts and stencil add up to 16.9 ms; its
  `ghosts_update` (19.3 ms, against 4.9 for its parts) is off too.
  Probably a disturbed run; it was not repeated.
- **Old `gpu_cells ()` against new:**
  - 1.3 to 2.7 times the step at m = 16;
  - 3 to 9 times at m = 128-256, where the stencil and prolongation
    dominate.

Next, on the GPU: the launches.
- **About 0.2 ms a target region** is far more than a CUDA kernel launch
  costs (microseconds). It is likely libgomp's own work at each region:
  its mappings looked up, its arguments copied, a synchronisation. Not
  measured apart yet: an empty target region, timed, would say.
- **The ghosts and restriction launch once per level and kind,** so the
  2D GPU step has a floor of about 3.5 ms.
- **Fewer launches would lower it:**
  - levels fused, as far as the order allows (a level jump reads the
    coarser level's ghosts, filled first);
  - the copies and prolongations of a level in one kernel;
  - `nowait`, with dependences.
- **That matters most in 2D and for small trees.**

### Benchmark: the level indices

    sh tests/bench_memindex.sh [RAW]

It builds `tests/bench_memindex.c` with gcc, in 2D and 3D, and runs each case
(dense, shell, river, random) once per implementation (memindex,
voxel-memindex, nano-vdb, nano-vdb compact), each in its own process. It
prints a table, and checks that every implementation returns the same values
as memindex (checksum). nano-vdb's dense nodes skip the scattered 3D cases:
about 300 KB per isolated point. Results and their reading are in
`tests/bench_memindex.md`, summed up below in
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex).
The sparse cases need about 1.5 GB for memindex. Cache and dTLB misses are
printed only if `kernel.perf_event_paranoid` is 2 or less.

    sh tests/bench_balanced.sh [RAW]

It does the same on 2:1 balanced trees, built as `tree.c` builds them. It
times a depth-first Morton traversal, a pass over a level list, and the
neighbour box of every leaf. Results are in `tests/bench_balanced.md`.


## Level index: memindex, voxel-memindex, nano-vdb

Each level of the tree has an index: a sparse (2^level)^d array of node
pointers. `level-index.h` picks the implementation when compiling, and
`tree.c` and `tree.h` use only its interface:

| `-DBAMR_INDEX=` | files | structure |
|---|---|---|
| `1` (default) | `memindex.[ch]` | per-level virtual arrays, as Basilisk: `d` loads, the system commits pages as they are written |
| `2` | `voxel-memindex.[ch]` | a tree of 64-way nodes (4^3, 8^2, 2^6), 512 bytes each, from a pool |
| `3` | `nano-vdb.[ch]` | a VDB tree: root hash, upper 32^d, lower 16^d, leaf 8^d nodes, each with a mask of its occupied slots; `-DNVDB_COMPACT` stores only the occupied children, found by rank (`popcount`) |

The interface, for a `LevelIndex * m` (IJK: `i`, `j`, `k` as the dimension
has them; `len`: `1 << level`):

| | |
|---|---|
| `lidx_get (m, IJK)` | the pointer, or `NULL`, out of range too: one walk |
| `lidx_allocated (m, IJK)` | `lidx_get () != NULL` |
| `lidx_assign (m, IJK, len, p)` | `p`, not `NULL`, at IJK |
| `lidx_free (m, IJK, len, p)` | IJK empty again; empty nodes go back at once |
| `lidx_destroy (m)` | everything freed; `m` empty and reusable |
| `lidx_get_cursor (m, &c, IJK)` | `lidx_get ()` through a `LevelCursor c` |

A cursor makes nearby lookups cheaper, and is reset (`{0}`) after any
`lidx_free ()` or `lidx_destroy ()` of its index:
- memindex: nothing to keep, `mem_get ()`;
- voxel-memindex: the last bottom node (an index in the same 64 is one
  load);
- nano-vdb: the last leaf, lower and upper node, as VDB's `ValueAccessor`
  (the same leaf is one load, the same lower node two).

Each implementation also keeps its own names (`mem_`, `vmem_`, `nvdb_`),
so that the three can be used side by side, as `tests/test_index.c` and the
benchmarks do. `BAMR_INDEX_NAME` names the one chosen.

**nano-vdb, in short.** An index's bits are cut as VDB cuts a coordinate:
bits 2..0 pick a leaf slot, 6..3 a lower-node slot, 11..7 an upper-node slot,
the first index most significant. While `len <= 4096` there is a single
upper node; above, the root is a hash table of upper nodes keyed by the
indices `>> 12` (open addressing, deletion by backward shift). Dense nodes
are tables of every slot, as OpenVDB's internal nodes: one load per level,
but 32^3 slots (256 KiB) for a 3D upper node. Compact nodes hold a mask, a
running count of set bits per 64-bit word, and only the occupied children;
slot s is at `prefix [s/64] + popcount (mask [s/64] & ((1 << s%64) - 1))`.
Build compact nodes with `-mpopcnt` or `-march=native`: otherwise gcc calls
a library function for each popcount.

# Future plans

## Block allocation

Ideas to make allocating blocks, and sweeping over them, faster. They are
ranked by the gain expected for libb-amr as it is now. The gains are
estimates, not measurements.

How allocation works now:
- Every block comes from a pool (`pool.c`) of 1 MiB `mmap` chunks.
- A large block that is freed gives its pages back to the system at once
  (`madvise`).
- Each node holds a two-level table of its blocks, `s[id/8][id%8]`.

### 1. Stop handing pages back on every free

**Done.** `pool.c` keeps freed objects on a hot list, and
`pool_alloc_uninit ()` skips the zeroing; `pool_trim ()` (or
`grid_scalar_trim ()`) gives pages back after a coarsening. Measured by
`tests/test_pool.c`, reusing 62 KiB blocks: about 26 µs a block before (15
page faults each), 2.4 µs zeroed and 1.3 µs uninitialised now.

Refinement fills children in full, so their blocks use
`pool_alloc_uninit ()`.

**Lazy release, to be decided** (`_BAMR_MADV_DONTNEED_RECOMMENDED_` in
`pool.h`):
- **1, the default:** objects over the hot limit go cold with
  `MADV_DONTNEED`.
- **0:** they go to a "muzzy" list with `MADV_FREE`, all pages but the first
  (which holds the list's link). The system takes those pages only under
  memory pressure, so reuse costs no fault while memory is plentiful.
  `pool_trim ()` still gives everything back for real. The drawback: muzzy
  pages count in the RSS until the system takes them.

Measured by `tests/test_pool.c` (62 KiB blocks, over the hot limit, then
reused uninitialised):

| | faults on reuse | free | alloc |
|---|---|---|---|
| `MADV_DONTNEED` | 15 per block | 8.6 µs | 27.5 µs |
| `MADV_FREE` | 0 | 6.3 µs | 3.9 µs |

That is without memory pressure. Under pressure (a cgroup with a small
`memory.max`), `MADV_FREE` should fault like `MADV_DONTNEED`; this has not
been measured yet.

The reasoning, as planned:

`pool_free` gives a large block's pages back to the system at once. The next
`pool_alloc` of that block then takes a page fault for every 4 KiB page, and
the kernel zeroes each one. A 3D block with m = 16 is 20³ doubles, 62 KiB:
that is 16 faults, very roughly 10–20 µs. A stencil over that block's 4096
cells takes about 4 µs. A grid that refines and coarsens every few steps could
spend more time in page faults than computing.

The fix:
- **A hot free list.** Freed blocks are kept as they are, up to a limit, such
  as the most blocks freed in one regrid. Only blocks beyond the limit go back
  to the system.
- **`pool_alloc_uninit ()`** for blocks that are overwritten straight away, as
  when refinement fills children from their parent. These skip the zeroing
  too.

### 2. Allocate the children of a node as one group

**Done.** `tree_refine ()` allocates a family at once: one object of
TREE_NCHILDREN nodes, and for each scalar one object of TREE_NCHILDREN
blocks, which `tree_coarsen ()` frees together. Children are filled by linear
prolongation, which writes every cell, so their blocks are not zeroed
(`pool_alloc_uninit ()`).

The reasoning, as planned:

`add_children` always creates all 2^d children at once, and in a 2:1-balanced
tree a family is also coarsened at once. So allocate one object of 2^d blocks
per scalar, instead of 2^d separate ones:
- one pool operation instead of 8 in 3D;
- siblings are adjacent in memory, and they are also adjacent in Morton
  order, so a traversal walks memory nearly sequentially;
- stencils that cross between siblings stay within the same pages.

This is the "quadrant family" of p4est and similar codes. `grid-scalar` would
allocate per family, and each child's pointer would be an offset into the
family's object.

### 3. Ghost cells are a large share of memory in 3D

**Tried, then reverted: blocks are stored with their ghosts again** (see
[Ghost cells](#ghost-cells)).

For a while, blocks stored m^d cells, and `block_gather ()` filled a
scratch block with ghosts before every stencil. It halved memory in 3D. It
was reverted because:
- each sweep repeated the gather, and a coarse-neighbour ghost had to be
  rebuilt every time it was read;
- interpolating at level jumps, which second-order stencils need, made that
  repeated work larger;
- custom prolongations and derived scalars fit a single fill per step much
  better.

Storage per block, relative to the interior alone:

| | m = 16, g = 2 | m = 32, g = 2 |
|---|---|---|
| 1D | 1.25× | 1.13× |
| 2D | 1.56× | 1.27× |
| 3D | **1.95×** | 1.42× |

Ways to reduce it, if memory becomes the limit:
- **Larger m**: cheaper ghosts, but coarser refinement.
- **Fewer ghost layers** for scalars whose stencils need only one: a
  per-scalar g would need a block pool per size.
- **Face ghosts only**: a 5- or 7-point stencil never reads edge or corner
  ghosts.

### 4. Huge pages

**Gain: moderate. Effort: small.**

A sweep over many blocks of several scalars touches far more 4 KiB pages than
the TLB can hold: 300 MB of blocks is about 75,000 pages. With 2 MiB pages it
is about 150.
- Make pool chunks 2 MiB-aligned multiples of 2 MiB.
- Call `madvise (MADV_HUGEPAGE)` on them.

Bandwidth-bound codes often gain 5–20%.

This conflicts with 1: giving back part of a huge page splits it. So for
pools of huge pages, rely on the hot free list, and give memory back only in
an explicit `pool_trim ()`.

The reverse holds for `memindex`: huge pages only hurt its sparse arrays. With
transparent huge pages set to `[always]`, one pointer written can commit 2 MiB
instead of 4 KiB, up to 512 times the memory. `mymap ()` therefore calls
`madvise (MADV_NOHUGEPAGE)` on every array it maps. Under `[always]`,
khugepaged can also collapse a pool chunk whose cold objects were given back.
If that shows up, give pool chunks `MADV_NOHUGEPAGE` too, or adopt the
huge-page pools described above.

### 5. Contiguous storage in Morton order

**Gain: large for a GPU and for prefetching. Effort: medium.**

After each regrid:
- number the leaves in Morton order;
- keep each scalar's blocks in one slab, `u[rank * cells + c]`.

A traversal then reads memory sequentially, and the hardware prefetcher
follows blocks across their boundaries. The slab is also the contiguous device
array that `gpu_cells` (`block-loops.h`) wants. It needs no gathering through
`block_table`, and no table of pointers to map.

The cost is copying all blocks at each regrid: one `memcpy` per block, cheap
next to many computation steps. Combined with 2, it becomes "families in
Morton order".

### 6. A flat scalar table per node

**Gain: small, and simpler code. Effort: small.**

`node->s[id/8][id%8]` follows two pointers. It also needs a group pool, and
the group bookkeeping in `grid_scalar_check`.

A flat table of `NSCALAR_MAX` pointers is 512 bytes per node, under 1% of a
single 3D block. It removes one pointer per access and the group pool, and
makes `grid_scalar_node_free` simpler. The groups only pay off when blocks are
tiny, and `tree ()` already requires m ≥ 16.

### 7. NUMA and threads

**Matters on machines with more than one socket.**

Fresh `mmap` pages are placed in the memory of the first thread that writes
them. If blocks are first written in the same `foreach_leaf_parallel` pattern
later used to compute on them, each thread's blocks end up in its own socket's
memory. The uninitialised allocation of 1 helps, because the allocator itself
then never writes the memory.

`pool.c` is not thread-safe. If refinement ever runs in parallel, it needs
per-thread free lists, or a lock around `pool_alloc` and `pool_free`.

### 8. Further options

- **Mixed precision.** Done in part: `real.h` gives `Real`, a double, or a
  float with `-DBAMR_SINGLE`, and blocks store `Real`. Geometry stays double.
  What remains:
  - **Store in float, compute in double.** Blocks hold float, and kernels
    convert each value to double, compute, and convert back. Stencils are
    limited by bandwidth, so this is nearly as fast as float, with double
    arithmetic. It needs a computation type beside `Real`.
  - **Precision per field.** Pressure in double and a tracer in float, in one
    run. It needs a size per scalar, a block pool per size, and a typed
    accessor.
- **Reserving ahead.** Before a regrid the tree knows how many blocks it will
  create, so `pool_reserve (n)` could map them in one call.
- **To be considered: blocks per node for the halo (MPI).** Under MPI a
  family has blocks iff one of its members is active
  (`BAMR_NODE_IS_INACTIVE`, the skeleton of tree-mpi.h), so a halo family
  with one active member carries 2^D blocks: in 3D, 7 of 20^3 Reals (64 KB
  each in double) for nothing. Measured on `test_mpi`'s small trees the
  nodes with blocks are 0.33 of the tree against 0.37 kept (2D, 16
  ranks), 0.43 against 0.45 (3D): little saved there, as halo families
  are mostly mixed. Blocks per node in mixed families would save the rest,
  at the cost of the one-object family. Kept per family for now: the
  author plans to cut the partition by edge cut (fewer halo nodes), which
  changes what the halo costs; to be weighed after that.

### Order of work

| Order | Change | Effort | Expected gain |
|---|---|---|---|
| 1 | hot free list + uninitialised allocation | small | large, if the grid adapts often |
| 2 | flat scalar table | small | small, simpler code |
| 3 | huge pages, aligned chunks | small | moderate |
| 4 | children allocated as a group | medium | moderate (locality, fewer allocations) |
| 5 | Morton-ordered slabs, renumbered at each regrid | medium | large for a GPU and for prefetching |
| 6 | ~~blocks without stored ghosts~~ | large | tried and reverted: see 3 |

**Measure before changing anything:** time a stencil sweep, and a refine and
coarsen cycle. Use `grid_scalar_check ()` and `pool_used ()` to make sure no
change leaks blocks.

## Tree structure: voxel-memindex

`voxel-memindex.[ch]` (F, compared under
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex))
beat memindex on memory and build time for every level that is not full. Next
steps, in order:

1. **Wire it into `tree.c` behind a switch.** Done, with nano-vdb too:
   `-DBAMR_INDEX=2` or `3` (see
   [Level index](#level-index-memindex-voxel-memindex-nano-vdb)).
   `BAMR_NODE_ENV ()` makes a single lookup. The tests pass on every index.
   What remains: time `adapt` and `ghosts_update ()` on a sparse 3D case
   with each index. That is the measure that decides, rather than isolated
   lookups.
2. **Use cursors in the hot loops.** `block-ghosts` gathers neighbours and
   `adapt` checks 2:1 balance: both visit nearby indices, so one cursor per
   level makes most lookups a single load. Cursors must be reset after a
   `vmem_free ()`, so reset them at the end of each regrid pass. Each OpenMP
   thread needs its own cursors.
   - For the neighbour box, add a box cursor that holds the up to 2^d bottom
     nodes a leaf's box can touch. On balanced trees the box is the one place
     where memindex is still faster, by up to 36% (see On 2:1 balanced
     trees).
3. **Compact nodes for sparse levels.** A node holds a 64-bit mask and only
   its occupied children, found by rank (`popcount` from libutils/bits.h).
   - It would cut the 2–4 kB per isolated point by up to 8 times.
   - Children sit in arrays of varying sizes, from a few size-classed pools.
   - Inserting moves the children of one node, so batch it per regrid.
   - Possibly mixed: full nodes near the bottom, compact ones above.
4. **Memory after freeing.** Freed nodes stay hot in the pool. Call
   `pool_trim ()` after a coarsening pass, or give voxel its own hot limit
   (`pool_set_hot_max ()`), so that memory returns after a large coarsening.
5. **Measure the misses.** Rerun the benchmark with
   `kernel.perf_event_paranoid=1`, to get LLC and dTLB misses per lookup.
   Rerun it on a cluster node with transparent huge pages `[always]`, and with
   memindex without `MADV_NOHUGEPAGE`, to quantify the blow-up seen elsewhere.
6. **Later, for MPI or GPU:** the nodes already give Morton order, so a
   linear octree (C) can be built in one traversal of the leaves.

## Ghost operators and adaptation

Deferred; none started.

1. **Periodic boundaries.** Kept for later when the per-face boundaries
   were done (`boundary.left (op)`, ..., `boundary.all (op)`). A periodic
   pair of faces needs more than an operator: `block_ghosts ()` would fill
   the regions beyond such a face from the node on the other side of the
   domain, copied or prolonged like any neighbour; `neighbors_of ()` and
   the 2:1 balance would wrap around too; and the walls of those faces
   would be skipped.
   - **Walls in the node flags.** Whether a block is on a wall is now
     read from its index: `block_upper_walls ()` (`block-loops.h`), i ==
     2^level - 1 along each axis, and `q [a] < 0 || q [a] >= nb` for the
     regions outside (`block-ghosts-impl.h`). A periodic axis has no
     walls, so the index alone no longer tells. Instead, `tree.c` would
     set wall bits in `node->flags` when it makes a node
     (`allocate_root ()`, `add_children ()`), from the index and the
     domain's periodic axes; coarsening never moves a node, so the bits
     never go stale.
   - **As planned:** six bits, the lower and the upper wall of each axis,
     in axis order, so that `(flags >> shift) & 7` is today's mask (bit
     d for axis d, as `walls & c` with a child number wants); a mask
     (`BAMR_NODE_UPPER_WALLS`), not a seventh bit, for "any upper wall";
     `block_upper_walls ()` kept as what callers use, reading the flags;
     code with an index and no node (GPU kernels from `block_table ()`)
     may keep the compare, with the periodic axes. A test checks the
     bits against the index on every node, and `tests/bit/compare.sh`
     stays bit for bit.
   - **Not for speed:** the compare is done once per block, before work
     on m^D cells, and every caller loads the node anyway.
2. **A refinement operator, `@refine`, apart from the prolongation.**
   Basilisk has both. Now one operator fills the children of a refined
   block and the ghosts at a level jump, so the two agree: a fine cell
   reads the same values beside a coarse leaf before and after it is
   refined. The defaults (`prolong_linear`, `prolong_minmod`) keep the
   mean of the children, which is the main reason Basilisk needs a
   separate one (its bilinear default does not).
   - **Where it would differ:** new cells must keep the integral, while a
     ghost need only be smooth and accurate. Refinement runs once per new
     cell, so it can be expensive: a reconstruction, a solenoidal fill of
     face fields. And it can reach further, as the parent's ghosts are all
     valid.
   - **The case that decides it for the flood model:** if `zb` is read
     again from the DEM when a cell is refined, the children's mean `zb`
     is not the parent's. A lake-at-rest prolongation of `h` then gains or
     loses water at each refinement; a refinement operator would keep the
     surface flat and correct `h` to the parent's mean depth.
   - **As planned:** an optional fourth kind, `ghost_ops (h, refine
     (op, zb))`, the same ids, dependencies and order as the others;
     `block_prolong ()` takes it when bound, the prolongation otherwise;
     the docs to say that two operators lose the agreement above.
3. **Operators of a scalar freed with no tree.** The hooks of
   `grid-scalar.c` reset a scalar's operators when it is made and freed,
   but only while a tree exists. A scalar freed without one keeps its
   bindings, and the next scalar given its id inherits them. The
   dependencies are rebuilt over the scalars in use only, so a free id's
   bindings count for nothing meanwhile. The check that no scalar depends
   on one being freed (`ghost_ops_freeing ()`) is in the same hook, so it
   too runs only while a tree exists. The fix: reset and check in the free
   path whether or not a tree exists.
4. **Prolongation's bookkeeping per cell.** `prolong_points` computes the
   coarse cell, its half and two `block_at ()` for each fine cell. Its
   range is checked once per call now; the indices could also be stepped
   along a row instead of recomputed. Only worth it if 2D prolongation
   shows up in a profile: after the range check, the wavelet in 2D is
   faster than the old in-block one for every operator.
5. **Multigrid:** a first version is done (`liblinear-algebra`): the
   loops over a level (`foreach_level (l)`, `foreach_level_or_leaf (l)`,
   and `foreach_level_or_leaf_colour (l, colour)` over one colour,
   `BAMR_RED` or `BAMR_BLACK`, or every cell in order, `BAMR_BOTH`) and
   `ghosts_update_level (s, l)`, called explicitly. Still open: lazy
   ghosts per level, `coarse ()` inside `foreach ()`: see `multigrid.md`.
6. **The dependency closure, in near-linear time** (discussed, to do
   later). `ghost_dependents_rebuild ()` (block-ghosts-impl.h) makes the
   transitive closure of "who depends on s" by passes of `foreach_scalar
   (t) foreach_scalar (s)` until nothing changes: O(n² W) a pass (W =
   ⌈n / 64⌉ words), O(n³ W) at worst, with n the scalars in use. Small at
   `NSCALAR_MAX` 64 (W = 1), and it runs only when a scalar is made,
   freed, or its operators set, not per step; the closure itself stays,
   since it makes marking a scalar stale one OR (block-ghosts.h).
   1. **One pass, O((n + E) W).** The dependencies are a DAG (cycles are
      refused), so a depth-first walk in topological order gives
      closure (t) = the union, over t's direct dependents s, of {s} and
      closure (s): each scalar visited once, each edge (E of them, at
      most GHOST_DEPS_MAX a scalar) read once. The walk finds a cycle
      too (a scalar met again while on its stack), which replaces
      `ghost_dependents_cycle ()` and the rebuild, check, undo of
      `ghost_ops_set ()`.
   2. **No rebuild when nothing can change, O(1).** A scalar with the
      default operators has no dependency and none on it: making or
      freeing it changes no closure but its own row and stale bit. A
      test of edges either way (its uses, its bit in the closure) skips
      the rebuild: almost every `grid_scalar_new ()` and
      `grid_scalar_free ()`, the reset at new included.
   3. **Incremental, not worth it at n <= 64.** Adding an edge "s reads
      d" can be O(n W): a cycle iff s is d or s is in d's closure, one
      bit; else every x that is d or holds d gets {s} and closure (s).
      Removing edges is the hard direction (a closure does not know
      which paths an edge carried): recompute with 1.

   Tests for it: those that bind coupled operators (`test_ghosts`,
   `ops_coupled` in libc-pp, libvof's `prolong_plane` and its
   dependencies) must give the same closures and refuse the same cycles,
   and a timing of `grid_scalar_new ()` and `free ()` before and after.

## Faces, fluxes and vertices

Deferred; none started. What is done is above ("Faces and fluxes",
"Vertices"). These would make it faster or more exact.

1. **A cheaper face ghost update for fluxes.** Reading a stale face vector
   makes every ghost layer of every component in every direction, and
   restricts every parent. A flux divergence reads only `F.x[]` and
   `F.x[1]`: each block's faces and its upper face m, made to agree at
   level jumps. Measured in liblinear-algebra's residual (2D, 268 leaves):
   the face ghost update is 0.95 ms of 1.48, against 0.39 for the whole
   cell-centred residual it replaced.
   - **A partial update:** only the shared faces and the reflux: face m
     from the node above, face 0 from a finer region below, and the
     parents restricted only where a coarse leaf touches a finer region.
   - **How a loop asks for it:** cwrap knows what a loop reads. If a face
     vector is read only at offset 0 or 1 along its own axis, it could
     call the partial update; any other read, the full one.
   - **Or no ghosts at all for a flux:** keep the shared faces in the
     owner's block only, and read face m from the neighbour. That is a
     larger change to how faces are stored.
2. **A flux-form smoother at level jumps (multigrid).** The residual is in
   flux form, but the smoothers relax the cell stencil, so at a jump the
   two operators differ. That costs cycles: 20 against 14 in 2D, 28
   against 20 in 3D (red-black, liblinear-algebra's benchmark). A smoother
   whose diagonal and neighbours come from the refluxed faces at a jump
   would match the residual. With the partial update above, flux form
   would be back near the cell-centred form's 47 ms (2D) and 602 ms (3D),
   from 131 and 1913.
3. **Solenoidal refinement of face velocities.** `prolong_face` gives a
   new face across a coarse cell's middle the mean of the coarse cell's
   two faces. A divergence-free face velocity can lose that property on
   refinement. Basilisk's `refine_face_solenoidal` fixes the new interior
   faces so every fine cell's divergence is the coarse one's (zero). This
   is needed for a projection method, such as a non-hydrostatic pressure.
   It would be a `@prolongation` bound per face vector with `ghost_ops`.
4. **Second order across the axes for smooth face fields.** Along the
   other axes `prolong_face` is piecewise constant: right for a flux (per
   unit area, summed back exactly), first order for a smooth face field
   such as a velocity. A variant with a limited slope across the other
   axes, keeping the mean of the fine faces equal to the coarse face's,
   would be second order and still conservative.
5. **Boundary conditions of faces and vertices.** The ghosts beyond a wall
   are made as a cell scalar's: a mirror of the cells.
   - **Faces:** the normal component's wall face is the block's own (made
     by `foreach_face ()`), so a flux through a wall is what the loop
     says. But the ghost faces beyond it (index -1, m + 1) are mirrored
     as cells, not as faces, and a stencil reaching them reads the wrong
     place. Proper ones would mirror about the wall face: odd for the
     normal component of a velocity, even for a tangential one.
   - **Vertices:** a vertex on a wall is its owner's; beyond it, mirrored
     about the wall vertex, not as cells.
   - Both would be `@boundary` operators, defaults by location as the
     restriction and prolongation are.
6. **Faces and vertices on levels (multigrid).** `ghosts_update_level (s,
   l)` treats a face or vertex scalar as a cell one. The solver needs only
   `alpha`, made once by `ghosts_update_levels ()`. A flux-form smoother
   (2) would need the faces of the depth-l grid to agree at its jumps:
   the reflux of the depth-l grid.
7. **Diagonal-only level jumps for vertices.** A vertex is the owner's
   where blocks meet along an axis or across a wall, and corners take what
   the cell machinery gives them. A corner that touches a finer region
   only diagonally (no axis neighbour finer) gets the coarse side's
   interpolation; the fine side's own vertex may differ by its own
   computation. The tests cover axis-adjacent jumps and walls; a test
   with a diagonal-only jump would tell whether it needs the owner's copy
   too.

## Restriction

What a ghost update restricts, and restrictions that read beyond a child.
Item 1 is done in part; the rest are not started.

**Done.** `ghosts_update ()` restricts the parents of the coarsest leaf's
level and finer only (`tree_get_depth_min ()` .. `depth - 1`, finest first),
and the restriction loops run over the tree's list of parents
(`tree_caches ()->parents_by_level`) instead of testing every node. Every
parent is still restricted by `ghosts_update_levels ()` and `tree_restrict
()`. Measured in libvof's benchmark (uniform trees, no parent read): the
restrictions were a third of a 3D step; a step went from 203-226 to 144 ns
per cell, 62 to 32 in 2D.

### 1. Restrict only the parents that are read

**Who reads a parent**, in a ghost update without the parents' own ghosts:
only a leaf of its own level beside it, among its 3^d - 1 neighbours,
diagonals included (corner ghosts copy too):
- its ghost copies take the parent's interior rows;
- refluxing takes the parent's face 0 (a face vector), vertices the
  parent's vertices;
- a prolongation into a finer leaf's ghosts reads a coarse leaf's interior
  and first ghost layer, which came from such a parent: the same rule;
- coupled operators read the same coarse cells.

Multigrid (`tree_restrict ()`, `ghosts_update_levels ()`), `adapt_wavelet ()`
(a parent's wavelet) and code that reads coarse levels (`foreach_level ()`,
output) need every parent, and ask for it.

**The recursion.** A parent's restriction reads all of its children, so a
child that is itself a parent must be restricted first: with whole-block
restrictions, the set read is the marked parents and their subtrees. A deep
refined region whose outermost parent touches one coarse leaf would be
restricted whole for that leaf's two ghost layers.

**What 2:1 balance gives.** The coarse leaf reads only a slab of the parent,
`g` coarse cells deep. The fine cells under it are in the children that
touch the coarse leaf, and by 2:1 balance those are leaves: a restriction of
that slab alone never recurses.

| | What | Gain | Cost |
|---|---|---|---|
| A | down to the coarsest leaf's level | all on a uniform tree; the levels below the coarsest leaf otherwise | done |
| B | the parents a leaf of their level is beside, and their subtrees: a list in the tree caches, made again when the tree changes | all but the parents along level jumps and below them | moderate |
| C | only the slab each coarse leaf reads | the least there is; no recursion | every restriction given a region (as `Prolongation` has `lo`, `hi`): mean, face, vertex and coupled ones (`restriction_plane`) |

**The contract.** After `ghosts_update ()`, only the parents read would be
valid (as now below the coarsest leaf's level): to be documented, and
`tests/bit/compare.sh` to hash the leaves and the parents read.

**The order of work:** B, then C only if a profile of a real mesh (a flood
on a coarse floodplain) shows B's subtrees large.

### 2. Restrictions that read beyond a child

**Now.** Every restriction reads its children's own cells only (and the
faces and vertices m a child owns): in `foreach_restriction ()` an access
moves to a child (`foreach_child ()`), never to an offset on the fine level.
So no fine ghost is needed: every restriction runs first, then the ghosts,
coarsest level first. Mean, faces, vertices and libvof's
`restriction_plane` are all compact.

**What would read beyond:** full weighting for vertex-centred multigrid
(libb-amr's vertex scalars), higher-order or wavelet restrictions, a smoothed
coarse field. Such a stencil reads fine cells across the children's
boundaries: the fine level's ghosts, which the order above has not made yet.

**Not by updating every ghost before every level.** Making level l + 1's
ghosts before level l is restricted, for every scalar of the update, would
make all of them pay for the few that need it.

**The reach, known per operator.**
- cwrap: fine-level offsets allowed in `foreach_child ()` (`s[1]`, or a
  `fine (s, i, j)`), checked against a bound at compile time as `coarse ()`
  is against `BAMR_COARSE_REACH`; each operator's largest offset, its reach
  r, recorded and given as `@restriction reach (id)`, as `uses (id)` is.
- libb-amr: `ghosts_update_list ()` splits the scalars it updates into
  compact (r = 0: every operator now) and wide (r > 0). The compact keep the
  schedule as it is.

**What a wide scalar needs** before the parents of level l are restricted:
the fine cells within r of their children, across the children's
boundaries with the nodes of level l + 1 beside them:
1. a node of level l + 1, a leaf or a parent restricted already (finest
   first): its interior, a copy;
2. a wall: the boundary condition;
3. no node of level l + 1, a coarser leaf: the region would come from its
   prolongation, which reads that leaf's first ghost layer, which may come
   from a parent of level l restricted in this very step. A cycle.

**Three ways to give it:**

| | How | Cost | Case 3 |
|---|---|---|---|
| a | interleaved partial ghosts: before level l, a pass over level l + 1 for the wide scalars only, r layers deep, copies and walls, no prolongation | one pass per level, of the wide scalars | left to b |
| b | compact near jumps: the operator falls back to its compact form (the mean) for a child whose stencil would reach a coarser leaf | none | the fallback; 2:1 balance keeps such children few, those beside a jump |
| c | the neighbours read directly: the restriction given the 3^d neighbourhood of fine blocks (the children and their nodes of level l + 1), whose interiors it reads; no ghost pass | proportional to use; every wide operator does the lookups | absent nodes reported, the fallback of b |

**The choice: a with b.** The ghost machinery stays the one place that knows
the neighbours; the pass copies only, r layers, so it is cheap and never
prolongs, which removes the cycle; case 3 is the compact fallback. c only if
wide restrictions become common and the extra passes show in a profile.

**And:**
- **The schedule.** Wide scalars are restricted level by level, each after
  its partial pass: a sequence of passes, a barrier between them; the
  compact scalars of the same update all at once, beside the first. Within a
  level the blocks stay independent (OpenMP; a GPU kernel per level).
- **With item 1.** A wide restriction too is needed only for the parents
  read; the partial pass of level l + 1, only for the families whose parent
  is restricted.
- **Coarsening.** `tree_coarsen ()` restricts one family (`block_restrict
  ()`): a wide operator there needs the family's surroundings, the same
  copy-only gather for one family.
- **Demand.** None yet: cell-centred finite volumes and their multigrid
  restrict by the mean. Full weighting is the first likely user, for a
  vertex-centred multigrid.

### 3. A test of face m beyond its owner

A parent's face m (and a vertex at m) is restricted from its children's
face m, in their ghost layer. Where a child does not own it (a node of its
level is above it), that face is a copy left from the previous update when
the restriction runs: stale. It reaches only faces no leaf reads, by 2:1
balance and the ownership of faces (a parent's face m is read by a coarse
leaf above it only when nothing of the children's level is above them, and
then they made it). `test_faces.c` shows the fluxes balance at every jump;
no test poisons every face m beyond its owner before an update and checks
that no leaf's face changes. One should, as the argument is the only thing
that holds it.

## MPI

Done: the partition (by rank 0, broadcast), the neighbour bit maps, the
plans, pruning that frees the subtrees, the exchange in
`ghosts_update_list ()`, loops over owned leaves with reductions
(milestone 1), and `tree ()` distributing at once, each process finding
its halo by queries to the owners and making its plan from what it holds
(milestone 2), and repartition with the blocks moved (milestone 3) ("MPI:
partition, plans, pruning" above). The distributed tree (plan B,
milestones 4-5) next:

1. **Migration (milestone 3), done:** `tree_mpi_repartition ()`.
2. **Adapting (milestone 4), done:** `adapt_wavelet ()` on a distributed
   tree.
3. **To be considered: a partition by edge cut** (the author's plan). The
   cut along the Morton curve balances the load but not the surface; a
   cut that also minimises the edges between processes (graph
   partitioning on the leaves' adjacency) would shrink the halo, the
   messages and the blocks halo families carry (see "Block allocation",
   "Further options"). The skeleton and the plans do not depend on the
   curve, but the first keys and `tree_partition_owner ()` do: a cut off
   the curve needs another way to find owners.
4. **Multigrid and the rest (milestone 5).** Level exchanges for
   `ghosts_update_level ()`, libvof under MPI, real OpenMPI runs, a first
   scaling measure.
5. **The bit maps in `block_ghosts_cells ()`.** Today it looks up the
   level index on every update, for each of the 3^D - 1 regions. The maps,
   with the neighbours' node pointers, cached per node at each tree change,
   would remove those lookups on one process too, and give the GPU a table.
   The face and vertex rules read from them: `same | fine` at +a keeps face
   m, `fine` at -a refluxes face 0.
6. **Dummy children for the halo** (discussed, not started). The owner of
   a coarse leaf next to a finer region would fill dummy children with its
   prolongation, then send them like any block of the same level. Every
   halo would then be a copy, and a prolongation's dependencies (`zb`,
   `n`, `alpha`) would stay with their owner. Only the GHOSTS-thick strip
   that fine neighbours read needs filling. Keep the dummies in a table
   sorted by `Index`, out of the level index (no API change).
7. **Fewer rounds.** There is one exchange per level, as a fine halo needs
   the coarse leaf's ghosts. A coarse buffer per fine block (Parthenon),
   filled from the coarse neighbours' interiors in the same round as the
   same-level copies, would make it two rounds.
8. **Periodic boundaries** with it: the wall bits of "Ghost operators and
   adaptation", item 1, become the maps' "none" on a non-periodic axis.

## Compile-time block size

Under consideration; not started. The block size `m` is now an argument of
`tree ()` (`blockdepth`, at least 4: m >= 16), and `tree_block_shape ()`
makes the shape again at each call: a call into tree.c and an assert,
about 1.6 ns, for every node of every loop over the tree.

**Measured first: a global shape is not worth it.** A global, set once by
`tree ()` and returned by `tree_block_shape ()`, cut the call to 0.26 ns,
but saved only about 2% of a light `foreach ()` (2D, m 16), and nothing
measurable in `tree_restrict ()`, prolongation, ghost updates or 3D: each
node does hundreds to thousands of cell operations. It was reverted.

**The lever is a constant m:** `-DBAMR_BLOCK_M=M`. Then m, the strides and
every loop bound are constants: index arithmetic folds, loops of known
trip counts unroll and vectorise (gcc 11 vectorises none of the cell loops
at `-O2` now), scratch blocks and GPU tiles get static sizes, and a GPU
kernel its shape. Without the option, a default fit to the target.

What the default depends on:

| factor | pushes m | example |
|---|---|---|
| ghost overhead, ((m+2g)^d - m^d)/m^d | larger | 2D: 56% at m 16, 27% at 32; 3D: 95% at 16, 42% at 32 |
| a block of each scalar a stencil reads, in L1 or L2 | smaller in 3D | 3D, m 16: 64 KiB a scalar, more than a 48 KiB L1 |
| vector width, rows on whole vectors | its multiples | AVX2, AVX-512 |
| GPU: m^d threads, 1024 at most a block; shared memory | 16 in 2D, tiles in 3D | 16^3 = 4096 needs tiles |
| adaptivity: the block is the unit of refinement | smaller | narrow rivers on a coarse floodplain |

The last is the problem's, not the machine's: the default can only be a
guess fit to the hardware, which the user overrides.

The plan, in order:
1. **Prototype and measure.** With `BAMR_BLOCK_M` set,
   `tree_block_shape ()` returns a constant; every function copies it into
   a local, so gcc folds the fields with no other change. Time `foreach
   ()`, prolongation and ghost updates in 2D and 3D, and go on only for
   10% or more.
2. **Defaults in the header,** from the build's own macros (`dimension`,
   `__AVX512F__`, ...), never by probing the machine that compiles: a
   build on a login node runs on compute nodes, and the same source and
   flags must give the same binary.
3. **cwrap's driver picks m when offloading.** gcc compiles the host and
   the device passes apart, each seeing its own macros: a default in the
   header could differ between the two, and the blocks with it. The
   driver sees every flag once, and passes the same `-DBAMR_BLOCK_M` to
   both. A small table of libb-amr's, documented as such: cwrap stays a
   C translator.
4. **A mismatch must not link.** tree.c defines `bamr_block_m_<M>`, which
   the header makes every file reference: a file built with another m is
   an undefined symbol, not wrong indices. The binary records m, so that a
   run reports it.
5. **`tree ()` without `blockdepth`** (or asserting it equals m for a
   while); `tree (5, 4)` in the tests and libpinn follows. `BAMR_BLOCK_M=0`
   could keep today's run-time m, through the same `tree_block_shape ()`,
   while the two are compared. The tests gain a block-size axis, 16 and 32
   at least.
6. **Optional, opt-in: measured tuning.** Time a ghost update and a
   `foreach ()` for each candidate m on the machine itself, and keep the
   fastest, as FFTW keeps its wisdom. Never the silent default.

# Comparison of AMR implementations

How to store the tree itself: its nodes, and the way to find a node or its
neighbours. This is about the tree's structure, not the blocks of field data.

## The starting idea: dense tiles of N levels

One idea is to store the tree densely, N levels at a time:
- a complete tree of depth N is laid out contiguously, as a dense array;
- the same is done again for each subtree whose root is at a depth that is a
  multiple of N.

OpenVDB, NanoVDB and other sparse-voxel trees are built this way.

**The risk is sparseness.** In 3D, a tile of N = 4 levels has
1 + 8 + 64 + 512 + 4096 = 4681 slots. A tile holding just one refined corner
path uses about 33 of them, under 1%.

## Keep the tree's structure separate from its data

A node's structure is a few bytes: flags, an index, a pointer to its scalar
table. A block of field data is about 62 KiB per scalar in 3D with m = 16.
So:
- the structure is roughly 0.01% of memory;
- walking it costs a few nanoseconds per block, against about 4 µs of stencil
  work on the block's 4096 cells.

What follows from that:
- **Sparse structure is cheap.** 37 KiB of mostly empty slots in a tile costs
  about what one block costs.
- **Sparse data storage would be ruinous.** Laying out blocks densely in
  tiles would reserve 62 KiB per empty slot, per scalar.
- **Contiguity matters most for the blocks.** Items 2 and 5 of
  [Block allocation](#block-allocation), sibling groups and Morton-ordered
  slabs, address that. The structure can be organised for fast lookups
  instead.

## Alternatives for the tree's structure

### A. Per-level sparse arrays (what libb-amr has now)

`memindex` gives each level a virtual `[i][j][k]` array, and the operating
system only commits the pages that are actually written.

- **Good:** finding a node or a neighbour is pure index arithmetic, O(1),
  with no pointer chasing. Basilisk uses this scheme.
- **Cost:** sparseness becomes virtual address space, plus one 4 KiB page per
  isolated row of nodes, as a page holds 512 pointers. Refined regions are
  usually clustered, so pages fill well in practice. Each touched page also
  needs its share of a page-table page, which stays after the page is given
  back. It depends on the system: page size, overcommit, and transparent huge
  pages (opted out with `MADV_NOHUGEPAGE`, see Huge pages).
- **Verdict:** this is already a version of the dense-tile idea in which the
  operating system handles the sparseness. Measured against a tree of 64-way
  nodes below: it is cheap only on levels refined fairly evenly.

### B. Dense tiles, with occupancy bitmasks and rank indexing (NanoVDB-style)

Each tile stores:
- a bitmask of its occupied slots: 8 bits per family in 3D, or 4096 bits
  (512 bytes) for a whole tile's bottom level;
- an array of only the occupied nodes, in Morton order.

To find slot `s`, count the set bits of the mask below `s` (`popcount`). With
a small table of running counts, one per 64-bit word, that is O(1). This
removes the sparseness completely, while keeping the tiles contiguous.

libb-amr allows a simplification. In a 2:1 tree whose families are always
created whole, a family's mask is all ones or all zeros. The per-family mask
is then a single bit, and a tile is just a record of which families exist.

- **Good:** contiguous and compact, with O(1) lookup. A tile copies to a GPU
  in one piece.
- **Cost:** refining inserts into a tile's compact array, so tiles should be
  rebuilt in a batch at each regrid, rather than edited one node at a time.

### C. Linear octree: a sorted array of Morton keys (p4est, Dendro)

Store only the leaves, as keys made of the level and the interleaved
coordinates, in sorted order.

- **Good:** fully contiguous, and never sparse. A traversal is a scan.
  Splitting the work across MPI ranks means cutting the array into pieces. It
  goes to a GPU as it is.
- **Cost:** finding a neighbour needs a binary search, O(log n), or a hash on
  the side. Each regrid rebuilds the array, an O(n) merge, which is fine when
  regrids are batched.
- **Verdict:** the standard choice once MPI is needed.

### D. Hashed octree (Warren–Salmon): Morton key → node, in a hash table

- **Good:** O(1) lookup at any level. Memory is proportional to the number of
  nodes, and never sparse.
- **Cost:** a traversal in hash order jumps around memory. It needs a sorted
  list of the leaves beside it for traversals.

### E. Child families plus compaction, keeping a pointer tree

Each node points to one group of 2^d children, allocated together. Whole
families are always created, so nothing is wasted. After each regrid, the
families are copied into Morton order.

- **Good:** the smallest change from the code as it is.
- **Cost:** finding a neighbour goes back to walking up and down the tree, so
  the arithmetic lookups of A are lost.

## Comparison

| | Lookup / neighbour | Traversal locality | Wasted space | Refining | GPU / MPI |
|---|---|---|---|---|---|
| A. per-level arrays (now) | O(1) arithmetic | good, in clustered regions | virtual, per page | O(1) per node | awkward |
| dense tiles | O(1) | excellent | high when refinement is sparse | O(1) | good |
| B. tiles + bitmask/rank | O(1) with `popcount` | excellent | none | batched rebuild | very good |
| C. linear octree | O(log n), or a hash | excellent | none | batched O(n) | best |
| D. hashed | O(1) expected | poor alone | none | O(1) | fair |
| E. families + compaction | walk the tree | good after compaction | none | O(1) | fair |
| F. 64-way node tree (`voxel-memindex`) | ceil(L/B) loads; one with a cursor | good | per node of 64 | O(depth) | fair |
| G. VDB tree (`nano-vdb`) | 4 loads (root hash, upper, lower, leaf); one with the accessor | good | dense: 32^d-slot tables; compact: masks only | O(1) per level; compact: a shift in a node | good (fixed depth) |

## Measured: memindex against voxel-memindex

`voxel-memindex.[ch]` is a prototype of the packed-node idea (F above). It is
not wired into `tree.c` or the `Makefile`.
- Each node is 64 pointers, 512 bytes, from a pool. It takes B bits of every
  index: B = 2 in 3D (two tree levels), 3 in 2D, 6 in 1D. A level-L array is
  ceil (L / B) nodes deep.
- Its interface is the memindex one, prefixed `vmem_`, plus:
  - `vmem_get ()`, which finds and reads an entry in one walk, where
    memindex needs `mem_allocated ()` and then `mem_data ()`;
  - `vmem_get_cursor ()`, which remembers the last bottom node, so that an
    index in the same group of 64 costs a single load.
- Empty nodes go back to the pool as soon as their last entry is freed.

`tests/bench_memindex.sh` compares them (see Tests). The numbers below come
from an i7-8650U laptop with about 2.3 GB free, Linux 6.8, transparent huge
pages set to `[madvise]`, and memindex with `MADV_NOHUGEPAGE`, on a quiet
machine (2026-11-08). Each case ran once. `tests/bench_memindex.md` now holds
the later run with nano-vdb too (see
[With nano-vdb](#with-nano-vdb)).

| case (points) | B/point mem / voxel | build ns/point mem / voxel | shuffled lookup ns mem / voxel | neighbours ns mem / voxel+cursor | absent entry ns mem / voxel |
|---|--:|--:|--:|--:|--:|
| 2D dense L11 (4.2M) | 8 / 8 | 9 / 13 | 16 / 56 | 4.2 / 4.8 | — |
| 2D river L16 (2.1M) | 168 / 19 | 184 / 23 | 62 / 105 | 6.8 / 4.8 | 934 / 8 |
| 2D random L20 (50k) | 8370 / 1907 | 6398 / 1351 | 118 / 268 | 44 / 67 | 103 / 49 |
| 3D dense L7 (2.1M) | 32 / 8 | 48 / 15 | 32 / 69 | 6.1 / 6.3 | — |
| 3D shell L9 (0.5M) | 1027 / 47 | 1915 / 33 | 54 / 72 | 13.6 / 9.8 | 42 / 19 |
| 3D shell L10 (2.1M) | out of memory / 47 | — / 33 | — / 114 | — / 10.2 | — / 17 |
| 3D river L12 (0.85M) | 344 / 26 | 594 / 34 | 86 / 137 | 8.7 / 12.2 | 26 / 15 |
| 3D random L20 (30k) | 16715 / 3636 | 12325 / 2974 | 291 / 629 | 76 / 172 | 72 / 37 |

What it shows:
- **Memory.** voxel takes 3–60 times less wherever a level is not full, and
  about the same, or less in 3D, on a full level. memindex pays a page, and
  part of a page-table page, per row touched. On the 3D shell at level 10 that
  came to more than the 2.3 GB free, where voxel needed 99 MB.
- **Building.** voxel is 5–60 times faster on sparse levels, because it takes
  far fewer page faults. Both are about the same on full levels.
- **Finding an entry that exists, in random order.** memindex is 1.3–3 times
  faster, more so on deeper levels. voxel makes ceil (L / B) dependent loads,
  where memindex makes one per dimension. On a full level in Morton order,
  memindex is about 3 times faster (2.5 against 7.4 ns in 2D).
- **Neighbours.** With the cursor, voxel is about as fast as memindex on full
  and clustered levels, and faster on the 3D shell and the 2D river. It is
  1.5–2.3 times slower on isolated points (random, deep levels).
- **Absent entries** (indices drawn uniformly over a sparse level). voxel is
  2–120 times faster. memindex reads untouched pages, and each read is a page
  fault on the zero page.
- **After freeing everything.** memindex keeps its page tables and its upper
  rows (about 570–590 MB on the random level-20 cases). voxel keeps its freed
  nodes, hot, in the pool, until `pool_trim ()` or `vmem_destroy ()`.
- **voxel's weak case:** isolated points on deep levels, at 2–4 kB per point.
  Each point owns most of the nodes on its path.

In libb-amr the tree indexes blocks, not cells. Stencils run inside padded
blocks, and `block-ghosts` gathers each block's neighbourhood once. So the
slower lookups of voxel are mostly off the hot path, while its memory and
build costs are on it, at every regrid. Confirming that needs voxel wired into
`tree.c` (see Future plans).

### On 2:1 balanced trees

The random cases above have isolated points, which a 2:1 balanced tree never
has. In such a tree, families are whole and lie inside one bottom node, and
every refined cell has neighbours around it. `tests/bench_balanced.sh` builds
balanced trees as `tree.c` does, each level in its own array, and times what
libb-amr actually does with them. The table below is from 2026-11-09, on a
quiet machine; `tests/bench_balanced.md` now holds the later run with
nano-vdb too.

| case (nodes) | B/node mem / voxel | Morton dfs ns/node: mem / voxel / voxel+cursor | neighbour box ns/leaf: mem / voxel / voxel+cursor |
|---|--:|--:|--:|
| 2D uniform, depth 10 (1.4M) | 9 / 8 | 5.5 / 8.1 / 4.7 | 38 / 61 / 52 |
| 2D circle, depth 15 (0.64M) | 615 / 25 | 20.4 / 13.8 / 7.1 | 106 / 101 / 77 |
| 2D river, depth 15 (1.3M) | 234 / 23 | 11.7 / 12.5 / 6.6 | 71 / 92 / 73 |
| 3D uniform, depth 6 (0.3M) | 77 / 8 | 32.9 / 25.2 / 7.3 | 272 / 310 / 274 |
| 3D sphere, depth 9 (2.3M) | 252 / 14 | 44.5 / 28.3 / 8.5 | 310 / 462 / 392 |
| 3D river, depth 10 (0.37M) | 446 / 18 | 38.9 / 28.7 / 8.6 | 321 / 472 / 393 |

- **Morton traversal.** With a cursor per level, voxel is faster in every
  case: 1.2x on the uniform 2D tree, 1.8-2.9x on the 2D fronts, and 4.5-5.3x
  in 3D. A family's children are in one voxel node, while memindex reads them
  from 2 (2D) or 4 (3D) rows. Over a level list, as `ghosts_update ()` walks,
  voxel with a cursor is 1.2-3x faster.
- **Neighbour box** (3^d - 1 neighbours per leaf, about a fifth of them
  falling back to a coarser level). memindex keeps an edge here. With
  cursors, voxel ranges from 28% faster (2D circle) to 36% slower (2D
  uniform), and is up to 26% slower in 3D. The worst gap is about 80 ns per
  leaf. That is an estimated percent or less of the ghost copies that
  follow, for blocks of 16^3.
- **Memory and build** are as above: 10-25x less memory, and 12-25x faster
  builds, on every non-uniform tree.

So in a Morton-ordered pass over a balanced tree, voxel's slower lookups do
not show; they show only in the neighbour box, and there by less than the
work done per block. A box cursor would narrow that gap: it would hold the
up to 2^d bottom nodes that a leaf's box can touch, fetched once per leaf.

### With nano-vdb

Both benchmarks were run again (2026-11-10) with nano-vdb, dense and compact
(`-DNVDB_COMPACT`), all built with `-mpopcnt`. A browser was busy on the
machine, so the times are slower than above; compare implementations within a
case. Full tables and readings: `tests/bench_memindex.md`,
`tests/bench_balanced.md`.

| balanced tree (nodes) | B/node mem / voxel / nvdb / nvdb-compact | dfs+cursor ns/node mem / voxel / nvdb / compact | neighbour box+cursor ns/leaf mem / voxel / nvdb / compact |
|---|--:|--:|--:|
| 2D circle, depth 15 (0.64M) | 615 / 25 / 30 / 13 | 54 / 11.6 / 9.8 / 11.4 | 228 / 135 / 103 / 107 |
| 2D river, depth 15 (1.3M) | 234 / 23 / 28 / 13 | 46 / 10.3 / 12.5 / 10.5 | 207 / 95 / 100 / 98 |
| 3D sphere, depth 9 (2.3M) | 252 / 14 / 26 / 12 | 53 / 20.7 / 27.2 / 25.0 | 629 / 614 / 561 / 604 |
| 3D river, depth 10 (0.37M) | 446 / 18 / 44 / 14 | 45 / 21.0 / 26.7 / 24.5 | 657 / 624 / 587 / 612 |

(memindex has no cursor: its dfs and box columns are plain lookups.)

- **Compact nano-vdb takes the least memory** on every non-uniform tree,
  12-14 B per node, and on sparse single levels (12-22 B per point on
  shells and rivers). Its one weak case is scattered 3D points on deep
  levels (7.4 KB a point at level 20: the 4 KiB of mask bits of a 3D upper
  node), which a balanced tree does not have.
- **Dense nano-vdb** costs 26-44 B per node (its 32^3-slot upper tables),
  frees fastest, and is hopeless for scattered 3D points (about 300 KB
  each).
- **The whole-path accessor makes nano-vdb the fastest on the 3D neighbour
  box**: 11-19% faster than memindex, where voxel's bottom-only cursor is
  2-13% faster.
- **Traversal** with cursors: the three trees are close (voxel usually
  first), all 1.7-5.5 times faster than memindex.

## Recommendation

1. **For now, keep A for the structure.** Spend the effort on the data
   instead: sibling-group allocation and Morton-ordered slabs for blocks.
   That is where the memory traffic is.
2. **When a GPU or MPI is needed, add C beside A.** C is a sorted array of
   leaf keys, rebuilt at each regrid. A stays the fast way to find
   neighbours, and C becomes the contiguous list for GPU launches and for
   splitting the work across ranks. Each regrid already traverses the whole
   tree, so keeping the two in step is cheap.
3. **For tiles, build them as B,** with bitmasks and rank indexing. Pick N by
   dimension: N = 2 in 3D (a tile of 73 slots), and N = 3 or 4 in 2D. Tiles
   suit regions that are refined fairly evenly, and the bitmask removes the
   cost of those that are not.

**To decide, measure:** prototype C beside the current tree, and compare the
two on a benchmark that refines and traverses.

**Measured since:** F (`voxel-memindex`) against A, in isolation (see
[Measured: memindex against voxel-memindex](#measured-memindex-against-voxel-memindex)).
F wins on memory and build time on every level that is not full. A wins on
random lookups, by 1.3–3 times. This revises point 1: replace A by F, once a
run of `adapt` and the ghosts on a sparse case confirms that the lookups stay
off the hot path (Future plans, Tree structure).
