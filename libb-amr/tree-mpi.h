#ifndef _BAMR_TREE_MPI_H_
#define _BAMR_TREE_MPI_H_

  #include <stdint.h>
  #include "tree.h"

  /*
  .. The tree shared among npe processes: a balanced partition, what each
  .. process sends and receives for the ghosts, and the subtrees a process
  .. does without. The first step to MPI (README, "MPI"): every process
  .. holds the whole tree, as now; the partition, the plans and the pruning
  .. are computed from it in the same way everywhere, with no MPI call, so
  .. they run and are tested in one process for any npe. With -D_MPI,
  .. tree_mpi_partition () applies them on MPI_COMM_WORLD, and
  .. ghosts_update () exchanges blocks between the processes.
  ..
  ..
  .. 1. The partition
  .. ----------------
  ..
  .. The leaves are cut along the Morton curve (tree_caches ()->leaves) into
  .. npe pieces, one a process, each a contiguous run: the curve keeps a
  .. piece compact, so its surface, what it exchanges, stays small. A leaf
  .. has a weight, its cost (1 each by default). Leaf q, its weight w_q, the
  .. weights before it W_q, the total W, goes to
  ..
  ..   pid = floor (npe (W_q + w_q / 2) / W)
  ..
  .. the process whose share of the curve holds the leaf's middle. So the
  .. pieces follow the curve in order (pid never decreases along it), and
  .. none weighs more than W / npe + the largest weight. With unit weights,
  .. leaf q of n goes to floor (npe (2 q + 1) / (2 n)), computed in
  .. integers: the pieces differ by one leaf at most, and none is empty if
  .. n >= npe.
  ..
  .. A parent goes to the process of its first leaf (in Morton order), the
  .. smallest pid below it. A node's process is its node->pid.
  ..
  .. The partition is the tree's at the time: a refinement or a coarsening
  .. makes it stale (tree_generation), and new nodes have pid 0. Make it
  .. again after adapting; tree_partition_valid () tells. It is made on the
  .. whole tree: a pruned tree is unpruned first.
  ..
  ..   tree_partition (npe, weight)  partitions; weight [q] for leaf q of
  ..                                 tree_caches ()->leaves, each > 0, or
  ..                                 NULL: 1 each. Returns 0, -1 for npe < 1
  ..                                 or a weight not > 0 (nothing changed)
  ..   tree_get_partition ()         the last one made, NULL if none
  ..   tree_partition_valid ()       1 if made since the tree last changed
  ..                                 (pruning and unpruning keep it valid)
  ..   tree_partition_free ()        frees it (tree_destroy () does not)
  ..
  ..
  .. 2. The neighbours of a node, as bit maps
  .. ----------------------------------------
  ..
  .. The 3^dimension - 1 directions around a node p of level l, numbered as
  .. block_ghosts () visits the regions of its ghosts: d = (d0, d1, d2), each
  .. -1, 0 or 1, is bit (d0 + 1) + 3 (d1 + 1) + 9 (d2 + 1), the centre's bit
  .. (1, 3, 13 in 1, 2, 3 dimensions) never set. Under 2:1 balance, the
  .. index p + d is one of four, a bit in one of three maps or in none:
  ..
  ..   same     a leaf of level l: p's ghosts there copy its cells
  ..   fine     a parent of level l, its children leaves of level l + 1:
  ..            the ghosts copy its restriction; a face's reflux reads its
  ..            children
  ..   coarse   no node of level l: the leaf of level l - 1 covering it,
  ..            from which the ghosts there are prolonged
  ..   none     outside the domain: p's boundary conditions
  ..
  .. Made from the tree's structure (the children a node has, pruned or
  .. not); the walls are those of no bit, ~(same | fine | coarse).
  ..
  ..   tree_neighbours (p)   the maps of the node at p
  ..
  ..
  .. 3. What is sent and received: the plan
  .. --------------------------------------
  ..
  .. A process r updates the ghosts of its own nodes (node->pid == r), and
  .. restricts its own parents, as ghosts_update_list () does on one
  .. process (block-ghosts.h), at each level in the same order. What it
  .. reads of other processes' nodes is sent to it whole, as blocks
  .. (interior and ghosts, every scalar of the update), before the step that
  .. reads them, by their owners. Three kinds of exchange, at a level l:
  ..
  ..   TREE_MPI_RESTRICT  before the parents of level l are restricted: the
  ..                      children (level l + 1) of r's parents that are
  ..                      not r's, made by their owners (restricted, if
  ..                      parents, at level l + 1 just before)
  ..   TREE_MPI_GHOSTS    before the ghosts of level l: for each of r's
  ..                      nodes of level l, the nodes of its same and fine
  ..                      bits (level l: their cells, restricted already),
  ..                      and those of its coarse bits (level l - 1: their
  ..                      ghosts made at level l - 1 just before, for the
  ..                      prolongation, which reads them)
  ..   TREE_MPI_FACES     with the ghosts of level l, face scalars only: for
  ..                      each of r's leaves of level l whose lower
  ..                      neighbour along an axis a is a parent, that
  ..                      parent's upper children along a, whose faces m
  ..                      the reflux reads (block-ghosts-impl.h)
  ..
  .. Each a list per peer, in the same order on both sides (by level, then
  .. k, j, i), so that what r sends to s is what s expects from r: a
  .. message per peer and exchange, of known length, with no size sent
  .. first. A plan for the ghosts of the leaves (ghosts_update ()) and one
  .. for every node too (ghosts_update_levels ()): parents [0] and [1].
  ..
  ..   tree_mpi_plan_new (rank)   the plan of process rank, from the
  ..                              partition (valid, the tree unpruned)
  ..   tree_mpi_plan_free (plan)
  ..
  ..
  .. 4. Pruning
  .. ----------
  ..
  .. A process needs its own nodes, the nodes it receives, and their
  .. ancestors. A subtree with none of them is pruned: its root gets
  .. BAMR_NODE_IS_PRUNED, so that tree_has_children () says no, and the
  .. tree's lists (tree_caches ()), foreach_node (), foreach_leaf () and
  .. morton () skip what is below. The nodes and their blocks stay (the
  .. tree is still whole on each process): a scalar made while the tree is
  .. pruned gets blocks below too (grid-scalar.c), so unpruning restores the
  .. tree as it was. Refining and coarsening need the tree unpruned.
  ..
  ..   tree_prune (plan)     prunes for plan's process; the roots pruned
  ..   tree_unprune ()       every node back; how many were pruned
  ..   tree_pruned ()        the roots pruned now
  ..
  ..
  .. 5. With MPI (-D_MPI)
  .. --------------------
  ..
  ..   tree_mpi_partition (weight)  on MPI_COMM_WORLD: the partition (the
  ..                                weights the same on every process),
  ..                                this process's plan, pruned, and the
  ..                                exchanges on: every scalar marked stale
  ..   tree_mpi_end ()              exchanges off, unpruned, plan freed
  ..   tree_mpi_plan ()             the plan applied, NULL if none
  ..   tree_mpi_exchange (parents, kind, level, n, s)
  ..                                one exchange of the plan, the blocks of
  ..                                scalars s [0..n): what ghosts_update_list
  ..                                () calls
  ..
  .. Not yet: foreach () and the other loops still visit every leaf left in
  .. the lists (the process's, the ones it receives, and those not pruned
  .. away), so a reduction counts some twice; adapting under MPI;
  .. ghosts_update_level () (multigrid), block_ghosts () alone.
  */

  /* 1. the partition */
  typedef struct
  {
    int npe;                /* the processes                             */
    int nleaves;            /* the leaves partitioned                    */
    int * leaf_first;       /* process r holds leaves leaf_first [r] ..
                               leaf_first [r + 1] - 1 of tree_caches ()->
                               leaves, of the tree unpruned; npe + 1     */
    double * load;          /* the weight of each process's leaves       */
    double imbalance;       /* the largest load over the mean, >= 1      */
    unsigned generation;    /* tree_generation when made                 */
  } TreePartition;

  int                   tree_partition (int npe, const double * weight);
  const TreePartition * tree_get_partition (void);
  int                   tree_partition_valid (void);
  void                  tree_partition_free (void);

  /* 2. the neighbours */
  #define TREE_NDIRS (dimension == 1 ? 3 : dimension == 2 ? 9 : 27)
  #define TREE_DIR_CENTRE (TREE_NDIRS / 2)
  typedef struct
  {
    uint32_t same, fine, coarse;
  } TreeNeighbours;

  TreeNeighbours tree_neighbours (Index p);

  /* 3. the plan */
  enum { TREE_MPI_RESTRICT, TREE_MPI_GHOSTS, TREE_MPI_FACES, TREE_MPI_KINDS };

  typedef struct
  {
    int peer;               /* the other process                         */
    int n;                  /* its nodes                                 */
    TreeItem * item;        /* by level, then k, j, i                    */
  } TreeMpiList;

  typedef struct
  {
    int nsend, nrecv;
    TreeMpiList * send, * recv;    /* by peer, ascending */
  } TreeMpiExchange;

  typedef struct
  {
    int rank, npe;
    int depth;              /* the deepest level of the tree unpruned    */
    int depth_min;          /* its coarsest leaf's                       */
    TreeMpiExchange ex [2][TREE_MPI_KINDS][BAMR_MAX_DEPTH + 2];
                            /* [parents][kind][level]                    */
  } TreeMpiPlan;

  TreeMpiPlan * tree_mpi_plan_new (int rank);
  void          tree_mpi_plan_free (TreeMpiPlan * plan);

  /* 4. pruning */
  int tree_prune (const TreeMpiPlan * plan);
  int tree_unprune (void);
  int tree_pruned (void);

  /* 5. MPI */
  #ifdef _MPI
  int                 tree_mpi_partition (const double * weight);
  void                tree_mpi_end (void);
  const TreeMpiPlan * tree_mpi_plan (void);
  void                tree_mpi_exchange (int parents, int kind, int level,
                                         int n, const int * s);
  #endif

#endif
