#ifndef _BAMR_TREE_MPI_H_
#define _BAMR_TREE_MPI_H_

  #include <stdint.h>
  #include "tree.h"

  /*
  .. The tree shared among npe processes: a balanced partition, what each
  .. process sends and receives for the ghosts, and the subtrees a process
  .. does without, pruned away. Every process builds the same whole tree
  .. first. The partition, the plans and the pruning are computed from it
  .. with no MPI call, so they run and are tested in one process for any
  .. npe. With -D_MPI, tree_mpi_partition () applies them on
  .. MPI_COMM_WORLD: rank 0 partitions and sends the result to ranks 1 to
  .. npe - 1, each process makes its plan and prunes its tree, and
  .. ghosts_update () exchanges blocks between the processes.
  ..
  .. The processes are tree_npe () and tree_rank () (tree.h), the owner of a
  .. node its node->pid.
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
  .. smallest pid below it.
  ..
  .. The partition is the tree's at the time: a refinement, a coarsening or
  .. a pruning makes it stale (tree_generation); tree_partition_valid ()
  .. tells. It needs the whole tree, not pruned.
  ..
  ..   tree_partition (npe, weight)  partitions: sets every node's pid;
  ..                                 weight [q] for leaf q of tree_caches
  ..                                 ()->leaves, each > 0, or NULL: 1 each.
  ..                                 Returns 0, -1 for npe < 1, a weight not
  ..                                 > 0 or a pruned tree (nothing changed)
  ..   tree_partition_owner (p)      the process whose piece of the curve
  ..                                 holds p, from the first keys alone
  ..                                 (key_first): a leaf's, its pid; a
  ..                                 node's, its first leaf's, so its pid
  ..                                 too. p need not be in this process's
  ..                                 tree. With no partition made,
  ..                                 tree_rank (): the tree is its maker's
  ..   tree_get_partition ()         the last one made, NULL if none
  ..   tree_partition_valid ()       1 if made since the tree last changed
  ..   tree_partition_free ()        frees it (tree_destroy () does not)
  ..
  ..
  .. 2. The neighbours of a node, as bit maps
  .. ----------------------------------------
  ..
  .. The 3^dimension - 1 directions around a node p of level l, numbered as
  .. block_ghosts () visits the regions of its ghosts: d = (d0, d1, d2), each
  .. -1, 0 or 1, is bit (d0 + 1) + 3 (d1 + 1) + 9 (d2 + 1), the centre's bit
  .. (1, 4, 13 in 1, 2, 3 dimensions) never set. Under 2:1 balance, the
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
  .. and pid [dir] the rank of the node read there (p + d, or the coarse
  .. leaf covering it): its node->pid, -1 outside. Made from the tree's
  .. structure: on a pruned tree a pruned node is "same", having no
  .. children; the walls are those of no bit, ~(same | fine | coarse).
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
  ..                              partition (valid: the tree whole)
  ..   tree_mpi_plan_free (plan)
  ..
  ..
  .. 4. Pruning
  .. ----------
  ..
  .. A process needs its own nodes, the nodes it receives, and their
  .. ancestors. Every subtree with none of them below its root is pruned:
  .. tree_prune () (tree.h) frees its nodes and their blocks. The root
  .. stays, with its blocks and its pid, as a pruned node: no children, not
  .. a leaf. Only the largest such subtrees are pruned, so no pruned node is
  .. below another; a family whose one child is needed keeps the others,
  .. families being allocated together. It cannot be undone, and the tree
  .. cannot be refined or coarsened after.
  ..
  ..   tree_mpi_prune (plan)   prunes for plan's process, right after the
  ..                           plan was made (the tree unchanged since):
  ..                           the subtrees pruned
  ..
  ..
  .. 5. With MPI (-D_MPI): the skeleton
  .. ----------------------------------
  ..
  .. A distributed process holds its skeleton: its own nodes, their
  .. ancestors, and N (p) of each own node p: p's same-level neighbours,
  .. the coarser leaf covering an absent one, for a leaf the children of a
  .. neighbour with children that touch it, p's parent and children. That
  .. is every node the plans read, both ways: N (x) holds p when N (p)
  .. holds x (2:1 balance), so the ranks of N () of its nodes, its peers,
  .. name it too, and every exchange is point to point with them, sizes
  .. first: no message to all. A node of the skeleton it neither owns nor
  .. receives is inactive (BAMR_NODE_IS_INACTIVE, tree.h): flags and pid,
  .. no blocks when its whole family is (blocks go by family). The
  .. owners of the neighbours are their pids, known here.
  ..
  ..   tree_mpi_partition (weight)  on MPI_COMM_WORLD, every process with
  ..                                the same whole tree: rank 0 partitions
  ..                                (its weights; the others' are not read)
  ..                                and sends every node's pid and the
  ..                                loads to ranks 1 to npe - 1 (MPI_Bcast);
  ..                                each prunes to its skeleton, makes its
  ..                                plan (receives from the tree here,
  ..                                sends as its peers ask) and its nodes
  ..                                active or not; the exchanges on, every
  ..                                scalar marked stale.
  ..                                Returns 0, -1 if rank 0's tree_partition
  ..                                () refused (nothing changed)
  ..   tree_mpi_repartition (weight)
  ..                                on the distributed tree, collective: a
  ..                                new cut by weight [q] for each leaf q
  ..                                of tree_caches ()->owned (NULL: 1
  ..                                each), by prefix sums (MPI_Exscan),
  ..                                each leaf where tree_partition () on
  ..                                the whole tree would put it; every
  ..                                node that changes hands sent to its
  ..                                new owner with its blocks (every
  ..                                scalar) and N () of it (structure),
  ..                                made there (tree_insert ()), between
  ..                                the old and new owners of each piece
  ..                                of the curve (from the first keys of
  ..                                the two partitions); pids by the new
  ..                                keys; then the skeleton and plan again.
  ..                                Returns the nodes this process sent,
  ..                                -1 if a weight is not > 0 anywhere or
  ..                                no partition was made (nothing done).
  ..                                Integer weights give the whole tree's
  ..                                cut exactly; others may differ at a
  ..                                tie, the sums being in another order
  ..   tree_mpi_end ()              exchanges off, the plan and partition
  ..                                freed; the tree stays pruned
  ..   tree_mpi_plan ()             the plan applied, NULL if none
  ..   tree_mpi_exchange (parents, kind, level, n, s)
  ..                                one exchange of the plan, the blocks of
  ..                                scalars s [0..n): what ghosts_update_list
  ..                                () calls
  ..
  .. Not yet: adapting under MPI; ghosts_update_level () (multigrid),
  .. block_ghosts () alone.
  */

  /* 1. the partition */
  typedef struct
  {
    int npe;                /* the processes                             */
    int nleaves;            /* the leaves partitioned                    */
    int * leaf_first;       /* process r holds leaves leaf_first [r] ..
                               leaf_first [r + 1] - 1 of tree_caches ()->
                               leaves, of the whole tree; npe + 1        */
    double * load;          /* the weight of each process's leaves       */
    uint64_t * key_first;   /* the Morton key of each process's first
                               leaf (tree_morton_key ()), the next one's
                               for a process with none; npe + 1, the
                               last UINT64_MAX                           */
    double imbalance;       /* the largest load over the mean, >= 1      */
    unsigned generation;    /* tree_generation when made                 */
  } TreePartition;

  int                   tree_partition (int npe, const double * weight);
  int                   tree_partition_owner (Index p);
  const TreePartition * tree_get_partition (void);
  int                   tree_partition_valid (void);
  void                  tree_partition_free (void);

  /* 2. the neighbours */
  #define TREE_NDIRS (dimension == 1 ? 3 : dimension == 2 ? 9 : 27)
  #define TREE_DIR_CENTRE (TREE_NDIRS / 2)
  typedef struct
  {
    uint32_t same, fine, coarse;
    int pid [TREE_NDIRS];   /* the rank of the node read; -1 outside     */
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
    int depth;              /* the deepest level of the whole tree       */
    int depth_min;          /* its coarsest leaf's                       */
    unsigned generation;    /* tree_generation when made                 */
    int npeers, * peers;    /* the ranks exchanged with, sorted: those
                               of N () of this process's nodes (the
                               distributed plan's; none in a whole-tree
                               plan of tree_mpi_plan_new ())             */
    TreeMpiExchange ex [2][TREE_MPI_KINDS][BAMR_MAX_DEPTH + 2];
                            /* [parents][kind][level]                    */
  } TreeMpiPlan;

  TreeMpiPlan * tree_mpi_plan_new (int rank);
  void          tree_mpi_plan_free (TreeMpiPlan * plan);

  /* 4. pruning */
  int tree_mpi_prune (const TreeMpiPlan * plan);

  /* 5. MPI */
  #ifdef _MPI
  int                 tree_mpi_partition (const double * weight);
  int                 tree_mpi_repartition (const double * weight);
  void                tree_mpi_end (void);
  const TreeMpiPlan * tree_mpi_plan (void);
  void                tree_mpi_exchange (int parents, int kind, int level,
                                         int n, const int * s);
  #endif

#endif
