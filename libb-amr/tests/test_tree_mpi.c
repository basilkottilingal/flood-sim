/* tree () on several processes (-D_MPI): MPI started, the tree refined at
 * once to the first level with 4 npe families of leaves at least, and
 * distributed. Against what that must be: every leaf of this process at
 * that level, the leaves over the processes the whole level once, within
 * one family of each other (the cut never splits one), owned as the first
 * keys say; a cell scalar, an integer of the
 * cell's global index, exact in every ghost inside the domain after
 * ghosts_update () (a ghost missing from the halo, or not sent, is wrong);
 * a reduction over the cells the whole level once. Build with -D_MPI. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "tree-mpi.h"
#include "block-ghosts.h"

scalar a;

static int checks, fails, rank, npe;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("rank %d: FAIL line %d: %s\n", rank, __LINE__, #c); } } while (0)

/* the value of the cell of global indices g: exact in float too */
static Real value (const long * g)
{
  return (Real) (g [0] + 37 * g [1] + 1369 * g [2]);
}

int main (void)
{
  tree (5, 4);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  CHECK (tree_npe () == npe && tree_rank () == rank);

  /* the level: the first with 4 npe families of leaves at least, below
   * depth 5; on one process tree () as always, level 1 */
  int L = 1;
  while (npe > 1 && L < 5 && (1L << (dimension * (L - 1))) < 4L * npe)
    L++;
  const TreeCaches * tc = tree_caches ();
  int level = 1, keyed = 1;
  for (int q = 0; q < tc->nowned; q++) {
    level &= tc->owned [q].p.level == L;
    keyed &= tree_partition_owner (tc->owned [q].p) == rank;
  }
  long n = tc->nowned, total, most, least;
  MPI_Allreduce (&n, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&n, &most, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&n, &least, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
  /* within one family of each other: the cut never splits one */
  CHECK (level && keyed && total == 1L << (dimension * L)
         && most - least <= TREE_NCHILDREN && least > 0);
  CHECK (npe == 1 ? !tree_mpi_plan () : tree_mpi_plan () != NULL);
  /* the counts the tree keeps, against its flags */
  int npr = 0, nlf = 0;
  for (int q = 0; q < tc->nall; q++) {
    npr += (tc->all [q].node->flags & BAMR_NODE_IS_PRUNED) != 0;
    nlf += (tc->all [q].node->flags & BAMR_NODE_LEAF) != 0;
  }
  CHECK (npr == tree_pruned () && nlf == tree_get_nleaves ()
         && tree_get_nblocks () == tc->nall);

  /* the cells by their global indices; then every ghost inside the domain */
  BlockShape sh = tree_block_shape ();
  int m = sh.m, g = sh.g;
  long side = (long) m << L;
  foreach ()
    a[] = value ((long []) { (long) floor (x / Delta), (long) floor (y / Delta),
                             (long) floor (z / Delta) });
  ghosts_update (a);
  long compared = 0, wrong = 0;
  foreach_leaf () {
    const Real * b = block_data (sblocks, a);
    long o [3] = { (long) i * m, 0, 0 };
#if dimension > 1
    o [1] = (long) j * m;
#endif
#if dimension > 2
    o [2] = (long) k * m;
#endif
    for (long c = 0; c < sh.cells; c++) {
      int x [3] = { 0, 0, 0 };
      long r = c;
      x [0] = (int) (r % sh.nx) - sh.ox;
#if dimension > 1
      r /= sh.nx;
      x [1] = (int) (r % (m + 2 * g)) - g;
#endif
#if dimension > 2
      r /= (m + 2 * g);
      x [2] = (int) r - g;
#endif
      long gi [3] = { 0, 0, 0 };
      int in = 1, ghost = 0;
      for (int ax = 0; ax < dimension; ax++) {
        in &= x [ax] >= -g && x [ax] < m + g;
        gi [ax] = o [ax] + x [ax];
        in &= gi [ax] >= 0 && gi [ax] < side;
        ghost |= x [ax] < 0 || x [ax] >= m;
      }
      if (!in || !ghost)
        continue;
      compared++;
      wrong += b [c] != value (gi);
    }
  }
  long w [2] = { compared, wrong }, ws [2];
  MPI_Allreduce (w, ws, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  CHECK (ws [1] == 0 && ws [0] > 0);

  /* the cells counted once over the processes, with the value before */
  long cells = 3;
  foreach (reduction (+:cells))
    cells++;
  CHECK (cells == 3 + (1L << (dimension * L)) * (long) pow (m, dimension));

  tree_destroy ();
  int f = fails, ftot;
  MPI_Allreduce (&f, &ftot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf ("dimension %d, %d processes: level %d, %ld leaves, %ld ghosts "
            "compared, %ld wrong; %d checks, %d failed\n", dimension, npe, L,
            total, ws [0], ws [1], checks, ftot);
  return ftot != 0;
}
