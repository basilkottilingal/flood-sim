/* Stored ghost cells: ghosts_update () at level jumps, at the walls and in
 * corners, on parents; the operations of a scalar (prolongation,
 * restriction, boundary condition); derived scalars on the halo; and
 * refining and coarsening keeping the ghosts valid. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cwrap/cache.h>
#include "adapt.h"
#include "block-ghosts.h"
/* operators of its own: this file dispatches, in place of block-ghosts.c */
#include "block-ops.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

scalar u, v, w, z, e, u2;
scalar zbl, hl;           /* a bed and a water depth on it: declared last,
                             hl takes the lower id, so that its bed coming
                             first is the order of dependence, not of ids */

static const double tol = 1e3 * REAL_EPSILON;

/* the block coordinates of p, as an array */
static void coords (Index p, int P [3])
{
  P [0] = p.i, P [1] = 0, P [2] = 0;
#if dimension > 1
  P [1] = p.j;
#endif
#if dimension > 2
  P [2] = p.k;
#endif
}

/* the interior coordinates of cell (ci, cj, ck) of a stored block */
#if dimension == 1
# define XOF(x) int x [3] = { ci - shape.ox, 0, 0 }
#elif dimension == 2
# define XOF(x) int x [3] = { ci - shape.ox, cj - shape.g, 0 }
#else
# define XOF(x) int x [3] = { ci - shape.ox, cj - shape.g, ck - shape.g }
#endif

/* the centre of interior cell x of block p, in [0, 1] (outside for ghosts) */
static void centre (Index p, const int x [3], double r [3])
{
  int P [3];
  coords (p, P);
  double n = (double) (tree_get_block_m () << p.level);
  for (int a = 0; a < 3; a++)
    r [a] = a < dimension ? (P [a] * tree_get_block_m () + x [a] + 0.5) / n : 0.5;
}

static int is_ghost (BlockShape shape, const int x [3])
{
  int gh = 0;
  for (int a = 0; a < dimension; a++)
    gh |= x [a] < 0 || x [a] >= shape.m;
  return gh;
}

static double linear (const double r [3]) { return 1 + 2 * r [0] + 3 * r [1] + 5 * r [2]; }
static double wavy (const double r [3])
{
  return sin (5 * r [0]) + (dimension > 1) * cos (3 * r [1]) + (dimension > 2) * r [2] * r [2];
}
static double kink (const double r [3]) { return fmax (0, r [0] - 0.37) + fmax (0, 0.61 - r [1]); }

/* scalar s = f on the interior of every leaf */
static void set (int s, double (* f) (const double [3]))
{
  BlockShape shape = tree_block_shape ();
  foreach_leaf () {
    Real * b = block_data (sblocks, s);
    block_cells (shape, 0) {
      XOF (x);
      double r [3];
      centre ((Index) {IJKL}, x, r);
      b [c] = (Real) f (r);
    }
  }
}

/* the largest |s - f| over the ghosts of every leaf, or every node */
static double ghost_error (int s, double (* f) (const double [3]), int parents)
{
  BlockShape shape = tree_block_shape ();
  double err = 0;
  morton (TREE_ROOT_INDEX) {
    if (!parents && !(node->flags & BAMR_NODE_LEAF))
      continue;
    const Real * b = block_data (sblocks, s);
    block_cells (shape, shape.g) {
      XOF (x);
      double r [3];
      centre ((Index) {IJKL}, x, r);
      err = fmax (err, fabs (b [c] - f (r)));
    }
  }
  return err;
}

/* a linear extrapolation at the walls: exact for a linear field */
@boundary bc_linear (const Boundary * bc)
{
  boundary_cells (bc) {
    Real w0 = bc->b [ghost + (layer + 1) * inward];
    Real w1 = bc->b [ghost + (layer + 2) * inward];
    bc->b [ghost] = w0 + (layer + 1) * (w0 - w1);
  }
}

/* restrictions counted, to check each parent gets one */
static long restricted;
@restriction restriction_counted (Index p, scalar s)
{
  BLOCK_PRAGMA (omp atomic)
  restricted++;
  @restriction (restriction_mean) (p, s);
}

/* lake at rest: the depth from the surface h + zb, injected, less the fine
.. bed, which is prolonged first */
@prolongation prolong_lake (const Prolongation * pr) (scalar zb)
{
  scalar s = pr->s;
  foreach_prolongation (pr)
    s[] = coarse (s) + coarse (zb) - zb[];
}

/* the core, dispatching to the operators above and those given */
#include "block-ghosts-impl.h"

/* the largest |h + b - surface| over every leaf cell, ghosts included */
static double surface_error (int h, int b, double surface)
{
  BlockShape shape = tree_block_shape ();
  double err = 0;
  foreach_leaf () {
    const Real * x = block_data (sblocks, h), * y = block_data (sblocks, b);
    block_cells (shape, shape.g)
      err = fmax (err, fabs ((double) x [c] + y [c] - surface));
  }
  return err;
}

/* the leaf holding the point r */
static Index leaf_at (const double r [3])
{
  for (int l = 0; ; l++) {
    int n = 1 << l;
    Index q = { (int) (r [0] * n),
#if dimension > 1
                (int) (r [1] * n),
#endif
#if dimension > 2
                (int) (r [2] * n),
#endif
                l };
    if (index_is_allocated (q) && (index_node (q)->flags & BAMR_NODE_LEAF))
      return q;
  }
}

/* the ghost cells of leaves next to a coarser leaf: those that test the
 * prolongation */
static long jump_ghosts (void)
{
  BlockShape shape = tree_block_shape ();
  long n = 0;
  foreach_leaf () {
    Index p = {IJKL};
    int P [3];
    coords (p, P);
    int nb = 1 << level;
    block_cells (shape, shape.g) {
      XOF (x);
      if (!is_ghost (shape, x))
        continue;
      int q [3] = { 0, 0, 0 }, in = 1;
      for (int a = 0; a < dimension; a++) {
        q [a] = P [a] + (x [a] < 0 ? -1 : x [a] >= shape.m ? 1 : 0);
        in &= q [a] >= 0 && q [a] < nb;
      }
      Index b = { q [0],
#if dimension > 1
                  q [1],
#endif
#if dimension > 2
                  q [2],
#endif
                  level };
      n += in && !index_is_allocated (b);
    }
  }
  return n;
}

/* every leaf's block of scalar s, ghosts included, copied */
static Real * snapshot (int s, long * n)
{
  long cells = tree_block_shape ().cells, k = 0;
  Real * all = malloc ((size_t) tree_get_nblocks () * cells * sizeof (Real));
  foreach_leaf () {
    memcpy (all + k * cells, block_data (sblocks, s), cells * sizeof (Real));
    k++;
  }
  *n = k * cells;
  return all;
}

/* the largest difference between the leaves of s, ghosts included, as
 * they are and after ghosts_update_scalars () */
static double kept_against_fresh (int s)
{
  long n;
  Real * kept = snapshot (s, &n);
  ghosts_update_scalars (0);
  Real * fresh = snapshot (s, &n);
  double d = 0;
  for (long q = 0; q < n; q++)
    d = fmax (d, fabs (kept [q] - fresh [q]));
  free (kept);
  free (fresh);
  return d;
}

int main (void)
{
  tree (5, 4);
  BlockShape shape = tree_block_shape ();
  CHECK (shape.g == GHOSTS && shape.m == 16);

  /* levels 1 to 4 around one point, 2 around another: level jumps of one,
   * in faces, edges and corners */
  {
    double a [3] = { 0.3, 0.3, 0.3 }, b [3] = { 0.9, 0.1, 0.6 };
    for (int l = 0; l < 3; l++)
      tree_refine (leaf_at (a));
    tree_refine (leaf_at (b));
    CHECK (leaf_at (a).level == 4 && tree_get_depth_min () < 4);
  }
  long jumps = jump_ghosts ();
  CHECK (jumps > 0);

  /* 1. a linear field, extrapolated at the walls: every ghost exact, at
   * level jumps, walls and corners; parents too with ghosts_update_levels */
  ghost_ops (u, boundary.all (bc_linear));
  set (u, linear);
  ghosts_update (u);
  CHECK (ghost_error (u, linear, 0) < tol);
  ghosts_update_levels (u);
  CHECK (ghost_error (u, linear, 1) < tol);

  /* 2. the defaults: zero gradient at the walls, mirrored over every wall a
   * ghost is beyond; bc_odd the same, with a sign per wall */
  ghost_ops (w, boundary.all (bc_odd));
  set (v, wavy);
  set (w, wavy);
  ghosts_update (v, w);
  {
    double dn = 0, dodd = 0;
    long walls = 0;
    foreach_leaf () {
      int P [3];
      coords ((Index) {IJKL}, P);
      int nb = 1 << level;
      const Real * bv = block_data (sblocks, v), * bw = block_data (sblocks, w);
      block_cells (shape, shape.g) {
        XOF (x);
        int y [3] = { x [0], x [1], x [2] }, flips = 0;
        for (int a = 0; a < dimension; a++) {
          if (P [a] == 0 && x [a] < 0)
            y [a] = -1 - x [a], flips++;
          if (P [a] == nb - 1 && x [a] >= shape.m)
            y [a] = 2 * shape.m - 1 - x [a], flips++;
        }
        if (!flips)
          continue;
        walls++;
        long mc = block_at (shape, y);
        dn = fmax (dn, fabs (bv [c] - bv [mc]));
        dodd = fmax (dodd, fabs (bw [c] - (flips & 1 ? -1 : 1) * bw [mc]));
      }
    }
    CHECK (walls > 0 && dn == 0 && dodd == 0);
  }

  /* 2b. per face: odd at every wall but left and right, zero gradient
   * there; a ghost is its mirror, the sign flipped once per wall crossed
   * but those two, corners too */
  {
    scalar pf;
    ghost_ops (pf, boundary.all (bc_odd), boundary.left (bc_neumann),
                   boundary.right (bc_neumann));
    CHECK (scalar_get_boundary (pf, GHOST_LEFT) == bc_neumann
           && scalar_get_boundary (pf, GHOST_RIGHT) == bc_neumann);
    for (int f = 2; f < GHOST_FACES; f++)
      CHECK (scalar_get_boundary (pf, f) == bc_odd);
    set (pf, wavy);
    ghosts_update (pf);
    double d = 0;
    long walls = 0, both = 0;
    foreach_leaf () {
      int P [3];
      coords ((Index) {IJKL}, P);
      int nb = 1 << level;
      const Real * b = block_data (sblocks, pf);
      block_cells (shape, shape.g) {
        XOF (x);
        int y [3] = { x [0], x [1], x [2] }, flips = 0, odd = 0;
        for (int a = 0; a < dimension; a++)
          if ((P [a] == 0 && x [a] < 0) || (P [a] == nb - 1 && x [a] >= shape.m)) {
            y [a] = x [a] < 0 ? -1 - x [a] : 2 * shape.m - 1 - x [a];
            flips++, odd += a > 0;
          }
        if (!flips)
          continue;
        walls++, both += flips > odd && odd > 0;
        d = fmax (d, fabs (b [c] - (odd & 1 ? -1 : 1) * b [block_at (shape, y)]));
      }
    }
    CHECK (walls > 0 && d == 0 && (dimension == 1 || both > 0));
  }

  /* 3. injection: a ghost next to a coarser leaf is that leaf's cell, as
   * it is; next to a node of its level, that node's cell */
  scalar_set_ops (z, (ScalarOps) { .prolongation = prolong_injection });
  set (z, linear);
  ghosts_update (z);
  {
    double err = 0;
    long coarse = 0;
    int m = shape.m;
    foreach_leaf () {
      int P [3];
      coords ((Index) {IJKL}, P);
      int nb = 1 << level;
      const Real * b = block_data (sblocks, z);
      block_cells (shape, shape.g) {
        XOF (x);
        int q [3] = { 0, 0, 0 }, in = 1, X [3] = { 0, 0, 0 };
        if (!is_ghost (shape, x))
          continue;
        for (int a = 0; a < dimension; a++) {
          q [a] = P [a] + (x [a] < 0 ? -1 : x [a] >= m ? 1 : 0);
          in &= q [a] >= 0 && q [a] < nb;
          X [a] = P [a] * m + x [a];
        }
        if (!in)
          continue;
        Index qb = { q [0],
#if dimension > 1
                    q [1],
#endif
#if dimension > 2
                    q [2],
#endif
                    level };
        double r [3] = { 0.5, 0.5, 0.5 };
        if (index_is_allocated (qb))
          centre ((Index) {IJKL}, x, r);
        else {
          coarse++;
          double n = (double) (m << (level - 1));
          for (int a = 0; a < dimension; a++)
            r [a] = ((X [a] >> 1) + 0.5) / n;
        }
        err = fmax (err, fabs (b [c] - linear (r)));
      }
    }
    CHECK (coarse == jumps && err < tol);
  }

  /* 4. minmod: exact for a linear field too; and no new minimum, where the
   * linear one undershoots a kink */
  scalar_set_ops (z, (ScalarOps) { .prolongation = prolong_minmod });
  ghost_ops (z, boundary.all (bc_linear));
  set (z, linear);
  ghosts_update (z);
  CHECK (ghost_error (z, linear, 0) < tol);
  ghost_ops (z, boundary.all (bc_neumann));
  set (z, kink);
  ghosts_update (z);
  {
    double lo = 1;
    foreach_leaf () {
      const Real * b = block_data (sblocks, z);
      block_cells (shape, shape.g) lo = fmin (lo, b [c]);
    }
    CHECK (lo >= 0);
  }

  /* 5. a restriction of one's own: once per parent, by ghosts_update */
  scalar_set_ops (v, (ScalarOps) { .restriction = restriction_counted });
  restricted = 0;
  ghosts_update (v);
  CHECK (restricted == tree_get_nblocks () - tree_get_nleaves ());
  restricted = 0;
  tree_restrict (v);
  CHECK (restricted == tree_get_nblocks () - tree_get_nleaves ());
  scalar_set_ops (v, (ScalarOps) {0});

  /* 6. a derived scalar, e = u + 2 u2, computed on the halo from the ghosts
   * of u and u2, is what filling its own ghosts gives: prolongation and
   * restriction are linear */
  ghost_ops (u2, boundary.all (bc_linear));
  ghost_ops (e, boundary.all (bc_linear));
  set (u, wavy);
  set (u2, kink);
  ghosts_update (u, u2);
  foreach_leaf () {
    Real * be = block_data (sblocks, e);
    const Real * a = block_data (sblocks, u), * b = block_data (sblocks, u2);
    block_cells (shape, shape.g) be [c] = a [c] + 2 * b [c];
  }
  {
    long n;
    Real * derived = snapshot (e, &n);
    ghosts_update (e);
    Real * filled = snapshot (e, &n);
    double d = 0;
    for (long q = 0; q < n; q++)
      d = fmax (d, fabs (derived [q] - filled [q]));
    CHECK (d < 10 * tol);
    free (derived);
    free (filled);
  }

  /* 7. refining and coarsening keep every leaf's ghosts what ghosts_update
   * gives, without calling it, for every scalar: first coarsening every
   * family that can be, from values set on the leaves, so that their
   * prolongation from the parent differs from what they held; then
   * refining and coarsening at random */
  set (v, wavy);
  set (u, linear);
  ghosts_update_scalars (0);
  {
    int coarsened = 0, refined = 0;
    cache (Index) parents;
    morton (TREE_ROOT_INDEX)
      if (tree_has_children (node))
        cache_push (parents, (Index) {IJKL});
    foreach (p in parents)              /* children before their parent */
      coarsened += tree_coarsen (p);
    double d = kept_against_fresh (v), du = kept_against_fresh (u);
    printf ("  %d families coarsened: kept against fresh %g, %g\n", coarsened, d, du);
    CHECK (coarsened > 0 && d < 10 * tol && du < 10 * tol);

    set (v, wavy);
    ghosts_update_scalars (0);
    coarsened = 0;
    srand (11);
    for (int it = 0; it < 60; it++) {
      double r [3] = { rand () / (double) RAND_MAX, rand () / (double) RAND_MAX,
                       rand () / (double) RAND_MAX };
      Index q = leaf_at (r);
      if (it % 3 == 2 && q.level > 1) {
        q.level--;
        q.i >>= 1;
#if dimension > 1
        q.j >>= 1;
#endif
#if dimension > 2
        q.k >>= 1;
#endif
        coarsened += tree_coarsen (q);
      }
      else
        refined += tree_refine (q);
    }
    d = kept_against_fresh (v), du = kept_against_fresh (u);
    printf ("  %d refined, %d coarsened: kept against fresh %g, %g\n", refined, coarsened, d, du);
    CHECK (refined > 0 && d < 10 * tol && du < 10 * tol);
    /* and a linear field stays exact through it all */
    CHECK (ghost_error (u, linear, 0) < 10 * tol);
  }

  /* 8. a freed scalar gives its id back with the default operations */
  {
    int id;
    {
      scalar t;
      id = t;
      scalar_set_ops (t, (ScalarOps) { .prolongation = prolong_injection });
      ghost_ops (t, boundary.right (bc_odd));
      CHECK (scalar_get_ops (t).prolongation == prolong_injection
             && scalar_get_boundary (t, GHOST_RIGHT) == bc_odd);
    }
    scalar t2 = scalar_new ();
    CHECK (t2 == id && scalar_get_ops (t2).prolongation == prolong_linear);
    for (int f = 0; f < GHOST_FACES; f++)
      CHECK (scalar_get_boundary (t2, f) == bc_neumann);
    scalar_free (t2);
  }

  /* 9. blocks freed and reused: the hot list, with the same accounting */
  {
    long before = grid_scalar_blocks_attached ();
    { scalar t; (void) t; }
    CHECK (grid_scalar_blocks_attached () == before && grid_scalar_check ());
    grid_scalar_trim (0);
    CHECK (grid_scalar_check ());
  }
  /* 9b. the range a prolongation is checked for, once per call: a region
   * whose ends map into the coarse block, and one cell past either end not;
   * an empty region always */
  {
    BlockShape sh = tree_block_shape ();
    int m = sh.m;
    Prolongation pr = { .shape = sh };
    for (int a = 0; a < 3; a++)
      pr.off [a] = 0, pr.lo [a] = 0, pr.hi [a] = a < dimension ? m : 1;
    pr.off [0] = m;                        /* the upper half: 2 m - 1 last */
    CHECK (prolongation_in_range (&pr));
    pr.lo [0] = -m, pr.hi [0] = m;         /* down to X = 0 */
    CHECK (prolongation_in_range (&pr));
    pr.lo [0] = -m - 1;                    /* X = -1 */
    CHECK (!prolongation_in_range (&pr));
    pr.lo [0] = 0, pr.hi [0] = m + 1;      /* X = 2 m */
    CHECK (!prolongation_in_range (&pr));
    pr.lo [0] = 5, pr.hi [0] = 5, pr.off [0] = 100;   /* empty */
    CHECK (prolongation_in_range (&pr));
  }

  /* 10. a coupled prolongation: the depth hl from the surface hl + zbl,
   * the bed zbl made first though hl's id is lower, at level jumps, walls
   * and in refined leaves; zbl stale makes hl stale; a cycle, a wrong
   * count, or scalar_set_ops () of it are refused, binding nothing */
  {
    const double surface = 7;           /* not 'level': foreach_leaf's */
    ghost_ops (hl, prolongation (prolong_lake, zbl));
    CHECK ((ghost_dependents [zbl][hl / 64] >> (hl % 64)) & 1);
    {
      int order [NSCALAR_MAX], n = ghost_order_all (order), ph = -1, pz = -1;
      for (int q = 0; q < n; q++)
        ph = order [q] == hl ? q : ph, pz = order [q] == zbl ? q : pz;
      CHECK (pz >= 0 && ph > pz);           /* the bed first, whatever ids */
      fprintf (stderr, "  (hl %d, zbl %d: %s order of ids)\n", (int) hl,
               (int) zbl, hl < zbl ? "against the" : "in the");
    }
    set (zbl, linear);
    BlockShape shape = tree_block_shape ();
    foreach_leaf () {
      Real * h = block_data (sblocks, hl);
      const Real * b = block_data (sblocks, zbl);
      block_cells (shape, 0)
        h [c] = (Real) (surface - b [c]);
    }
    ghosts_update (hl, zbl);              /* hl given first: zbl done first */
    CHECK (surface_error (hl, zbl, surface) < tol * surface);
    double r [3] = { 0.62, 0.41, 0.33 };
    CHECK (tree_refine (leaf_at (r)) > 0);
    CHECK (surface_error (hl, zbl, surface) < tol * surface);
    ghosts_update (hl, zbl);
    CHECK (!ghosts_is_dirty (hl) && !ghosts_is_dirty (zbl));
    ghosts_mark_dirty (zbl);
    CHECK (ghosts_is_dirty (hl));
    ghosts_update (hl);                   /* zbl stale: updated with it */
    CHECK (!ghosts_is_dirty (zbl));
    /* the wavelet by each scalar's own prolongation: a lake at rest over a
     * wavy bed is none for hl, though the bed is not linear, and the
     * blocks are as they were after */
    set (zbl, wavy);
    foreach_leaf () {
      Real * h = block_data (sblocks, hl);
      const Real * b = block_data (sblocks, zbl);
      block_cells (shape, 0)
        h [c] = (Real) (surface - b [c]);
    }
    ghosts_update_levels (hl, zbl);
    {
      double wh = 0, wz = 0;
      foreach_leaf ()
        if (level > 0)
          wh = fmax (wh, block_wavelet ((Index) {IJKL}, hl)),
          wz = fmax (wz, block_wavelet ((Index) {IJKL}, zbl));
      fprintf (stderr, "  (wavelets: lake %g, its bed %g)\n", wh, wz);
      CHECK (wh < tol * surface && wz > 1e-3);    /* the bed's: not small */
      CHECK (surface_error (hl, zbl, surface) < tol * surface);
    }
    ghost_ops_abort = 0;
    fprintf (stderr, "  (three refusals expected:)\n");
    CHECK (ghost_ops_set (zbl, GHOST_PROLONGATION, GHOST_ALL_FACES, prolong_lake,
                          (const scalar []) { hl }, 1) == 1);
    CHECK (scalar_get_ops (zbl).prolongation == prolong_linear);
    CHECK (ghost_ops_set (u, GHOST_PROLONGATION, GHOST_ALL_FACES, prolong_lake, NULL, 0) == 1);
    scalar_set_ops (u, (ScalarOps) { .prolongation = prolong_lake });
    CHECK (scalar_get_ops (u).prolongation == prolong_linear);
    ghost_ops_abort = 1;
    scalar_set_ops (hl, (ScalarOps) {0});
    CHECK (!((ghost_dependents [zbl][hl / 64] >> (hl % 64)) & 1));
  }

  printf ("  %d leaves, %ld ghost cells at level jumps\n", tree_get_nleaves (), jumps);
  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
