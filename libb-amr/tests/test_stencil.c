/* cwrap's foreach () and s[...] on the leaf blocks, against the same sums
 * written by hand with foreach_leaf, block_cells and cell_at */
#include <stdio.h>
#include <math.h>
#include "block-ghosts.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

static long ipow (long x, int e) { long r = 1; while (e--) r *= x; return r; }

scalar s, t, r, q;
#if dimension > 1
vector u;
#endif
scalar a0, a1, gp, gq, gr;
scalar list [2];          /* an array: no ids of its own, a0 and a1 below */

/* t over every cell of every leaf, ghosts included: a value per cell */
static void fill_t (void)
{
  BlockShape shape = tree_block_shape ();
  foreach_leaf () {
    Real * a = block_data (sblocks, t);
    block_cells (shape, shape.g)
      a [c] = (Real) sin (0.37 * c + level);
  }
}

/* the largest difference between two scalars over the interior cells */
static double diff (scalar a, scalar b)
{
  BlockShape shape = tree_block_shape ();
  double e = 0;
  foreach_leaf () {
    const Real * x = block_data (sblocks, a), * y = block_data (sblocks, b);
    block_cells (shape, 0)
      e = fmax (e, fabs ((double) x [c] - (double) y [c]));
  }
  return e;
}

/* a scalar given as a parameter */
static void set_all (scalar a, Real v)
{
  foreach ()
    a[] = v;
}

int main (void)
{
  tree (4, 3);
  BlockShape shape = tree_block_shape ();
  long cells = tree_get_nleaves () * ipow (shape.m, dimension);
  list [0] = a0;
  list [1] = a1;

  /* 1. every interior cell once, and a sum over them */
  long n = 0;
  foreach () {
    s[] = 2;
    n++;
  }
  CHECK (n == cells);
  double total = 0;
  foreach ()
    total += s[];
  CHECK (total == 2. * cells);

  /* the ghosts are not written */
  long ghosts = 0;
  foreach_leaf () {
    const Real * a = block_data (sblocks, s);
    block_cells (shape, shape.g)
      ghosts += a [c] == 0;
  }
  CHECK (ghosts == tree_get_nleaves () * (ipow (shape.m + 2 * shape.g, dimension)
                                          - ipow (shape.m, dimension)));

  /* 2. a stencil, offsets left out being 0, against cell_at */
  fill_t ();
  foreach ()
#if dimension == 1
    r[] = t[1] + t[-1] - 2 * t[] + t[2];
#elif dimension == 2
    r[] = t[1] + t[-1,0] + t[0,1] + t[0,-1] - 4 * t[] + t[1,-1];
#else
    r[] = t[1] + t[-1,0] + t[0,1] + t[0,-1,0] + t[0,0,1] + t[0,0,-1]
          - 6 * t[] + t[1,-1,2];
#endif
  foreach_leaf () {
    const Real * a = block_data (sblocks, t);
    Real * b = block_data (sblocks, q);
    block_cells (shape, 0)
#if dimension == 1
      b [c] = cell_at (a, 1) + cell_at (a, -1) - 2 * cell_at (a, 0) + cell_at (a, 2);
#elif dimension == 2
      b [c] = cell_at (a, 1, 0) + cell_at (a, -1, 0) + cell_at (a, 0, 1)
              + cell_at (a, 0, -1) - 4 * cell_at (a, 0, 0) + cell_at (a, 1, -1);
#else
      b [c] = cell_at (a, 1, 0, 0) + cell_at (a, -1, 0, 0) + cell_at (a, 0, 1, 0)
              + cell_at (a, 0, -1, 0) + cell_at (a, 0, 0, 1) + cell_at (a, 0, 0, -1)
              - 6 * cell_at (a, 0, 0, 0) + cell_at (a, 1, -1, 2);
#endif
  }
  CHECK (diff (r, q) == 0);

  /* offsets that are constant expressions, and +=, ++ on a cell */
  foreach () {
    r[] = t[1 - 2];
    r[] += 1;
    r[]++;
  }
  foreach_leaf () {
    const Real * a = block_data (sblocks, t);
    Real * b = block_data (sblocks, q);
    block_cells (shape, 0)
    {
      Real v = a [c - 1];
      v += 1;
      v++;
      b [c] = v;
    }
  }
  CHECK (diff (r, q) == 0);

  /* 3. a vector's components, an array of scalars, a pointer to one */
  set_all (s, 2);
#if dimension > 1
  foreach () {
    u.x[] = t[];
    u.y[] = s[];
  }
  CHECK (diff (u.x, t) == 0 && diff (u.y, s) == 0);
#endif
  foreach ()
    for (int k = 0; k < 2; k++)       /* list[k]: found again at each cell */
      list[k][] = k + s[];
  scalar * p = &list [1];
  foreach ()
    r[] = (*p)[] - 1;
  CHECK (diff (r, list [0]) == 0);
  foreach ()
    r[] = list[0][] - p[0][] + 1;
  set_all (q, 0);
  CHECK (diff (r, q) == 0);

  /* 4. the body sees its own names, not the loops': c, shape, i, level, ... */
  {
    double c = 0.25, level = 0.5;
    int i = 3, node = 4, sblocks = 5;
    foreach ()
      s[] = c + level + i + node + sblocks;
    set_all (q, (Real) (0.25 + 0.5 + 3 + 4 + 5));
    CHECK (diff (s, q) == 0);
  }
  /* and a name of its own, declared in the body, may be c */
  foreach () {
    Real c = t[];
    s[] = c;
  }
  CHECK (diff (s, t) == 0);

  /* 5. x, y, z and Delta: the centre and size of the cell, on a moved and
   * resized domain, with two levels of leaves; y and z are 0 beyond the
   * dimension, an outer x is hidden, the body's own x is its own */
  {
    Index first = { 0 };
    int got = 0;
    foreach_leaf ()
      if (!got)
        first = (Index) {IJKL}, got = 1;
    CHECK (tree_refine (first) > 0);
  }
  tree_shift_resize (-1, 2, 0.5, 4);
  double x = 99;
  foreach ()
    s[] = x + 10 * y + 100 * z + 1000 * Delta;
  foreach_leaf () {
    Real * b = block_data (sblocks, q);
    double D = 4. / ((double) shape.m * (double) (1L << level));
    block_cells (shape, 0) {
      double X = -1 + ((double) i * shape.m + (ci - shape.ox) + 0.5) * D, Y = 0, Z = 0;
#if dimension > 1
      Y = 2 + ((double) j * shape.m + (cj - shape.g) + 0.5) * D;
#endif
#if dimension > 2
      Z = 0.5 + ((double) k * shape.m + (ck - shape.g) + 0.5) * D;
#endif
      b [c] = (Real) (X + 10 * Y + 100 * Z + 1000 * D);
    }
  }
  CHECK (diff (s, q) < (sizeof (Real) == 4 ? 1e-3 : 1e-10));
  CHECK (x == 99);
  double dmin = 1e30, dmax = 0;          /* the refined leaf: Delta / 2 */
  foreach () {
    dmin = fmin (dmin, Delta);
    dmax = fmax (dmax, Delta);
  }
  CHECK (dmax == 2 * dmin);
  foreach () {
    double x = 7;
    s[] = x;
  }
  set_all (q, 7);
  CHECK (diff (s, q) == 0);

  /* 6. the ghosts of a scalar written, updated before a foreach () reads
   * it at an offset, as ghosts_update () by hand would, level jump and
   * walls included; and only then */
  foreach ()
    gp[] = x * x + 0.5 * y - z;
  CHECK (ghosts_is_dirty (gp));
  foreach ()
    gr[] = gp[1] - gp[-1];
  CHECK (!ghosts_is_dirty (gp) && ghosts_is_dirty (gr));
  foreach ()
    gq[] = x * x + 0.5 * y - z;
  ghosts_update (gq);
  CHECK (!ghosts_is_dirty (gq));
  foreach_leaf () {
    const Real * a = block_data (sblocks, gq);
    Real * b = block_data (sblocks, q);
    block_cells (shape, 0)
      b [c] = a [c + 1] - a [c - 1];
  }
  CHECK (diff (gr, q) == 0);

  /* 7. foreach_neighbor (1): the box around each cell, read after its
   * ghosts are updated (gp written just before), against the same sum by
   * hand; x of the neighbours, whose mean is the cell's */
  foreach ()
    gp[] = x * y + z - x;
  CHECK (ghosts_is_dirty (gp));
  foreach () {
    double sum = 0;
    foreach_neighbor (1)
      sum += gp[];
    gr[] = sum;
  }
  CHECK (!ghosts_is_dirty (gp));
  foreach_leaf () {
    const Real * a = block_data (sblocks, gp);
    Real * b = block_data (sblocks, q);
    int lo [3] = { -1, dimension > 1 ? -1 : 0, dimension > 2 ? -1 : 0 };
    int hi [3] = { 1, dimension > 1 ? 1 : 0, dimension > 2 ? 1 : 0 };
    block_cells (shape, 0) {
      double sum = 0;
      for (int dk = lo [2]; dk <= hi [2]; dk++)
      for (int dj = lo [1]; dj <= hi [1]; dj++)
      for (int di = lo [0]; di <= hi [0]; di++)
        sum += a [c + di + dj * shape.sy + dk * shape.sz];
      b [c] = sum;
    }
  }
  CHECK (diff (gr, q) == 0);
  foreach () {
    double xs = 0;
    int n = 0;
    foreach_neighbor (1) {
      xs += x;
      n++;
    }
    gr[] = xs / n - x;
  }
  set_all (q, 0);
  CHECK (diff (gr, q) < 1e-12);

  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
