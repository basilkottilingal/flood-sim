/* cwrap's foreach () and s[...] on the leaf blocks, against the same sums
 * written by hand with foreach_leaf, block_cells and cell_at */
#include <stdio.h>
#include <math.h>
#include "block-loops.h"

static int checks, fails;
#define CHECK(c) do { checks++; if (!(c)) { fails++; \
  printf ("FAIL line %d: %s\n", __LINE__, #c); } } while (0)

static long ipow (long x, int e) { long r = 1; while (e--) r *= x; return r; }

scalar s, t, r, q;
#if dimension > 1
vector u;
#endif
scalar a0, a1;
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

  /* offsets that are expressions, and +=, ++ on a cell */
  int one = 1;
  foreach () {
    r[] = t[one - 2];
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

  tree_destroy ();
  printf ("dimension %d: %d checks, %d failed\n", dimension, checks, fails);
  return fails != 0;
}
