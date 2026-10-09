#include "flow-utils.h"
#include "pgm.h"
#include "grid.h"

#include <stdio.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

# if defined(__GNUC__) || defined(__clang__)
#   if defined(__has_builtin)
#     if __has_builtin(__builtin_ctz)
#       define HAVE_CTZ
#     endif
#   endif
# endif

int pop_bit (uint8_t * v)
{
  if (*v == 0)
    return -1;

  #ifdef HAVE_CTZ
  int idx = __builtin_ctz ((unsigned int)*v);
  #else
  static const int8_t debruijn_table[8] =
    {
      0, 1, 6, 2, 7, 5, 4, 3
    };
  uint8_t isolated = (uint8_t) (*v & (-(int8_t)*v));
  int idx = debruijn_table [(uint8_t)(isolated * 0x1DU) >> 5];
  #endif

  *v = (uint8_t)(*v & (*v - 1));
  return idx;
}

/*
.. validates
.. (a) dir and invDir are complementary
.. (b) dir has only one bit
.. (c) and there is no cycle in the network
*/
void DAG_validity (uint8_t ** dir, uint8_t ** invDir)
{
  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height, hw = w * h;

  for (int y=-1; y<=h; ++y)
    for(int x=-1; x<=w; ++x)
      if (dir [y][x])
      {
        assert (! (y<0 || x<0 ||x>=w || y>=h));
        uint8_t code = dir [y][x];
        int nbr = pop_bit (&code);
        assert (!code);
        int inv = (nbr+4)%8;
        int xnbr = x + flow_neighbor [nbr].x;
        int ynbr = y + flow_neighbor [nbr].y;
        assert (invDir [ynbr][xnbr] & ((uint8_t) 1 << inv));
      }

  for (int y=0; y<h; ++y)
    for(int x=0; x<w; ++x)
      if (invDir [y][x])
      {
        uint8_t code = invDir [y][x];
        while (code)
        {
          int nbr = pop_bit (&code);
          int inv = (nbr+4)%8;
          int xnbr = x + flow_neighbor [nbr].x;
          int ynbr = y + flow_neighbor [nbr].y;
          assert (dir [ynbr][xnbr] & ((uint8_t) 1 << inv));
        }
      }

  /* making sure that there is no cycle in the graph */
  uint8_t ** visited = grid (uint8_t, w, h);
  assert (visited != NULL);
  for (int y=0; y<h; ++y)
    for(int x=0; x<w; ++x)
      if (dir [y][x])
      {
        int depth = 0;
        int xb = x;
        int yb = y;
        while ( !visited [y][x] && !(y<0 || x < 0 || x>= w || y>=h) && dir [y][x] )
        {
          visited [y][x] = 1u;
          uint8_t code = dir [y][x];
          int nbr = pop_bit ( &code );
          x += flow_neighbor [nbr].x, y += flow_neighbor [nbr].y;
          assert (! (x == xb && y == yb) );
          assert(depth++ < hw);
        }
        x = xb, y = yb;
      }
  grid_free (visited);
  printf ("DAG looks fine\n");
  fflush (stdout);
}

/* Forest of Directed Acyclic Graph corresponding to D8 flow network */
void D8 (double ** raster, uint8_t ** dir, uint8_t ** invDir)
{

  assert (raster && dir && invDir);
  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height;
  memset (& dir [-2][-2], 0, (w+4) * (h+4) * sizeof (uint8_t));
  memset (& invDir [-2][-2], 0, (w+4) * (h+4) * sizeof (uint8_t));

  double a = 1., b = sqrt (2.);
  const double delta [8] = {a, b, a, b, a, b, a, b};
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      double max = 0., val = raster [y][x];
      int nbr = -1;
      if ( isnan (val) )
        continue;
      for (int c=0; c<8; ++c)
      {
        double nbrVal = raster [y + flow_neighbor[c].y][x + flow_neighbor[c].x];
        if (isnan (nbrVal))
          continue;
        double downslope = (val - nbrVal) / delta [c];
        if (downslope > max)
          nbr = c, max = downslope;
      }
      if (nbr >= 0)
      {
        dir [y][x] = (uint8_t) 1 << nbr;
        int xnbr = x + flow_neighbor [nbr].x;
        int ynbr = y + flow_neighbor [nbr].y;
        /* DAG where every direction is inverted */
        invDir [ynbr][xnbr] |= (uint8_t) 1 << ((nbr+4)%8);
      }
    }
}

/* may not be equivalent to D infty algorithm in the literature*/
void DInfty (double ** raster, uint8_t ** dir, uint8_t ** invDir)
{

  assert (raster && dir && invDir);
  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height;
  memset (& dir [-2][-2], 0, (w+4) * (h+4) * sizeof (uint8_t));
  memset (& invDir [-2][-2], 0, (w+4) * (h+4) * sizeof (uint8_t));

  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      double val = raster [y][x];
      if ( isnan (val) )
        continue;
      for (int c=0; c<8; ++c)
      {
        int xnbr = x + flow_neighbor [c].x;
        int ynbr = y + flow_neighbor [c].y;
        double nbrVal = raster [ynbr][xnbr];
        if (isnan (nbrVal))
          continue;
        if (nbrVal < val) /* positive downslope */
        {
          dir [y][x]    = (uint8_t) 1 << c;
          invDir [y][x] = (uint8_t) 1 << ((c+4)%8);
        }
      }
    }
}

int drain_validity (double ** elevation, int ** tag, int ntags)
{

  assert (elevation && tag && ntags > 0);
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** dir     = grid (uint8_t, w, h);
  uint8_t ** invDir  = grid (uint8_t, w, h);
  assert (dir && invDir);
  D8 (elevation, dir, invDir);

  int * isdrain = malloc (ntags * sizeof (int));
  assert (isdrain);
  memset (isdrain, 0, ntags * sizeof (int));

  int err1 = 0, err2 = 0;

  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      if (!tag [y][x])
      {
        if (dir [y][x] || (x == 0 || x == w-1 || y == 0 || y == h-1))
          continue;
        err1++;
      }
      else if (dir [y][x] || (x == 0 || x == w-1 || y == 0 || y == h-1))
        isdrain [tag [y][x]] = 1;
    }

  if (err1 || err2)
  {
    uint8_t ** stagnant = grid (uint8_t, w, h);
    assert (stagnant);
    for (int y=0; y<h; ++y)
      for (int x=0; x<w; ++x)
      {
        if (!tag [y][x])
        {
          if (dir [y][x] || (x == 0 || x == w-1 || y == 0 || y == h-1))
            continue;
          stagnant [y][x] = 255u;
        }
        else if (!isdrain [tag [y][x]])
          stagnant [y][x] = 255u;
      }
    pgm (stagnant, "stagnant.pgm");
    grid_free (stagnant);
  }

  free (isdrain);
  grid_free (dir);
  grid_free (invDir);

  return !(err1 || err2);
}
