#include "geotiff.h"
#include "filemap.h"
#include "flow.h"
#include "min-heap.h"
#include "pgm.h"
#include "grid.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define error(e) filemap_close_all(e)

typedef enum
{
  SUCCESS = 0,
  ERR_FILE_ACCESS,
  ERR_NOT_IMPLEMENTED,
  ERR_MALLOC,
} ERR;

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

int flow_network (DEM dem, FLOW_DRAIN type, FlowNetwork * network)
{
  if (dem.raster == NULL)
    return -1;

  int w = dem.w, h = dem.h;
  assert (w > 0 && h > 0);
  uint8_t ** dir = grid (uint8_t, w, h);
  float ** raster = dem.raster;

  *network = (FlowNetwork)
    {
      .dir    = dir,
      .accumulation = NULL,
      .w      = w,
      .h      = h,
      .type   = type,
    };

  float a = 1., b = 1. / sqrt (2.);
  const float delta [8] = {a, b, a, b, a, b, a, b};

  switch (type)
  {
    case FLOW_D8  :
      for (int y=0; y<h; ++y)
        for (int x=0; x<w; ++x)
        {
          float max = 0.f, val = raster [y][x];
          uint8_t code = 0;
          if ( isnan (val) )
            continue;
          for (int c=0; c<8; ++c)
          {
            float nbrVal = raster [y + flow_neighbor[c].y][x + flow_neighbor[c].x];
            if (isnan (nbrVal))
              continue;
            float downslope = (val - nbrVal) / delta [c];
            if (downslope > max)
              code = (uint8_t) 1 << c, max = downslope;
          }
          dir [y][x] = code;
        }
      return 0;

    case FLOW_DINFTY :
    case FLOW_MFD :
      error ("flow drain type not implemented");
      break;

    default :
      error ("unknown flow drain type"); 
  }

  return -1;
}

void flow_network_free (FlowNetwork network)
{
  assert (network.dir != NULL);
  grid_free (network.dir);
  if (network.accumulation == NULL)
    return;
  grid_free (network.accumulation);
}

int flow_network_grayscale (FlowNetwork network, const char * out, int width, int height)
{
  if (network.dir == NULL || out == NULL || out [0] == '\0' || width < 0 || height < 0)
    return -1;

  width  = width  > network.w ? network.w > 1024 ? 1024 : network.w : width;
  height = height > network.h ? network.h > 1024 ? 1024 : network.h : height;

  FILE *fp = fopen (out, "wb");
  if (fp == NULL)
    return -1;

  /* PGM header */
  fprintf (fp, "P5\n%d %d\n255\n", width, height);
  for (int y = 0; y < height; y++)
    fwrite (network.dir [y], 1, width, fp);
  fclose(fp);

  return 0;  
}

int flow_accumulation_grayscale (FlowNetwork network, const char * out, int width, int height)
{
  width  = width  > network.w ? network.w > 1024 ? 1024 : network.w : width;
  height = height > network.h ? network.h > 1024 ? 1024 : network.h : height;

  return pgm_grayscale (network.accumulation, out, width, height);
}

int flow_accumulation (FlowNetwork * network, DEM dem, MinPQ * sink)
{
  if (network->dir == NULL || network->accumulation != NULL)
    return -1;

  int w = network->w, h = network->h, hw = w * h;
  float ** acc = network->accumulation = grid (float, w, h);
  float ** raster = dem.raster;
  uint8_t ** in_degree = grid (uint8_t, w, h);
  uint8_t ** dir = network->dir;
  if (sink && pq_create (sink))
    error ("flow_accumulation () : pq_create () : malloc");

  /* find the number of neighbors who are source to (i,j) */
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      acc [y][x] = 1.; /* you first add the "rainfall" at (i,j) to the flow[i][j]*/
      uint8_t code = dir [y][x];
      if ( code == 0 )
        continue;
      uint8_t nbr = pop_bit (&code); /* (i, j) is a source to nbr */
      int xnbr = x + flow_neighbor[nbr].x, ynbr = y + flow_neighbor[nbr].y;
      in_degree [ynbr][xnbr] |= (uint8_t) 1 << nbr;
    }

  /* first in first out queue */
  struct index { int x, y; } * queue = malloc (hw * sizeof (struct index));
  if (queue == NULL)
    error ("flow_accumulation () : malloc () failed");
  int push_at = 0, pop_at = 0;
  int processed = 0;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      if (isnan (raster [y][x]))
      {
        processed ++;
        continue;
      }
      if (in_degree [y][x] == 0)
        queue [push_at++] = (struct index) {x,y};
    }

  /* Kahn's topological propogation */
  while (push_at != pop_at)
  {
    processed ++;

    /* pop */
    struct index Idx = queue [pop_at]; pop_at = (pop_at+1) % hw;
    int x = Idx.x, y = Idx.y;
    uint8_t code = dir [y][x];
    if (code == 0)
      continue;
    
    uint8_t nbr = pop_bit (&code); /* (i, j) is a source to nbr */
    int xnbr = x + flow_neighbor[nbr].x, ynbr = y + flow_neighbor[nbr].y;

    if (xnbr<0 || xnbr>=w || ynbr<0 || ynbr>=h)
    {
      /* add this point to sink */
      if (sink && pq_push (sink, raster [ynbr][xnbr], xnbr, ynbr))
        error ("flow_accumulation () : pq_push () : realloc");
      continue;
    }

    /* update downstream accumulation*/
    acc [ynbr][xnbr] += acc [y][x];

    assert ( in_degree [ynbr][xnbr] & ((uint8_t) 1 << nbr) );
    /* push */
    if ( (in_degree [ynbr][xnbr] &= ~((uint8_t) 1 << nbr)) == 0 )
      queue [push_at] = (struct index) {xnbr, ynbr}, push_at = (push_at + 1) % hw;
    else if (sink && pq_push (sink, raster [ynbr][xnbr], xnbr, ynbr))
      error ("flow_accumulation () : pq_push () : realloc");
  }

  if (processed != hw)
    fprintf (stderr, "Kahn's topological propagation inconsistency");

  free (queue);
  grid_free (in_degree);
  return 0;
}

/*
.. depression filling/priority flood fill algorithm
*/
int flow_remove_pits (DEM * dem)
{
  int w = dem->w, h = dem->h;
  float ** elevation = dem->raster;
  uint8_t ** visited = grid (uint8_t, w, h);
  
  MinPQ pq;
  if (pq_create (&pq))
    return ERR_MALLOC;

  #define push(X,Y)                                                    \
    visited [Y][X] = 1;                                                \
    if (!isnan (elevation [Y][X]))                                     \
      if (pq_push (&pq, elevation [Y][X], X, Y))                       \
        return ERR_MALLOC;

  /* set points outside the box as visited */
  memset (&visited [-1][-1], -1, (w+2) * sizeof (uint8_t));
  memset (&visited [ h][-1], -1, (w+2) * sizeof (uint8_t));
  for (int y=0; y<h; ++y)
    visited [y][-1] = visited [y][w] = 1;

  for (int y=0; y<h; ++y)
  {
    push (0,y); push (w-1,y);
  }
  for (int x=0; x<w; ++x)
  {
    push (x,0); push (x,h-1);
  }

  #undef push

  float elev; int x, y;
  while (pq_pop (&pq, &elev, &x, &y))
  {
    for (int c=0; c<8; ++c)
    {
      int xnbr = x + flow_neighbor[c].x;
      int ynbr = y + flow_neighbor[c].y;
      if (visited [ynbr][xnbr])
        continue;
      visited [ynbr][xnbr] = 1;
      if (isnan (elevation [ynbr][xnbr]))
        continue;
      if (elev > elevation [ynbr][xnbr])
        elevation [ynbr][xnbr] = elev;
      if (pq_push (&pq, elevation [ynbr][xnbr], xnbr, ynbr))
        return ERR_MALLOC;
    }
  }

  grid_free (visited);
  pq_free   (&pq);

  return 0;
}


void flow_network_error ( int type )
{
  switch (type)
  {
    default :
      fprintf (stderr, "unknonw flow network error\n");
  }
}

#undef error
