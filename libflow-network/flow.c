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

#define error(e) filemap_close_all(e)

typedef enum
{
  SUCCESS = 0,
  ERR_FILE_ACCESS,
  ERR_NOT_IMPLEMENTED,
  ERR_MALLOC,
  ERR_DATA_MISSING,
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

static void D8 (double ** raster, int w, int h, uint8_t ** dir, uint8_t ** invDir)
{
  /* Forest of Directed Acyclic Graph corresponding to D8 flow network */

  double a = 1., b = sqrt (2.);
  const double delta [8] = {a, b, a, b, a, b, a, b};

  //for (int y=-1; y<h+1; ++y)
  //  for (int x=-1; x<w+1; ++x)
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

int flow_network (DEM dem, FLOW_DRAIN type, FlowNetwork * network)
{
  if (dem.raster == NULL)
    return ERR_DATA_MISSING;

  int w = dem.w, h = dem.h;
  assert (w > 0 && h > 0);
  uint8_t ** dir = grid (uint8_t, w, h);
  uint8_t ** invDir = grid (uint8_t, w, h);
  if (dir == NULL || invDir == NULL)
    return ERR_MALLOC;
  double ** raster = dem.raster;

  *network = (FlowNetwork)
    {
      .dir    = dir,
      .invDir = invDir,
      .accumulation = NULL,
      .w      = w,
      .h      = h,
      .type   = type,
    };

  switch (type)
  {
    case FLOW_D8  :
      D8 (raster, w, h, dir, invDir);
      return 0;

    case FLOW_DINFTY :
    case FLOW_MFD :
      error ("flow drain type not implemented");
      break;

    default :
      error ("unknown flow drain type"); 
  }

  assert (0);
  return -1;
}

void flow_network_free (FlowNetwork network)
{
  assert (network.dir != NULL && network.invDir != NULL);
  grid_free (network.dir);
  grid_free (network.invDir);
  if (network.accumulation == NULL)
    return;
  grid_free (network.accumulation);
}

int flow_accumulation (FlowNetwork * network, DEM dem, MinPQ * sink)
{
  if (network->dir == NULL || network->accumulation != NULL || dem.raster ==  NULL)
    return ERR_DATA_MISSING;

  int w = network->w, h = network->h, hw = w * h;
  double ** acc = network->accumulation = grid (double, w, h);
  double ** raster = dem.raster;
  uint8_t ** in_degree = grid (uint8_t, w, h);
  uint8_t ** dir = network->dir;

  if (acc == NULL || in_degree == NULL)
    return ERR_MALLOC;

  if (sink && pq_create (sink))
    error ("flow_accumulation () : pq_create () : malloc");

  /* every cell is a source to itself */
  memset (&acc [-2][-2], 1, (w+4)*(h*4) *sizeof (uint8_t));

  /* neighbors who are source to (x,y) are codified as '1' bits of in_degree*/
  grid_copy (uint8_t, network->invDir, in_degree, w, h);

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

    #define collect_sink(X,Y) \
      if (sink && pq_push (sink, raster [Y][X], X, Y)) \
        return ERR_MALLOC

    /* pop */
    struct index Idx = queue [pop_at]; pop_at = (pop_at+1) % hw;
    int x = Idx.x, y = Idx.y;
    uint8_t code = dir [y][x];
    if (code == 0)
    {
      collect_sink (x, y);
      continue;
    }
    
    uint8_t nbr = pop_bit (&code), inv = (nbr+4)%8;
    int xnbr = x + flow_neighbor[nbr].x, ynbr = y + flow_neighbor[nbr].y;

    if (xnbr<0 || xnbr>=w || ynbr<0 || ynbr>=h)
    {
      collect_sink (xnbr, ynbr);
      continue;
    }

    /* update downstream accumulation*/
    acc [ynbr][xnbr] += acc [y][x];

    assert ( in_degree [ynbr][xnbr] & ((uint8_t) 1 << inv) );
    if ( ! (in_degree [ynbr][xnbr] &= ~((uint8_t) 1 << inv))  )
      queue [push_at] = (struct index) {xnbr, ynbr}, push_at = (push_at + 1) % hw;
    else
      collect_sink (xnbr, ynbr);

    #undef collect_sink
  }
  assert ("Kahn's topological consistency" && processed == hw );

  free (queue);
  grid_free (in_degree);
  return 0;
}

/*
.. depression filling algorithm. also called as priority flood fill algorithm
*/
int flow_remove_pits (DEM * dem)
{
  int w = dem->w, h = dem->h;
  assert (w > 0 && h > 0);
  double ** elevation = dem->raster;
  if (elevation == NULL)
    return ERR_DATA_MISSING;
  uint8_t ** visited = grid (uint8_t, w, h);
  if (visited == NULL)
    return ERR_MALLOC;
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

  /* boundary points are pushed to the priority queue */
  for (int y=0; y<h; ++y)
  {
    push (0,y); push (w-1,y);
  }
  for (int x=0; x<w; ++x)
  {
    push (x,0); push (x,h-1);
  }

  #undef push

  double elev; int x, y;
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

static void DAG_validity (uint8_t ** dir, uint8_t ** invDir)
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
  printf ("DAG looks fine"); fflush (stdout);
}

/* creates a drain network in flattened reservoirs and filled pits */
static int fix_flattened_patches (uint8_t ** dir, uint8_t ** invDir, uint8_t ** type)
{

  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height, hw = w * h;
  uint8_t ** bit_stack = grid (uint8_t, w, h);
  uint8_t ** visited = grid (uint8_t, w, h);
  if (bit_stack == NULL || visited == NULL)
    return ERR_MALLOC;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      /* if there are multiple drains to a flattened patch, we will drain all
      .. those flattened pixels into (x,y). Actually there should be some
      .. tie-breaking algorithm */
      if ( (type [y][x] & FN_FLAT) && dir [y][x] && !visited [y][x])
      {
        /* now connect all the neighboring reservoir points. using tree traversal */
        int depth = 0;
        do
        {

          do
          {
            /*
            .. all the neighbors at same elevation which doesn't have a drain 'dir'
            .. now will drain to (x,y);
            */
            if (!visited [y][x])
            {
              visited [y][x] = 1u;
              assert (!bit_stack[y][x]);
              uint8_t code = dir [y][x];
              int nbr = pop_bit (&code);
              for (int c=0; c<8; ++c)
              {
                int xnbr = x+flow_neighbor[c].x;
                int ynbr = y+flow_neighbor[c].y;
                if ( c != nbr                         &&
                   (type [ynbr][xnbr] & FN_FLAT) &&
                   !visited [ynbr][xnbr]              &&
                   !dir [ynbr][xnbr]                  &&
                   !(ynbr < 0 || xnbr < 0 || xnbr >= w || ynbr >= h) )
                {
                  bit_stack [y][x] |= (uint8_t) 1 << c;
                  dir [ynbr][xnbr] = (uint8_t) 1 << ((c+4)%8);
                }
              }
              invDir [y][x] |= bit_stack [y][x];
            }

            if (!bit_stack [y][x])
              break;

            depth++;
            /* remove the popped bit so that you don't traverse it again */
            int nbr = pop_bit ( &bit_stack [y][x] );
            x += flow_neighbor [nbr].x;
            y += flow_neighbor [nbr].y;
            assert (depth < hw);

          } while (1);

          if (!depth--)
            break;

          /* don't edit dir[y][x]. Instead use a copy of it in pop_bit () */
          uint8_t code = dir [y][x]; assert (code);
          int nbr = pop_bit ( &code );
          x += flow_neighbor [nbr].x;
          y += flow_neighbor [nbr].y;
        } while (1);
      }
  #ifdef _FLOW_DEBUG_
  DAG_validity (dir, invDir);
  #endif
  grid_free (bit_stack);
  grid_free (visited);
  return 0;
}

static
uint8_t ** classify_nodes (double ** elevation, int w, int h, uint8_t ** dir, uint8_t ** invDir)
{
  uint8_t ** type = grid (uint8_t, w, h);
  uint8_t ** invDirCopy = grid (uint8_t, w, h);
  grid_copy (uint8_t, invDir, invDirCopy, w, h);
  if (type == NULL || invDirCopy == NULL)
    return NULL;
  /* fixme : change this to -1:y, -1:w */
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      if (y == 0|| x == 0 || y == h-1 || x == w-1)
        type [y][x] |= FN_BOUNDARY;

      if (elevation [y][x] <= 0.)
      {
        type [y][x] |= FN_OCEAN;
        continue; 
      }

      /* no source & no sink : most likely flattened reservoirs.
      .. This happens usually when the neighbors are also flattened
      .. to the same elevation by removing noises to be specifically identified
      .. as reservoirs. */
      if (!dir [y][x] && !invDir [y][x])
      {
        type [y][x] |= FN_FLAT;
    #if 0
        for (int dx =-1; dx <=1; ++dx)
          for (int dy =-1; dy <=1; ++dy)
          {
            //type [y+dy][x+dx] |= FN_RESERVOIR;
            if (invDir [y+dy][x+dx])
            {
    /* tree traversal to identify the catchement area */
    int depth = 0;
    int xbackup = x, ybackup = y;
    x += dx, y += dy;
    do
    {
      while ( invDirCopy [y][x] )
      {
        depth++;
        /* remove the popped bit so that you don't traverse it again */
        int nbr = pop_bit ( &invDirCopy [y][x] );
        x += flow_neighbor [nbr].x;
        y += flow_neighbor [nbr].y;
        assert (depth < hw);
      }

      type [y][x] |= FN_RESERVOIR_CATCHMENT;

      if (!depth--)
        break; 

      /* don't edit dir[y][x]. Instead use a copy of it in pop_bit () */
      uint8_t code = dir [y][x]; assert (code);
      int nbr = pop_bit ( &code );
      x += flow_neighbor [nbr].x;
      y += flow_neighbor [nbr].y;
    } while (1);
    /* end of tree traversal. iterators x and y are updated with the back up*/
    x = xbackup, y= ybackup;
            }
          }
    #endif
      }
    }
  return type;
}


int flow_routine (DEM dem)
{
  if (dem.raster == NULL || dem.w < 2 || dem.h < 2)
    return ERR_DATA_MISSING;
  int w = dem.w, h = dem.h;
  const double * const * raster = (const double * const *) dem.raster;
  double ** elevation = grid (double, w, h); /* for a copy of raster */
  if (elevation == NULL)
    return ERR_MALLOC;
  grid_copy (double, raster, elevation, w, h);
  uint8_t ** dir = grid (uint8_t, w, h);
  uint8_t ** invDir = grid (uint8_t, w, h);
  uint8_t ** color = grid (uint8_t, w, h);

  if (dir == NULL || invDir == NULL || color == NULL)
    return ERR_MALLOC;

  D8 (elevation, w, h, dir, invDir);

  uint8_t ** type = classify_nodes (elevation, w, h, dir, invDir);
  if (type == NULL)
    return ERR_MALLOC;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      color [y][x] =
        (type [y][x] & FN_OCEAN) ? 255u :
        (type [y][x] & FN_FLAT) ? 220u :
        //(type [y][x] & FN_RESERVOIR_CATCHMENT) ? 128u:
        (type [y][x] & FN_BOUNDARY) ? 32u : 0u;

  if (fix_flattened_patches (dir, invDir, type))
    return ERR_MALLOC;

  FlowNetwork network = {.w = w, .h = h, .dir = dir, .invDir = invDir};
  int err = flow_accumulation (&network, dem, NULL);
  if (err)
    return err;
  //uint8_t ** acc = network.accumulation;

//test
pgm (dir, "d8.pgm");
pgm (invDir, "invDir.pgm");
pgm (color, "reservoir-ocean.pgm");
pgm (network.accumulation, "accumulation.pgm");

  //MinPQ sinks;

  return 0;
}

void flow_network_error ( int type )
{
  if (type == SUCCESS)
    return;

  fprintf (stderr, "flow.h : ");
  switch (type)
  {
    case ERR_FILE_ACCESS : 
      fprintf (stderr, "file access error\n");
      break;
    case ERR_NOT_IMPLEMENTED :
      fprintf (stderr, "implementation error\n");
      break;
    case ERR_MALLOC :
      fprintf (stderr, "malloc/realloc failed\n");
      break;
    case ERR_DATA_MISSING :
      fprintf (stderr, "required data missing\n");
      break;
    default :
      fprintf (stderr, "unknown flow network error\n");
  }
}

#undef error
