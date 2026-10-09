#include "geotiff.h"
#include "filemap.h"
#include "flow.h"
#include "min-heap.h"
#include "pgm.h"
#include "grid.h"
#include "fifo.h"

#include <stdio.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#define error(e)         filemap_close_all(e)
#define IS_FLAT(x,y)     (!dir [y][x] && !invDir [y][x])
#define IS_OUTSIDE(x,y)  ((x)<0 || (y)<0 || (x)>=w || (y)>=h)
#define FLAT_INLET       128u
#define FLAT_OUTLET      64u
#define FLAT             16u

typedef enum
{
  SUCCESS = 0,
  ERR_FILE_ACCESS,
  ERR_NOT_IMPLEMENTED,
  ERR_MALLOC,
  ERR_DATA_MISSING,
} ERR;

void flow_network_error (int type)
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

/* depth first traversal*/
#define tree_cell_start(invDirCopy)               \
do                                                \
{                                                 \
  assert (!dir[y][x]);                            \
  int xb = x, yb = y;                             \
  int depth = 0;                                  \
  do {                                            \
    while (invDirCopy[y][x])                      \
    {                                             \
      int nbr = pop_bit (&invDirCopy [y][x]);     \
      x += flow_neighbor [nbr].x;                 \
      y += flow_neighbor [nbr].y;                 \
      assert (depth ++ < hw);                     \
    }
#define tree_cell_end()                           \
    if (!depth--)                                 \
      break;                                      \
    uint8_t code = dir [y][x];                    \
    int nbr = pop_bit (&code);                    \
    x += flow_neighbor [nbr].x;                   \
    y += flow_neighbor [nbr].y;                   \
  } while (1);                                    \
  assert (x == xb);                               \
  assert (y == yb);                               \
}                                                 \
while (0)

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

/* Forest of Directed Acyclic Graph corresponding to D8 flow network */
static void D8 (double ** raster, int w, int h, uint8_t ** dir, uint8_t ** invDir)
{
  /* warning : assumes dir and invDir are set to 0 */
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
static void DInfty (double ** raster, int w, int h, uint8_t ** dir, uint8_t ** invDir)
{
  /* warning : assumes dir and invDir are set to 0 */
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

/* create a flow network */
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
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      acc [y][x] = 1.;

  /* neighbors who are source to (x,y) are codified as '1' bits of in_degree*/
  grid_copy (uint8_t, network->invDir, in_degree, w, h);

  /* first in first out queue */
  Queue queue; if (queue_init (&queue)) return ERR_MALLOC;
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
        if (queue_push (&queue, x, y)) return ERR_MALLOC;
    }

  /* Kahn's topological propogation */
  int x, y;
  while (queue_pop (&queue, &x, &y))
  {
    processed ++;

    #define collect_sink(X,Y) \
      if (sink && pq_push (sink, raster [Y][X], X, Y)) \
        return ERR_MALLOC

    /* pop */
    uint8_t code = dir [y][x];
    if (code == 0)
    {
      collect_sink (x, y);
      continue;
    }
    
    uint8_t nbr = pop_bit (&code), inv = (nbr+4)%8;
    int xnbr = x + flow_neighbor[nbr].x, ynbr = y + flow_neighbor[nbr].y;

    if (IS_OUTSIDE (xnbr,ynbr))
    {
      collect_sink (xnbr, ynbr);
      continue;
    }

    /* update downstream accumulation*/
    acc [ynbr][xnbr] += acc [y][x];
    assert ( in_degree [ynbr][xnbr] & ((uint8_t) 1 << inv) );
    if ( ! (in_degree [ynbr][xnbr] &= ~((uint8_t) 1 << inv)) )
    {
      if (queue_push (&queue, xnbr, ynbr))
        return ERR_MALLOC;
    }
    else
      collect_sink (xnbr, ynbr);

    #undef collect_sink
  }
  assert ("Kahn's topological consistency" && processed == hw );

  queue_destroy (&queue);
  grid_free (in_degree);
  return 0;
}

static void log_accumulation (double ** acc, double ** l)
{
  GridData gd = grid_data (acc);
  int w = gd.width, h = gd.height;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      l [y][x] = log (acc [y][x]);
}

/*
.. identify connected flat patches and give them integer tags to identify.
.. Tag ID 0 is reserved for non-flat points.
*/
static int flat_tag (double ** elevation, int ** tag, int * ntags)
{
  #define push(x, y) do {                 \
    visited [y][x] = 1;                   \
    if (queue_push (&q, x,y))             \
      return ERR_MALLOC;                  \
  } while (0)

  assert ( !(tag == NULL || elevation == NULL) );
  Queue q;
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** visited = grid (uint8_t, w, h);
  if (!visited || queue_init (&q))
    return ERR_MALLOC;

  /* mark all points as tag 0*/
  memset (& tag  [-2][-2], 0, (w+4) * (h+4) * sizeof (int));

  /* assigning tag and type */
  int itag = 1;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      if (visited [y][x])
        continue;

      const int yb = y, xb = x;
      int count = 0;
      double e = elevation [y][x];

      push (x,y);

      while (queue_pop (&q, &x, &y))
      {
        tag [y][x] = itag;
        if (IS_OUTSIDE (x,y))
          continue;
        const int xb = x, yb = y;
        for (y = yb-1; y <= yb+1; ++y)
          for (x = xb-1; x <= xb+1; ++x)
          {
            if (visited [y][x])
              continue;
            if (!(e == elevation [y][x]))
              continue;
            count++;

            push (x,y);
          }
      }
      x = xb, y = yb;
      if (!count) /* not flat since no neighbors found at same elevation*/
        tag  [y][x] = 0;
      else
        itag++;
    }
  *ntags = itag;

  grid_free (visited);
  queue_destroy (&q);

  return 0;

  #undef push
}

static void accumulation_bug (int ** tag, uint8_t ** dir, uint8_t ** invDir)
{
  GridData gd = grid_data (tag);
  int w = gd.width, h = gd.height;
  double ** elevation = grid (double, w, h);
  uint8_t ** bug    = grid (uint8_t, w, h);
  D8 (elevation, w, h, dir, invDir);
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      bug [y][x] = !dir [y][x] ? 255 : tag [y][x] ? 16 : 0;
  pgm (bug, "no-drain-bug.pgm");
  grid_free (bug);
}

static int flat_type (double ** elevation, int ** tag, uint8_t ** type)
{
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** dir     = grid (uint8_t, w, h);
  uint8_t ** invDir  = grid (uint8_t, w, h);
  if (!dir || !invDir)
    return ERR_MALLOC;
  //DInfty (elevation, w, h, dir, invDir);
  D8 (elevation, w, h, dir, invDir);
//accumulation_bug (tag, dir, invDir);
  for (int y=-2; y<h+2; ++y)
    for (int x=-2; x<w+2; ++x)
    {
      type [y][x] =
        !tag [y][x]      ? 0           :
        dir [y][x]       ? FLAT_OUTLET :
        invDir [y][x]    ? FLAT_INLET  :
        IS_OUTSIDE (x,y) ? FLAT_OUTLET :
                           FLAT        ;
    }
  grid_free (dir);
  grid_free (invDir);
  return 0;
}

static int drain_validity (double ** elevation, int ** tag, int ntags)
{
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** dir     = grid (uint8_t, w, h);
  uint8_t ** invDir  = grid (uint8_t, w, h);
  assert (dir && invDir);
  //DInfty (elevation, w, h, dir, invDir);
  D8 (elevation, w, h, dir, invDir);

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
    uint8_t ** stagnant = grid (uint8_t, w, h); assert (stagnant);
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

/*
.. depression filling algorithm. also called as priority flood fill algorithm
*/
int flow_remove_pits (DEM * dem)
{
  int w = dem->w, h = dem->h;
  assert (w > 0 && h > 0);
  double ** elevation = dem->raster;
  assert (elevation != NULL);

  /* 'reservoir' non-zero tag to identify flat patches (possibly reservoirs) */
  int ** reservoir = grid (int, w, h);
  int ** visited = reservoir;
  MinPQ pq;
  if (reservoir == NULL || pq_create (&pq))
    return ERR_MALLOC;

  #define push(X,Y)                                                    \
    if (!visited [Y][X])                                               \
    {                                                                  \
      visited [Y][X] = 1;                                              \
      if (!isnan (elevation [Y][X]))                                   \
        if (pq_push (&pq, elevation [Y][X], X, Y))                     \
          return ERR_MALLOC;                                           \
    }

  /*
  .. tag each flat patches as reservoir and tag with a non-zero ID.
  .. Since we reuse reservoir [] as visited [], any points tagged with
  .. non-zero ID (corresponding to reservoir) will be skipped
  */
//fixme : not skipping reservoirs.
int ntags;
if (flat_tag (elevation, reservoir, &ntags))
  return ERR_MALLOC;

pgm (reservoir, "reservoirs.pgm");

  /* boundary points are pushed to the priority queue */
  for (int y=0; y<h; ++y)
  {
    push (0, y);
    push (w-1, y);
  }
  for (int x=0; x<w; ++x)
  {
    push (x, 0);
    push (x, h-1);
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
      if (IS_OUTSIDE (xnbr, ynbr))
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

  grid_free (reservoir);
  pq_free   (&pq);

  return 0;
}

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


#if 0
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

#endif


static int handle_flat (double ** elevation)
{

  #define push(q) do {                              \
      visited [y][x] = 1;                           \
      if (queue_push (&q, x, y)) return ERR_MALLOC; \
    } while (0)

  assert ( elevation != NULL );
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height, hw = w * h;

  uint8_t ** visited = grid (uint8_t, w, h);
  uint8_t ** type    = grid (uint8_t, w, h);
  int     ** tag     = grid (int,     w, h);
  if (type == NULL || visited == NULL || tag == NULL)
    return ERR_MALLOC;

  int ntags;
  /* first run the flood fill algorithm */
  if (flow_remove_pits ( & (DEM) {.raster = elevation, .w = w, .h = h} ) )
    return ERR_MALLOC;
  if (flat_tag (elevation, tag, &ntags))
    return ERR_MALLOC;
drain_validity (elevation, tag, ntags);
  if (flat_type (elevation, tag, type))
    return ERR_MALLOC;

pgm (tag, "flat-tag.pgm");
pgm (type, "flat-type.pgm");


  /* Garbrecht & Martz (1997) algorithm 
  .. fixme (1) need special case if the flat patch corresponds to reservoir
  .. (2) Cannot handle cases where flat patches shares boundary
  .. (3) I have an idea that you may use a Min heap for outlets such that
  .. a priority may be established for each drains out of the flat. (warning it is O(N log N))
  */
  Queue in, out;
  if (queue_init (&in) || queue_init (&out))
    return ERR_MALLOC;

  /*
  .. [creating] Gradient (downslope) away from the inlets.
  .. We use fifo queue to BFS search in each flat patches, starting with
  .. inlet points (which has some higher elevation neighors) as the seed.
  .. BFS using bits & pop_bit() is tricky. FIFO is the best option. 
  */
  int ** g1 = grid (int, w, h);
  if (!g1) return ERR_MALLOC;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      if (type [y][x] == FLAT_INLET)
      {
        assert (tag [y][x]);
        g1 [y][x] = 1;
        push (in);
      }
  int x, y;
  while (queue_pop (&in, &x, &y))
  {
    const int xb = x, yb = y, itag = tag [y][x], g = g1 [y][x] + 1;
    for (y = yb-1; y <= yb+1; ++y)
      for (x = xb-1; x <= xb+1; ++x)
        if (!visited [y][x] && tag [y][x] == itag)
        {
          push (in);
          g1 [y][x] = g;
        }
  }
  /* find max of g1 [] corresponding to each tag */
  int * max = malloc (ntags * sizeof (int));
  if (!max) return ERR_MALLOC;
  memset (max, 0, ntags * sizeof (int));
  for (int y=-1; y<h+1; ++y)
    for (int x=-1; x<w+1; ++x)
      if (tag [y][x])
      {
        int t = tag [y][x], g = g1 [y][x];
        if (g > max [t])
          max [t] = g;
      }
  /*
  .. inverting g1 so that the inlets are at higher elevation and
  .. flattens as you move inwards
  */
  for (int y=-1; y<h+1; ++y)
    for (int x=-1; x<w+1; ++x)
      if (tag [y][x])
      {
        g1 [y][x] = max [tag [y][x]] - g1 [y][x] + 1;
      }


  /* resetting visited [] for reuse */
  memset (& visited [-2][-2], 0, (w+4)*(h+4)*sizeof (uint8_t));

  /*
  .. [creating] Gradient toward the outlets.
  .. We use fifo queue to BFS search in each flat patches, starting with
  .. outlet points (which has some lower elevation neighors) as the seed
  .. includes one layer of ghost cells as outlet drains may fall outside [0:w-1][0:h-1]
  */
  
  int ** g2 = grid (int, w, h);
  if (!g2) return ERR_MALLOC;
  for (int y=-1; y<h+1; ++y)
    for (int x=-1; x<w+1; ++x)
      if (type [y][x] == FLAT_OUTLET)
      {
        assert (tag [y][x]);
        g2 [y][x] = 1;
        push (out);
      }
  while (queue_pop (&out, &x, &y))
  {
    const int xb = x, yb = y, itag = tag [y][x], g = g2 [y][x] + 1;
    for (y = yb-1; y <= yb+1; ++y)
      for (x = xb-1; x <= xb+1; ++x)
        if (!visited [y][x] && tag [y][x] == itag)
        {
          push (out);
          g2 [y][x] = g;
        }
  }

pgm (g1,"g1.pgm");
pgm (g2,"g2.pgm");


  /* adding perturbation weighted by g1 and g2
  double eps = 1.e-3, W1 = 1., W2 = 2.;
  for (int y=-1; y<h+1; ++y)
    for (int x=-1; x<w+1; ++x)
    {
      elevation [y][x] += eps * (W1 * (double) g1 [y][x] + W2 * (double) g2 [y][x]);
    }
  */
  uint8_t ** dir     = grid (uint8_t, w, h);
  uint8_t ** invDir  = grid (uint8_t, w, h);
  if (!dir || !invDir)
    return ERR_MALLOC;
  D8 (elevation, w, h, dir, invDir);
  int W1 = 0, W2 = 2;
  double a = 1., b = sqrt (2.);
  const double delta [8] = {a, b, a, b, a, b, a, b};
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      /* if (flood filled) elevation can find a direction, we rely on it */
      if (dir [y][x])
        continue;
      /* for non-tagged (i.e not flat) points at boundary, we skip */
      if (!tag [y][x] && (x == 0 || x == w-1 || y == 0 || y == h-1))
        continue;
      /* fallback. use g1 and g2 to resolve dir */
      double max = 0.;
      int itag   = tag [y][x];
      int val    = W1 * g1 [y][x] + W2 * g2 [y][x];
      int nbr    = -1;
      for (int c=0; c<8; ++c)
      {
        int xnbr = x + flow_neighbor [c].x;
        int ynbr = y + flow_neighbor [c].y;
        if (tag [ynbr][xnbr] != itag)
          continue; /* warning : skipping this will create cycles */
        int nbrVal = W1 * g1 [ynbr][xnbr] + W2 * g2 [ynbr][xnbr];
        double downslope = (double) (val - nbrVal) / delta [c];
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
pgm (dir, "d8.pgm");
pgm (invDir, "invDir.pgm");
  FlowNetwork network = {.w = w, .h = h, .dir = dir, .invDir = invDir};
  flow_accumulation (&network, (DEM) {.raster = elevation}, NULL);
  double ** acc = network.accumulation;
  double ** log_acc = grid (double, w, h);
  log_accumulation (acc, log_acc);
  pgm (log_acc, "accumulation.pgm");
grid_free (log_acc);
grid_free (network.accumulation);

    DAG_validity (dir, invDir);
    accumulation_bug (tag, dir, invDir);

  free (max);
  queue_destroy (&in);
  queue_destroy (&out);
  grid_free (type);
  grid_free (visited);
  grid_free (tag);
  grid_free (dir);
  grid_free (invDir);
  grid_free (g1);
  grid_free (g2);

  #undef push

  return 0;
}

static int pits (MinPQ * pq, double ** raster, uint8_t ** dir, uint8_t ** invDir)
{

  assert ( ! (pq == NULL || raster == NULL || dir == NULL || invDir == NULL) );

  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height, hw = w * h;
  uint8_t ** visited = grid (uint8_t, w, h);
  uint8_t ** bit_stack = grid (uint8_t, w, h);
  if (visited == NULL || bit_stack == NULL || pq_create (pq))
    return ERR_MALLOC;

  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      if (visited [y][x])
        continue;

      int ybackup = y;
      int xbackup = x;
      do
      {
        visited [y][x] = 1u;
        if (!dir [y][x])
        {

          if (!invDir [y][x])
          {
  /* grouping connected flat points as a single pit */
  int depth = 0;
  visited [y][x] = 0;
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
          int xnbr = x + flow_neighbor[c].x;
          int ynbr = y + flow_neighbor[c].y;
          if ( c != nbr && IS_FLAT (xnbr,ynbr) && !visited [ynbr][xnbr] && !IS_OUTSIDE (xnbr, ynbr))
          {
            bit_stack [y][x] |= (uint8_t) 1 << c;
            dir [ynbr][xnbr]  = (uint8_t) 1 << ((c+4)%8);
          }
        }
        invDir [y][x] |= bit_stack [y][x];
      }

      if (!bit_stack [y][x])
        break;

      /* remove the popped bit so that you don't traverse it again */
      int nbr = pop_bit ( &bit_stack [y][x] );
      x += flow_neighbor [nbr].x;
      y += flow_neighbor [nbr].y;
      assert (depth++ < hw);

    } while (1);

    if (!depth--)
      break;

    uint8_t code = dir [y][x]; assert (code);
    int nbr = pop_bit ( &code );
    x += flow_neighbor [nbr].x;
    y += flow_neighbor [nbr].y;
  } while (1);
  /* end of grouping connected flat points */
          }

          /* warning : the (x,y) may fall outside [0:w-1]x[0:h-1] */
          if (pq_push (pq, raster [y][x], x, y))
            return ERR_MALLOC;
          break;
        }
        uint8_t code = dir [y][x];
        int nbr = pop_bit (&code);
        x += flow_neighbor [nbr].x;
        y += flow_neighbor [nbr].y;
      } while (!visited [y][x]);
      y = ybackup;
      x = xbackup;

    }  /* end of x, y loops */

  grid_free (visited);
  grid_free (bit_stack);
  DAG_validity (dir, invDir);
  pq_truncate (pq); 

  return 0;
}

static int pit_color (MinPQ pq, uint8_t ** dir, uint8_t ** invDirConst)
{
  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height, hw = w * h;
  uint8_t ** color_tag  = grid (uint8_t, w, h);
  uint8_t ** color_elev = grid (uint8_t, w, h);
  if (color_tag == NULL || color_elev == NULL)
    return ERR_MALLOC;
  uint8_t ** invDir = grid (uint8_t, w, h);
  if (!invDir)
    return ERR_MALLOC;
  grid_copy (uint8_t, invDirConst, invDir, w, h);
  PQNode * nodes = pq.data;
  double min = DBL_MAX, max = DBL_MIN;
  for (int i=0; i<pq.size; ++i)
  {
    double elev = nodes [i].key;
    if (elev > max) max = elev;
    if (elev < min) min = elev;
  }
  double den = max - min;
  uint8_t tag = 0;
  for (int i=0; i<pq.size; ++i)
  {
    double elev = nodes [i].key;
    int x = nodes[i].x;
    int y = nodes[i].y;
    tag = ((int) tag + 1) % 255;
    tree_cell_start (invDir);
      color_elev [y][x] = (uint8_t) ((elev-min) * 255. /den); 
      color_tag  [y][x] = tag;
    tree_cell_end ();
  }
  /* pgm plots */
  pgm (color_tag, "pit-tag.pgm");
  pgm (color_elev, "pit-elev.pgm");
  /* free */
  grid_free (invDir);
  grid_free (color_tag);
  grid_free (color_elev);
  return 0;
}

int flow_routine (DEM dem)
{
  if (dem.raster == NULL || dem.w < 2 || dem.h < 2)
    return ERR_DATA_MISSING;
  int w = dem.w, h = dem.h;
  const double * const * raster = (const double * const *) dem.raster;
  double ** elevation           = grid (double, w, h); /* for a copy of raster */
  uint8_t ** dir                = grid (uint8_t, w, h);
  uint8_t ** invDir             = grid (uint8_t, w, h);
  int err;
  if (dir == NULL || invDir == NULL || elevation == NULL)
    return ERR_MALLOC;
  grid_copy (double, raster, elevation, w, h);

  #if 0
  /* create a forest of DAGs using elevation gradient*/
  D8 (elevation, w, h, dir, invDir);

  #endif

  /* flood fill test */ 
  handle_flat (elevation);
  //D8 (elevation, w, h, dir, invDir);

  /* create a min heap with key corresponding to pit centers (lowest point of each pit) */
  MinPQ pq;
  err = pits (&pq, elevation, dir, invDir);
  if (err) return err;
  pit_color (pq, dir, invDir);
  pq_free (&pq);


  //test
  pgm_hillshade (elevation, "dem-corrected.pgm", w, h, 315, 45, 30, 1);
  //pgm (dir, "d8.pgm");
  //pgm (invDir, "invDir.pgm");

  /* freeing all memory associated */
  grid_free (elevation);
  grid_free (dir);
  grid_free (invDir);

  return 0;
}

#undef error
#undef IS_FLAT
#undef IS_OUTSIDE
#undef FLAT_INLET
#undef FLAT_OUTLET
#undef FLAT
