#include "flow.h"
#include "flow-utils.h"
#include "geotiff.h"
#include "filemap.h"
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

void flow_network_free (FlowNetwork * network)
{
  assert (network->dir && network->invDir && network->accumulation && network->elevation);
  grid_free (network->dir);
  grid_free (network->invDir);
  grid_free (network->accumulation);
  grid_free (network->elevation);
  *network = (FlowNetwork) {0};
}

int flow_accumulation (double ** acc, uint8_t ** dir, uint8_t ** invDir)
{

  assert (acc && dir && invDir);
  GridData gd = grid_data (dir);
  int w = gd.width, h = gd.height, hw = h * w;
  uint8_t ** in_degree = grid (uint8_t, w, h);
  assert (in_degree); /* malloc failed ? */
  grid_copy (uint8_t, invDir, in_degree, w, h);

  /* every cell is a source to itself */
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      acc [y][x] = 1.;

  /* first in first out queue */
  Queue queue;
  assert (!queue_init (&queue));  /* malloc failed ? */
  int processed = 0;
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      #if 0
      if (isnan (raster [y][x]))
      {
        processed ++;
        continue;
      }
      #endif
      if (in_degree [y][x] == 0)
        if (queue_push (&queue, x, y))
          return ERR_MALLOC;
    }

  /* Kahn's topological propogation */
  int x, y;
  while (queue_pop (&queue, &x, &y))
  {
    processed ++;
    uint8_t code = dir [y][x];
    if (code == 0)
      continue;
    uint8_t nbr = pop_bit (&code);
    uint8_t inv = (nbr + 4) % 8;
    int xnbr    = x + flow_neighbor[nbr].x;
    int ynbr    = y + flow_neighbor[nbr].y;
    if (IS_OUTSIDE (xnbr,ynbr))
      continue;
    /* update downstream accumulation*/
    acc [ynbr][xnbr] += acc [y][x];

    assert ( in_degree [ynbr][xnbr] & ((uint8_t) 1 << inv) );
    if ( ! (in_degree [ynbr][xnbr] &= ~((uint8_t) 1 << inv)) )
      if (queue_push (&queue, xnbr, ynbr))
        return ERR_MALLOC;
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

  Queue q;
  assert (tag && elevation);
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** visited = grid (uint8_t, w, h);
  assert (visited && !queue_init (&q)); /* malloc failed ? */

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

static void accumulation_bug (int ** tag, uint8_t ** dir)
{
  GridData gd = grid_data (tag);
  int w = gd.width, h = gd.height;
  uint8_t ** bug    = grid (uint8_t, w, h);
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
      bug [y][x] = !dir [y][x] ? 255 : tag [y][x] ? 16 : 0;
  pgm (bug, "no-drain-bug.pgm");
  grid_free (bug);
}

/*
.. identify cells as flat, flat inlet, flat outlet and none (0)
*/
static int flat_type (double ** elevation, int ** tag, uint8_t ** type)
{
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  uint8_t ** dir     = grid (uint8_t, w, h);
  uint8_t ** invDir  = grid (uint8_t, w, h);
  if (!dir || !invDir)
    return ERR_MALLOC;
  D8 (elevation, dir, invDir);
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

/*
.. depression filling algorithm. also called as priority flood fill algorithm
*/
int flow_remove_pits (double ** elevation)
{
  assert (elevation);
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;
  assert (w > 0 && h > 0);

  /* 'reservoir' non-zero tag to identify flat patches (possibly reservoirs) */
  int ** reservoir = grid (int, w, h);
  int ** visited   = reservoir;
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

  /* warning : need additional info to distinguish real reservoirs from
  .. other flattened patches */
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

/*
.. Note : input should be the output of priority fill corrected elevation
*/
static int Garbrecht_Martz (double ** elevation, uint8_t ** dir, uint8_t ** invDir)
{

  #define push(q) do {                              \
      visited [y][x] = 1;                           \
      if (queue_push (&q, x, y)) return ERR_MALLOC; \
    } while (0)

  assert ( elevation != NULL );
  GridData gd = grid_data (elevation);
  int w = gd.width, h = gd.height;

  uint8_t ** visited = grid (uint8_t, w, h);
  uint8_t ** type    = grid (uint8_t, w, h);
  int     ** tag     = grid (int,     w, h);
  assert (type && visited && tag); /* malloc failed ? */

  int ntags;
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


  #if 0
  double eps = 1.e-3, W1 = 1., W2 = 2.;
  for (int y=-1; y<h+1; ++y)
    for (int x=-1; x<w+1; ++x)
TanDEM-X    {
      elevation [y][x] += eps * (W1 * (double) g1 [y][x] + W2 * (double) g2 [y][x]);
    }
  #endif
  /*
  .. G-M algorithm add perturbation weighted by g1 and g2 to elevation
  .. and recomupute "dir". We instead use a fall back mechanism,
  .. when dir [] == 0 (where dir is computed by D8), we use g1[] + g2[]
  .. to compute best local downslope
  */
  D8 (elevation, dir, invDir);
pgm (dir, "dir-0.pgm");
pgm (invDir, "invDir-0.pgm");
  int W1 = 1, W2 = 2;
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

DAG_validity (dir, invDir);
accumulation_bug (tag, dir);

  free (max);
  queue_destroy (&in);
  queue_destroy (&out);
  grid_free (type);
  grid_free (visited);
  grid_free (tag);
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

/* create a flow network */
int flow_network (DEM dem, FLOW_DRAIN type, FlowNetwork * network)
{
  assert (dem.raster && network);
  assert ("not implemented" && type == FLOW_D8);

  int w = dem.w, h = dem.h;
  const double * const * raster = (const double * const *) dem.raster;
  double  ** elevation          = grid (double,  w, h);
  uint8_t ** dir                = grid (uint8_t, w, h);
  uint8_t ** invDir             = grid (uint8_t, w, h);
  double  ** accumulation       = grid (double,  w, h);
  assert (dir && invDir && elevation && accumulation); /* malloc failed ? */
  grid_copy (double, raster, elevation, w, h);

  * network = (FlowNetwork)
    {
      .raster = raster,
      .dir    = dir,
      .invDir = invDir,
      .accumulation = accumulation,
      .elevation = elevation, 
      .type   = type,
      .w = w,
      .h = h
    };

  /* 1. priority flood filling algorthm */
  pgm_hillshade (elevation, "dem.pgm", w, h, 315, 45, 30, 1);
  if (flow_remove_pits (elevation))
    return ERR_MALLOC;
  pgm_hillshade (elevation, "dem-corrected.pgm", w, h, 315, 45, 30, 1);

  /* 2. (modified) Grabrecht Martz algorithm */
  if (Garbrecht_Martz (elevation, dir, invDir))
    return -1;

pgm (dir, "dir.pgm");
pgm (invDir, "invDir.pgm");

  /* 3. accumulation based on the new dir */
  if (flow_accumulation (accumulation, dir, invDir))
    return -1;

  /* following are only for debugging */
  double ** acc_log = grid (double,  w, h);
  assert (acc_log);
  log_accumulation (accumulation, acc_log);
  pgm (acc_log, "accumulation.pgm");

  /* create a min heap with key corresponding to pit centers (lowest point of each pit) 
  MinPQ pq;
  int err = pits (&pq, elevation, dir, invDir);
  if (err) return err;
  pit_color (pq, dir, invDir);
  pq_free (&pq);
*/
  
  return 0;
}

#undef error
#undef IS_FLAT
#undef IS_OUTSIDE
#undef FLAT_INLET
#undef FLAT_OUTLET
#undef FLAT
