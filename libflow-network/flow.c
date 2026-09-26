#include "geotiff.h"
#include "filemap.h"
#include "flow.h"

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
} ERR;

# if defined(__GNUC__) || defined(__clang__)
#   if defined(__has_builtin)
#     if __has_builtin(__builtin_ctz)
#       define HAVE_CTZ
#     endif
#   endif
# endif

static inline int pop_bit (uint8_t * v)
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

static void ** grid (size_t s, int w, int h, int ng)
{
  /* warning : use this only for uint8_t and float */
  assert (w > 0 && h > 0 && ng >= 0);
  size_t size = (w + 2*ng) * (h + 2*ng) * s;

  #define _grid_(type)                                          \
    if (s == sizeof (type))                                     \
    {                                                           \
      type * mem = malloc (size);                               \
      if (mem == NULL)                                          \
        error ("flow.h : grid () : malloc failed");             \
      memset (mem, 0, size);                                    \
      type ** memptr = malloc ((h + 2*ng) * sizeof (type *));   \
      if (memptr == NULL)                                       \
        error ("flow.h : grid () : malloc failed");             \
      for (int y = -ng; y < h + ng; ++y)                        \
        memptr [y + ng] = & mem [(y + ng) * (w + 2*ng) + ng];   \
      return (void **) (memptr + ng);                           \
    }

  _grid_ (uint8_t)
  _grid_ (float)

  #undef _grid_
  
  error ("flow.h : grid () : unknown type");
  return NULL;
}

#define grid_free(g,nghost) free (& g[-nghost][-nghost]), free (& g[-nghost])

/*
.. Neighbors in this order
..
..  3 2 1
..  4 . 0
..  5 6 7               
*/
const struct { int x, y; } neighbor [] = 
  { {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1} };

int flow_network (DEM dem, FLOW_DRAIN type, FlowNetwork * network)
{
  if (dem.raster == NULL)
    return -1;

  int w = dem.w, h = dem.h;
  assert (w > 0 && h > 0);
  uint8_t ** dir = (uint8_t **) grid (sizeof (uint8_t), w, h, 1);
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
            float nbrVal = raster [y + neighbor[c].y][x + neighbor[c].x];
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
  grid_free ( network.dir, 1 );
  if (network.accumulation == NULL)
    return;
  grid_free ( network.accumulation, 1 );
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
  if (network.accumulation == NULL || out == NULL || out [0] == '\0' || width < 0 || height < 0)
    return -1;

  width  = width  > network.w ? network.w > 1024 ? 1024 : network.w : width;
  height = height > network.h ? network.h > 1024 ? 1024 : network.h : height;

  FILE *fp = fopen (out, "wb");
  if (fp == NULL)
    return -1;

  /* Find range */
  float min = FLT_MAX, max = FLT_MIN;
  float ** raster = network.accumulation;
  assert (raster != NULL);

  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
    {
      if (raster [y][x] < min)
        min = raster [y][x];
      if (raster [y][x] > max)
        max = raster [y][x];
    }

  /* PGM header */
  fprintf(fp, "P5\n%d %d\n255\n", width, height);

  /* Convert to grayscale */
  if (max == min)
  {
    float pixel = 0.0f;
    for (int i=0; i<width*height; ++i)
      fwrite (&pixel, 1, 1, fp);
    fclose (fp);
    return 0;
  }

  float den = max - min;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
    {
      float v = (raster [y][x] - min) / den;
      unsigned char pixel = (unsigned char)(v * 255.0);
      fwrite (&pixel, 1, 1, fp);
    }

  fclose (fp);
  return 0;  
}

int flow_accumulation (FlowNetwork * network, DEM dem)
{
  if (network->dir == NULL || network->accumulation != NULL)
    return -1;

  int w = network->w, h = network->h, hw = w * h;
  float ** acc = network->accumulation = (float **) grid (sizeof (float), w, h, 1);
  float ** raster = dem.raster;
  uint8_t ** in_degree = (uint8_t **) grid (sizeof (uint8_t), w, h, 1);
  uint8_t ** dir = network->dir;

  /* find the number of neighbors who are source to (i,j) */
  for (int y=0; y<h; ++y)
    for (int x=0; x<w; ++x)
    {
      acc [y][x] = 1.; /* you first add the "rainfall" at (i,j) to the flow[i][j]*/
      uint8_t code = dir [y][x];
      if ( code == 0 )
        continue;
      uint8_t nbr = pop_bit (&code); /* (i, j) is a source to nbr */
      int xnbr = x + neighbor[nbr].x, ynbr = y + neighbor[nbr].y;
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
    int xnbr = x + neighbor[nbr].x, ynbr = y + neighbor[nbr].y;

    if (xnbr<0 || xnbr>=w || ynbr<0 || ynbr>=h)
      continue;

    /* update downstream accumulation*/
    acc [ynbr][xnbr] += acc [y][x];

    assert ( in_degree [ynbr][xnbr] & ((uint8_t) 1 << nbr) );
    /* push */
    if ( (in_degree [ynbr][xnbr] &= ~((uint8_t) 1 << nbr)) == 0 )
      queue [push_at] = (struct index) {xnbr, ynbr}, push_at = (push_at + 1) % hw;

  }

  if (processed != hw)
    fprintf (stderr, "Kahn's topological propagation inconsistency");

  free (queue);
  grid_free (in_degree, 1);
  return 0;
}

#if 0
int flow_accumulation_pit_removed (FlowNetwork * network, DEM dem)
{
  int n = network->n, m = network->m, hw = n * m;
  float ** h = grid (sizeof (float), n, m, 2);
  for (int j=-2; j<m+2; ++j)
    memcpy (& h[
  uint8_t ** in_degree = (uint8_t **) grid (sizeof (uint8_t), n, m, 1);
  uint8_t ** dir = network->dir;

  #define VISITED 0

  DEM pit_removed = {.raster = h};
  return flow_accumulation (network, pit_removed)
  grid_free (h, 2);

  #undef VISITED
}
#endif


void flow_network_error ( int type )
{
  switch (type)
  {
    default :
      fprintf (stderr, "unknonw flow network error\n");
  }
}

#undef error
