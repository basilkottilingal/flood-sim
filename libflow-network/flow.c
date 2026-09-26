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

static inline int pop_bit (uint8_t * v)
{
  if (*v == 0)
    return -1;
#if defined(__GNUC__) || defined(__clang__)
  int idx = __builtin_ctz((unsigned int)*v);
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

static void ** grid (size_t s, int n, int m, int ng)
{
  /* warning : use this only for uint8_t and float */
  assert (n > 0 && m > 0 && ng >= 0);
  size_t size = (n + 2*ng) * (m + 2*ng) * s;

  if (s == sizeof (uint8_t))
  {
    uint8_t * mem = malloc (size);
    if (mem == NULL)
      error ("flow.h : grid () : malloc failed");
    memset (mem, 0, size);
    uint8_t ** memptr = malloc ((m + 2*ng) * sizeof (uint8_t *));
    if (memptr == NULL)
      error ("flow.h : grid () : malloc failed");
    for (int j = -ng; j < m + ng; ++j)
      memptr [j + ng] = & mem [(j + ng) * ( n + 2*ng ) + ng];
    return (void **) (memptr + ng);
  }
  if (s == sizeof (float))
  {
    float * mem = malloc (size);
    if (mem == NULL)
      error ("flow.h : grid () : malloc failed");
    memset (mem, 0, size);
    float ** memptr = malloc ((m + 2*ng) * sizeof (float *));
    if (memptr == NULL)
      error ("flow.h : grid () : malloc failed");
    for (int j = -ng; j < m + ng; ++j)
      memptr [j + ng] = & mem [(j + ng) * ( n + 2*ng ) + ng];
    return (void **) (memptr + ng);
  }
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

  int n = dem.n, m = dem.m;
  assert (n > 0 && m > 0);
  uint8_t ** dir = (uint8_t **) grid (sizeof (uint8_t), n, m, 1);
  float ** raster = dem.raster;

  *network = (FlowNetwork)
    {
      .dir    = dir,
      .accumulation = NULL,
      .n      = n,
      .m      = m,
      .type   = type,
    };

  float a = 1., b = 1. / sqrt (2.);
  const float delta [8] = {a, b, a, b, a, b, a, b};

  switch (type)
  {
    case FLOW_D8  :
      for (int j=0; j<m; ++j)
        for (int i=0; i<n; ++i)
        {
          float max = 0.f, val = raster [i][j];
          uint8_t code = 0;
          if ( isnan (val) )
            continue;
          for (int c=0; c<8; ++c)
          {
            float nbrVal = raster [i + neighbor[c].x][j + neighbor[c].y];
            if (isnan (nbrVal))
              continue;
            float downslope = (val - nbrVal) / delta [c];
            if (downslope > max)
              code = (uint8_t) 1 << c, max = downslope;
          }
          dir [i][j] = code;
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

  width  = width  > network.n ? network.n > 1024 ? 1024 : network.n : width;
  height = height > network.m ? network.m > 1024 ? 1024 : network.m : height;

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

  width  = width  > network.n ? network.n > 1024 ? 1024 : network.n : width;
  height = height > network.m ? network.m > 1024 ? 1024 : network.m : height;

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

  int n = network->n, m = network->m, nm = n * m;
  float ** acc = network->accumulation = (float **) grid (sizeof (float), n, m, 1);
  float ** raster = dem.raster;
  uint8_t ** in_degree = (uint8_t **) grid (sizeof (uint8_t), n, m, 1);
  uint8_t ** dir = network->dir;

  /* find the number of neighbors who are source to (i,j) */
  for (int j=0; j<m; ++j)
    for (int i=0; i<n; ++i)
    {
      acc [i][j] = 1.; /* you first add the "rainfall" at (i,j) to the flow[i][j]*/
      uint8_t code = dir [i][j];
      if ( code == 0 )
        continue;
      uint8_t nbr = pop_bit (&code); /* (i, j) is a source to nbr */
      int inbr = i + neighbor[nbr].x, jnbr = j + neighbor[nbr].y;
      in_degree [inbr][jnbr] |= (uint8_t) 1 << nbr;
    }

  /* first in first out queue */
  struct index { int i, j; } * queue = malloc (nm * sizeof (struct index));
  if (queue == NULL)
    error ("flow_accumulation () : malloc () failed");
  int push_at = 0, pop_at = 0;
  int processed = 0;
  for (int j=0; j<m; ++j)
    for (int i=0; i<n; ++i)
    {
      if (isnan (raster [i][j]))
      {
        processed ++;
        continue;
      }
      if (in_degree [i][j] == 0)
        queue [push_at++] = (struct index) {i,j};
    } 

  /* Kahn's topological propogation */
  while (push_at != pop_at)
  {
    processed ++;

    /* pop */
    struct index Idx = queue [pop_at]; pop_at = (pop_at+1) % nm;
    int i = Idx.i, j = Idx.j;
    uint8_t code = dir [i][j];
    if (code == 0)
      continue;
    
    uint8_t nbr = pop_bit (&code); /* (i, j) is a source to nbr */
    int inbr = i + neighbor[nbr].x, jnbr = j + neighbor[nbr].y;

    if (inbr<0 || inbr>=n || jnbr<0 || jnbr>=m)
      continue;

    /* downstream */
    acc [inbr][jnbr] += acc [i][j];

    assert ( in_degree [inbr][jnbr] & ((uint8_t) 1 << nbr) );
    /* push */
    if ( (in_degree [inbr][jnbr] &= ~((uint8_t) 1 << nbr)) == 0 )
      queue [push_at] = (struct index) {inbr, jnbr}, push_at = (push_at + 1) % nm;

  }

  if (processed != nm)
    fprintf (stderr, "Kahn's topological propagation inconsistency");

  free (queue);
  grid_free (in_degree, 1);
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
