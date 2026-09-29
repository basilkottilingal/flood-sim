#include "grid.h"
#include "pgm.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define not_unused(v) (void)v

static int limits (double ** raster, int w, int h, double * min, double * max)
{
  int foundnan = 0;
  *min = DBL_MAX, *max = DBL_MIN;
  for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
      if (isnan (raster [y][x]))
      {
        foundnan = 1;
        continue;
      }
      if (raster [y][x] < *min)
        *min = raster [y][x];
      if (raster [y][x] > *max)
        *max = raster [y][x];
    }
  return foundnan;
}

int pgm_grayscale (double ** raster, const char * out, int width, int height)
{
  if (raster == NULL || width < 1 || height < 1 || out == NULL || out [0] == '\0')
    return -1;

  FILE *fp = fopen (out, "wb");
  if (fp == NULL)
    return -1;

  /* Find range */
  double min, max;
  int foundNan = limits (raster, width, height, &min, &max);
  not_unused (foundNan);

  /* PGM header */
  fprintf(fp, "P5\n%d %d\n255\n", width, height);

  /* Convert to grayscale */
  if (max == min)
  {
    double pixel = 0.0f;
    for (int i=0; i<width*height; ++i)
      fwrite (&pixel, 1, 1, fp);
    fclose (fp);
    return 0;
  }

  double den = max - min;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
    {
      double v = isnan (raster [y][x]) ? 0. : (raster [y][x] - min) / den;
      unsigned char pixel = (unsigned char)(v * 255.0);
      fwrite (&pixel, 1, 1, fp);
    }

  fclose (fp);
  return 0;  
}

static unsigned char *
contour_pgm (double ** raster, int w, int h, double *levels, int nlevels)
{
  unsigned char * img = malloc ((size_t)w * h);
  if (img == NULL)
    return NULL;

  /* white background */
  memset (img, 255, (size_t)w * h);

  for (int k = 0; k < nlevels; ++k)
  {
    const double level = levels [k];
    const unsigned char color = 128u - 128u * (unsigned char) k / nlevels;
    for (int y = 0; y < h - 1; ++y)
    {
      for (int x = 0; x < w - 1; ++x)
      {

        /*
        ..   v0 -- e0 -- v1
        ..   |           |
        ..   e3          e1
        ..   |           |
        ..   v3 -- e2 -- v2
        */
        const double v0 = raster [y][x];
        const double v1 = raster [y][x+1];
        const double v2 = raster [y+1][x+1];
        const double v3 = raster [y+1][x];

        if (isnan (v0) || isnan (v1) || isnan (v2) || isnan (v3))
          continue;

        int c = 0;
        if (v0 >= level) c |= 1 << 0;
        if (v1 >= level) c |= 1 << 1;
        if (v2 >= level) c |= 1 << 2;
        if (v3 >= level) c |= 1 << 3;
        if (c == 0 || c == 15)
          continue;

        double ex[4], ey[4];
        int edge[4]; not_unused (edge);
        int n = 0;

        /* intesections */
        #define intersection(C,N,M) (((C >> M) ^ (C >> N)) & 1u)
        if (intersection (c, 0, 1))
        {
          double t = (level - v0) / (v1 - v0);
          ex[n] = x + t;
          ey[n] = y;
          edge[n++] = 0;
        }

        if (intersection (c, 1, 2))
        {
          double t = (level - v1) / (v2 - v1);
          ex[n] = x + 1;
          ey[n] = y + t;
          edge[n++] = 1;
        }

        if (intersection (c, 2, 3))
        {
          double t = (level - v2) / (v3 - v2);
          ex[n] = x + 1 - t;
          ey[n] = y + 1;
          edge[n++] = 2;
        }

        if (intersection (c, 3, 0))
        {
          double t = (level - v3) / (v0 - v3);
          ex[n] = x;
          ey[n] = y + 1 - t;
          edge[n++] = 3;
        }
        #undef intersection

        /*
        .. Normally there are 2 intersections.
        .. Saddle cases have 4.
        ..
        .. Connect pairs. For a visualization this is
        .. sufficient; for mathematically exact saddle
        .. handling, use the cell-center value to decide
        .. the connectivity.
        */
        for (int i = 0; i + 1 < n; i += 2)
        {
          int x0 = (int)lround (ex[i]);
          int y0 = (int)lround (ey[i]);
          int x1 = (int)lround (ex[i + 1]);
          int y1 = (int)lround (ey[i + 1]);

          /*
          .. draw the segment using Bresenham.
          */
          int dx = abs (x1 - x0);
          int sx = x0 < x1 ? 1 : -1;
          int dy = -abs (y1 - y0);
          int sy = y0 < y1 ? 1 : -1;
          int err = dx + dy;

          for (;;)
          {
            if ((unsigned) x0 < (unsigned) w && (unsigned) y0 < (unsigned) h)
              img [y0 * w + x0] = color;

            if (x0 == x1 && y0 == y1)
              break;

            int e2 = 2 * err;

            if (e2 >= dy)
            {
              err += dy;
              x0 += sx;
            }

            if (e2 <= dx)
            {
              err += dx;
              y0 += sy;
            }
          }
        }
      }   /* for (x = 0:w-1)     */
    }     /* for (y = 0:h-1)     */
  }       /* for (k = 0:nlevels) */

  return img;
}

int pgm_contour_grayscale (double ** raster, const char * out, int w, int h, int ncontours)
{
  if (raster == NULL || w < 2 || h < 2 || ncontours < 2 || out == NULL || out [0] == '\0')
    return -1;

  FILE *fp = fopen (out, "wb");
  if (fp == NULL)
    return -1;

  /* Find contour levels*/
  double min, max;
  int foundNan = limits (raster, w, h, &min, &max);
  not_unused (foundNan);
  double * levels = malloc (ncontours * sizeof (double));
  if (levels == NULL)
  {
    fclose (fp);
    return -1;
  }
  double d = (max - min) / ncontours; 
  for (int i=0; i<ncontours; ++i)
    levels [i] = min + d * (double) i;

  unsigned char * img = contour_pgm (raster, w, h, levels, ncontours);
  if (img == NULL)
  {
    free (levels);
    fclose (fp);
    return -1;
  }

  /* PGM header */
  fprintf (fp, "P5\n%d %d\n255\n", w, h);
  fwrite (img, 1, (size_t) w*h, fp);

  free (img);
  free (levels);
  fclose (fp);
  return 0;
}

int pgm_hillshade (
  double ** dem,
  const char * out,
  int width,
  int height,
  double azimuth,
  double altitude,
  double cellsize,
  double z_factor )
{
  if (dem == NULL || width < 3 || height < 3)
    return -1;

  double ** hillshade = grid (double, width, height);
  if (hillshade == NULL)
    return -1;

  if (azimuth < 0 || azimuth > 360)
    azimuth = 315;
  if (altitude < 0 || altitude > 90)
    altitude = 45;
  if (cellsize <= 0.)
    cellsize = 30.; /* assumes 30 m x 30 m raster resolution */
  double zenith = (90-altitude) * M_PI / 180.;
  double az     = (double) ((int) (360 - azimuth + 90) % 360) * M_PI / 180.;

  /* Horn 1981 algo */
  for (int y=0; y<height; ++y)
    for (int x=0; x<width; ++x)
    {
      double a = dem [y-1][x-1];
      double b = dem [y-1][x  ];
      double c = dem [y-1][x+1];
      double d = dem [y  ][x-1];
      double f = dem [y  ][x+1];
      double g = dem [y+1][x-1];
      double h = dem [y+1][x  ];
      double i = dem [y+1][x+1];
        
      double dzdx = ((c + 2*f + i) - (a + 2*d + g)) / (8*cellsize);
      double dzdy = ((g + 2*h + i) - (a + 2*b + c)) / (8*cellsize);

      double slope  = atan (z_factor * sqrt (dzdx*dzdx + dzdy*dzdy));
      double aspect = atan2 (dzdy, -dzdx);

      hillshade [y][x] = 255 * 
        fmax (0., cos (zenith) * cos (slope) + sin (zenith) * sin (slope) * cos (az - aspect));
    }

  int err = pgm_grayscale (hillshade, out, width, height);
  grid_free (hillshade);
  return err;
}
