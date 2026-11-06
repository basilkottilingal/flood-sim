#include "image.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* M_PI is POSIX, not ISO C99/C11 */
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

int image_check (const void * raster, const char * out, int w, int h, int minsize)
{
  return raster != NULL && out != NULL && out [0] != '\0' &&
    w >= minsize && h >= minsize;
}

double image_value (const void * g, size_t size, int x, int y)
{
  switch (size)
  {
    case sizeof (uint8_t) : return ((uint8_t * const *) g) [y][x];
    case sizeof (int)     : return ((int     * const *) g) [y][x];
    case sizeof (double)  : return ((double  * const *) g) [y][x];
    default               : return NAN;
  }
}

long image_limits (const void * g, size_t size, int w, int h, double * min, double * max)
{
  long n = 0;
  *min = DBL_MAX, *max = -DBL_MAX;
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
    {
      double v = image_value (g, size, x, y);
      if (!isfinite (v))
        continue;
      if (v < *min) *min = v;
      if (v > *max) *max = v;
      n++;
    }
  return n;
}

double * image_normalize (const void * g, size_t size, int w, int h)
{
  double * t = malloc ((size_t) w * h * sizeof (double));
  if (t == NULL)
    return NULL;

  double min, max;
  image_limits (g, size, w, h, &min, &max);
  double scale = max > min ? 1. / (max - min) : 0.;

  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
    {
      double v = image_value (g, size, x, y);
      t [(size_t) y * w + x] = isfinite (v) ? (v - min) * scale : NAN;
    }
  return t;
}

double * image_hillshade (double ** dem, int w, int h, double azimuth,
  double altitude, double cellsize, double z_factor)
{
  if (azimuth < 0. || azimuth > 360.)
    azimuth = 315.;
  if (altitude < 0. || altitude > 90.)
    altitude = 45.;
  if (cellsize <= 0.)
    cellsize = 30.;   /* assumes 30 m x 30 m raster resolution */

  double * shade = malloc ((size_t) w * h * sizeof (double));
  if (shade == NULL)
    return NULL;

  const double zenith = (90. - altitude) * M_PI / 180.;
  const double az     = fmod (360. - azimuth + 90., 360.) * M_PI / 180.;
  const double cz = cos (zenith), sz = sin (zenith);

  for (int y = 0; y < h; ++y)
  {
    const int ym = y > 0 ? y - 1 : y, yp = y < h - 1 ? y + 1 : y;
    for (int x = 0; x < w; ++x)
    {
      const int xm = x > 0 ? x - 1 : x, xp = x < w - 1 ? x + 1 : x;
      const double e = dem [y][x];
      if (isnan (e))
      {
        shade [(size_t) y * w + x] = NAN;
        continue;
      }

      /*
      ..   a b c
      ..   d e f
      ..   g h i
      .. missing (NaN) neighbours take the central value e.
      */
      #define z(Y,X) (isnan (dem [Y][X]) ? e : dem [Y][X])
      double a = z (ym, xm), b = z (ym, x), c = z (ym, xp);
      double d = z (y,  xm),                f = z (y,  xp);
      double g = z (yp, xm), hh = z (yp, x), i = z (yp, xp);
      #undef z

      double dzdx = ((c + 2*f + i) - (a + 2*d + g)) / (8*cellsize);
      double dzdy = ((g + 2*hh + i) - (a + 2*b + c)) / (8*cellsize);

      double slope  = atan (z_factor * sqrt (dzdx*dzdx + dzdy*dzdy));
      double aspect = atan2 (dzdy, -dzdx);

      shade [(size_t) y * w + x] =
        fmax (0., cz * cos (slope) + sz * sin (slope) * cos (az - aspect));
    }
  }
  return shade;
}

/* Bresenham line from (x0,y0) to (x1,y1) tagging pixels with k */
static void contour_segment (int * img, int w, int h, int x0, int y0, int x1, int y1, int k)
{
  int dx = abs (x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -abs (y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy;

  for (;;)
  {
    if ((unsigned) x0 < (unsigned) w && (unsigned) y0 < (unsigned) h)
      img [(size_t) y0 * w + x0] = k;
    if (x0 == x1 && y0 == y1)
      break;
    int e2 = 2 * err;
    if (e2 >= dy) { err += dy; x0 += sx; }
    if (e2 <= dx) { err += dx; y0 += sy; }
  }
}

int * image_contour (double ** raster, int w, int h, int nlevels)
{
  int * img = malloc ((size_t) w * h * sizeof (int));
  if (img == NULL)
    return NULL;
  for (size_t p = 0; p < (size_t) w * h; ++p)
    img [p] = -1;

  double min, max;
  if (image_limits (raster, sizeof (double), w, h, &min, &max) == 0 || max <= min)
    return img;

  for (int k = 0; k < nlevels; ++k)
  {
    const double level = min + (max - min) * (k + 1) / (nlevels + 1);
    for (int y = 0; y < h - 1; ++y)
      for (int x = 0; x < w - 1; ++x)
      {
        /*
        ..   v0 -- e0 -- v1
        ..   |           |
        ..   e3          e1
        ..   |           |
        ..   v3 -- e2 -- v2
        */
        const double v [4] = { raster [y][x],   raster [y][x+1],
                               raster [y+1][x+1], raster [y+1][x] };
        /* corner offsets of v0..v3 */
        static const int cx [4] = {0, 1, 1, 0}, cy [4] = {0, 0, 1, 1};

        if (isnan (v[0]) || isnan (v[1]) || isnan (v[2]) || isnan (v[3]))
          continue;

        int c = 0;
        for (int j = 0; j < 4; ++j)
          if (v [j] >= level)
            c |= 1 << j;
        if (c == 0 || c == 15)
          continue;

        /* intersection on edge j (between corners j and j+1) */
        double ex [4], ey [4];
        int n = 0;
        for (int j = 0; j < 4; ++j)
        {
          int jn = (j + 1) & 3;
          if (!(((c >> j) ^ (c >> jn)) & 1))
            continue;
          double t = (level - v [j]) / (v [jn] - v [j]);
          ex [n] = x + cx [j] + t * (cx [jn] - cx [j]);
          ey [n] = y + cy [j] + t * (cy [jn] - cy [j]);
          n++;
        }

        /*
        .. Normally 2 intersections, saddles have 4. Pairing consecutive
        .. ones is sufficient for visualization (an exact treatment would
        .. use the cell-centre value to resolve the saddle).
        */
        for (int i = 0; i + 1 < n; i += 2)
          contour_segment (img, w, h,
            (int) lround (ex [i]),     (int) lround (ey [i]),
            (int) lround (ex [i + 1]), (int) lround (ey [i + 1]), k);
      }
  }
  return img;
}

int image_write (const char * out, int channels, int w, int h, const uint8_t * pixels)
{
  if (out == NULL || pixels == NULL || (channels != 1 && channels != 3))
    return -1;

  FILE * fp = fopen (out, "wb");
  if (fp == NULL)
    return -1;

  size_t n = (size_t) w * h * channels;
  int err = fprintf (fp, "P%d\n%d %d\n255\n", channels == 1 ? 5 : 6, w, h) < 0;
  err |= fwrite (pixels, 1, n, fp) != n;
  err |= fclose (fp) != 0;
  return err ? -1 : 0;
}
