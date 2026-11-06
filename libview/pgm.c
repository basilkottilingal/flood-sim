#include "pgm.h"
#include "image.h"

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/* t in [0,1] (NaN => black) to an 8-bit gray image */
static int pgm_unit (const double * t, const char * out, int w, int h)
{
  if (t == NULL)
    return -1;
  uint8_t * img = malloc ((size_t) w * h);
  if (img == NULL)
    return -1;
  for (size_t p = 0; p < (size_t) w * h; ++p)
    img [p] = isnan (t [p]) ? 0 : (uint8_t) lround (255. * t [p]);
  int err = image_write (out, 1, w, h, img);
  free (img);
  return err;
}

static int pgm_typed (const void * g, size_t size, const char * out, int w, int h)
{
  if (!image_check (g, out, w, h, 1))
    return -1;
  double * t = image_normalize (g, size, w, h);
  int err = pgm_unit (t, out, w, h);
  free (t);
  return err;
}

int pgm_grayscale (double ** raster, const char * out, int w, int h)
{
  return pgm_typed (raster, sizeof (double), out, w, h);
}

int pgm_general (void * g, const char * out, GridData gd)
{
  int w = gd.width  > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.width;
  int h = gd.height > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.height;
  switch (gd.size)
  {
    case sizeof (uint8_t) :
    case sizeof (int) :
    case sizeof (double) :
      return pgm_typed (g, gd.size, out, w, h);
    default :
      return -1;
  }
}

int pgm_contour_grayscale (double ** raster, const char * out, int w, int h, int ncontours)
{
  if (!image_check (raster, out, w, h, 2) || ncontours < 1)
    return -1;

  int * label = image_contour (raster, w, h, ncontours);
  uint8_t * img = malloc ((size_t) w * h);
  int err = -1;
  if (label != NULL && img != NULL)
  {
    /* white background, lines from mid gray (lowest) to black (highest) */
    for (size_t p = 0; p < (size_t) w * h; ++p)
      img [p] = label [p] < 0 ? 255 :
        (uint8_t) (128 - 128 * label [p] / ncontours);
    err = image_write (out, 1, w, h, img);
  }
  free (img);
  free (label);
  return err;
}

int pgm_hillshade (double ** dem, const char * out, int w, int h,
  double azimuth, double altitude, double cellsize, double z_factor)
{
  if (!image_check (dem, out, w, h, 1))
    return -1;
  double * shade = image_hillshade (dem, w, h, azimuth, altitude, cellsize, z_factor);
  int err = pgm_unit (shade, out, w, h);
  free (shade);
  return err;
}
