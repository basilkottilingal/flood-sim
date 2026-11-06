#include "ppm.h"
#include "image.h"

#include <stdlib.h>
#include <math.h>

typedef struct { double t; uint8_t r, g, b; } Stop;

/* piecewise linear colormaps. Each table ends with t = 1 */
static const Stop viridis [] = {
  {0.000,  68,   1,  84}, {0.125,  71,  44, 122}, {0.250,  59,  81, 139},
  {0.375,  44, 113, 142}, {0.500,  33, 144, 141}, {0.625,  39, 173, 129},
  {0.750,  92, 200,  99}, {0.875, 170, 220,  50}, {1.000, 253, 231,  37}
};
/* approximation of Google's turbo */
static const Stop turbo [] = {
  {0.000,  48,  18,  59}, {0.100,  67,  70, 171}, {0.200,  70, 117, 237},
  {0.300,  57, 162, 252}, {0.400,  27, 207, 212}, {0.500,  36, 236, 166},
  {0.600, 114, 254,  94}, {0.700, 175, 250,  55}, {0.800, 249, 186,  56},
  {0.900, 230, 100,  25}, {1.000, 122,   4,   3}
};
static const Stop terrain [] = {
  {0.000,  51,  51, 153}, {0.150,   0, 153, 255}, {0.250,   0, 204, 102},
  {0.500, 255, 255, 153}, {0.750, 128,  92,  84}, {1.000, 255, 255, 255}
};
static const Stop blues [] = {
  {0.000, 247, 251, 255}, {0.250, 198, 219, 239}, {0.500, 107, 174, 214},
  {0.750,  33, 113, 181}, {1.000,   8,  48, 107}
};
static const Stop gray [] = {
  {0.000,   0,   0,   0}, {1.000, 255, 255, 255}
};

static const uint8_t nan_rgb [3] = {255, 255, 255};

void ppm_color (PPMColormap cmap, double t, uint8_t rgb [3])
{
  if (isnan (t))
  {
    rgb [0] = nan_rgb [0], rgb [1] = nan_rgb [1], rgb [2] = nan_rgb [2];
    return;
  }
  const Stop * s;
  switch (cmap)
  {
    case PPM_TURBO   : s = turbo;   break;
    case PPM_TERRAIN : s = terrain; break;
    case PPM_BLUES   : s = blues;   break;
    case PPM_GRAY    : s = gray;    break;
    default          : s = viridis; break;
  }
  t = t < 0. ? 0. : t > 1. ? 1. : t;
  while (s [1].t < t)
    s++;
  double u = (t - s [0].t) / (s [1].t - s [0].t);
  rgb [0] = (uint8_t) lround (s [0].r + u * (s [1].r - s [0].r));
  rgb [1] = (uint8_t) lround (s [0].g + u * (s [1].g - s [0].g));
  rgb [2] = (uint8_t) lround (s [0].b + u * (s [1].b - s [0].b));
}

/* t in [0,1] (NaN => white) to an RGB image */
static uint8_t * ppm_unit (const double * t, int w, int h, PPMColormap cmap)
{
  if (t == NULL)
    return NULL;
  uint8_t * img = malloc ((size_t) w * h * 3);
  if (img == NULL)
    return NULL;
  for (size_t p = 0; p < (size_t) w * h; ++p)
    ppm_color (cmap, t [p], img + 3 * p);
  return img;
}

static int ppm_typed (const void * g, size_t size, const char * out,
  int w, int h, PPMColormap cmap)
{
  if (!image_check (g, out, w, h, 1))
    return -1;
  double * t = image_normalize (g, size, w, h);
  uint8_t * img = ppm_unit (t, w, h, cmap);
  int err = img == NULL ? -1 : image_write (out, 3, w, h, img);
  free (img);
  free (t);
  return err;
}

int ppm_colormap (double ** raster, const char * out, int w, int h, PPMColormap cmap)
{
  return ppm_typed (raster, sizeof (double), out, w, h, cmap);
}

int ppm_general (void * g, const char * out, GridData gd, PPMColormap cmap)
{
  int w = gd.width  > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.width;
  int h = gd.height > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.height;
  switch (gd.size)
  {
    case sizeof (uint8_t) :
    case sizeof (int) :
    case sizeof (double) :
      return ppm_typed (g, gd.size, out, w, h, cmap);
    default :
      return -1;
  }
}

/* distinct hue per label (golden ratio walk on the hue circle) */
static void label_color (double v, uint8_t rgb [3])
{
  if (isnan (v))
  {
    rgb [0] = nan_rgb [0], rgb [1] = nan_rgb [1], rgb [2] = nan_rgb [2];
    return;
  }
  long long label = llround (v);
  if (label == 0)
  {
    rgb [0] = rgb [1] = rgb [2] = 0;
    return;
  }
  double hue = fmod ((double) label * 0.618033988749895, 1.);
  hue = (hue < 0. ? hue + 1. : hue) * 6.;
  const double s = 0.65, val = 0.95;
  int i = (int) hue;
  double f = hue - i;
  double p = val * (1 - s), q = val * (1 - s * f), r = val * (1 - s * (1 - f));
  double c [6][3] = { {val, r, p}, {q, val, p}, {p, val, r},
                      {p, q, val}, {r, p, val}, {val, p, q} };
  for (int k = 0; k < 3; ++k)
    rgb [k] = (uint8_t) lround (255. * c [i % 6][k]);
}

int ppm_labels_general (void * g, const char * out, GridData gd)
{
  int w = gd.width  > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.width;
  int h = gd.height > IMAGE_MAX_SIZE ? IMAGE_MAX_SIZE : gd.height;
  if (!image_check (g, out, w, h, 1) ||
      (gd.size != sizeof (uint8_t) && gd.size != sizeof (int) && gd.size != sizeof (double)))
    return -1;

  uint8_t * img = malloc ((size_t) w * h * 3);
  if (img == NULL)
    return -1;
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
      label_color (image_value (g, gd.size, x, y), img + 3 * ((size_t) y * w + x));
  int err = image_write (out, 3, w, h, img);
  free (img);
  return err;
}

int ppm_contour (double ** raster, const char * out, int w, int h,
  int ncontours, PPMColormap cmap)
{
  if (!image_check (raster, out, w, h, 2) || ncontours < 1)
    return -1;

  double * t = image_normalize (raster, sizeof (double), w, h);
  int * label = image_contour (raster, w, h, ncontours);
  uint8_t * img = label == NULL ? NULL : ppm_unit (t, w, h, cmap);
  int err = -1;
  if (img != NULL)
  {
    for (size_t p = 0; p < (size_t) w * h; ++p)
      if (label [p] >= 0)
        img [3*p] = img [3*p + 1] = img [3*p + 2] = 0;
    err = image_write (out, 3, w, h, img);
  }
  free (img);
  free (label);
  free (t);
  return err;
}

int ppm_hillshade (double ** dem, const char * out, int w, int h,
  double azimuth, double altitude, double cellsize, double z_factor, PPMColormap cmap)
{
  if (!image_check (dem, out, w, h, 1))
    return -1;

  double * t = image_normalize (dem, sizeof (double), w, h);
  double * shade = image_hillshade (dem, w, h, azimuth, altitude, cellsize, z_factor);
  uint8_t * img = shade == NULL ? NULL : ppm_unit (t, w, h, cmap);
  int err = -1;
  if (img != NULL)
  {
    /* keep some ambient light so that shadows don't turn black */
    for (size_t p = 0; p < (size_t) w * h; ++p)
      if (!isnan (shade [p]))
        for (int k = 0; k < 3; ++k)
          img [3*p + k] = (uint8_t) lround (img [3*p + k] * (0.3 + 0.7 * shade [p]));
    err = image_write (out, 3, w, h, img);
  }
  free (img);
  free (shade);
  free (t);
  return err;
}
