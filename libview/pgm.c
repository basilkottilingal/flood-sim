#include "geotiff.h"
#include "pgm.h"

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>

int pgm_grayscale (DEM dem, const char * out, int width, int height)
{
  if (width < 0 || height < 0 || out == NULL || out [0] == '\0')
    return -1;
  width  = width  > dem.n ? dem.n > 1024 ? 1024 : dem.n : width;
  height = height > dem.m ? dem.m > 1024 ? 1024 : dem.n : height;

  /* Find range */
  float min = FLT_MAX, max = FLT_MIN;
  float ** raster = dem.raster;
  assert (raster != NULL);

  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
    {
      if (raster [y][x] < min)
        min = raster [y][x];
      if (raster [y][x] > max)
        max = raster [y][x];
    }

  FILE *fp = fopen(out, "wb");
  if (fp == NULL)
    return -1;

  /* PGM header */
  fprintf(fp, "P5\n%d %d\n255\n", width, height);

  /* Convert to grayscale */
  if (max == min)
  {
    float pixel = 0.0f;
    for (int i=0; i<width*height; ++i)
      fwrite(&pixel, 1, 1, fp);
    fclose (fp);
    return 0;
  }

  float den = max - min;
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
    {
      float v = (raster [y][x] - min) / den;
      unsigned char pixel = (unsigned char)(v * 255.0);
      fwrite(&pixel, 1, 1, fp);
    }

  fclose(fp);
  return 0;  
}
