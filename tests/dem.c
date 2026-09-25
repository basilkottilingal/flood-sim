/*
.. Convert geotiff file to 2d array of points (mmap() is used to avoid RAM overload)
.. # 
.. gcc -D_GNU_SOURCE -O2 -o run geotiff.c && ./run ./../topography_ESA_Copernicus_30m_resolution.tif
.. # strict POSIX compliance
.. gcc -D_XOPEN_SOURCE=700 -O2 -o run geotiff.c && ./run ./../topography_ESA_Copernicus_30m_resolution.tif
*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

/*
.. One of the limitation of this script is the reliance on POSIX
*/
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
  
#ifndef MADV_DONTNEED        /* In case of strict posix compliance */
  #define MADV_DONTNEED POSIX_MADV_DONTNEED
  #define madvise(a,b,c) posix_madvise (a,b,c)
#endif

/*
.. related to handling (geotiff) file like open (), mmap () to 
.. virtual reading
*/

#include "filemap.h"
#include "geotiff.h"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf (stderr, "Usage: %s file.tif\n", argv[0]);
    return 1;
  }

  geotiff_map (argv [1]);

  /*
  .. kozhikode 11.258753 N, 75.780411 E
  */
  CoordG c = (CoordG) {.lat = 11.258753, .lon = 75.780411};
  DEM dem;
  if (geotiff_dem_window (c.lon, c.lat, 0., 0., &dem))
    filemap_close_all ("dem window failed");

  printf ("coord %g %g\n", c.lat, c.lon);
  printf ("DEM  tiepoint pixel (%g %g), coord (%g %g)\n",
   dem.tiepoint [0], dem.tiepoint [1], 
   dem.tiepoint [4], dem.tiepoint [3]);
  printf ("elevation grid\n dem.raster[-2:%d][-2:%d]\n", dem.n+1, dem.m+1);
  for (int j=-2; j<=dem.m+1; ++j)
  {
    printf ("\t");
    for (int i=-2; i<=dem.n+1; ++i)
      printf ("%#4.5g ", (double) dem.raster [i][j]);
    printf ("\n");
  }

  geotiff_dem_free (dem);

  geotiff_map_destroy ();

  return 0;
}
  
/*
.. 
~~~gnuplot
.. plot this
~~~
*/
