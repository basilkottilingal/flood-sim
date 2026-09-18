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

  /* coord of kozhikode 11.258753 N, 75.780411 E
  */
  double y = 11.258753, x = 75.780411;
  printf ("elevation at (%g N %g E) = %g\n",
    y, x, geotiff_elevation_at ( (double []) {x, y} ) );

  geotiff_map_destroy ();

  return 0;
}
  
/*
.. 
~~~gnuplot
.. plot this
~~~
*/
