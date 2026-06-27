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
#include "lzw-decode.h"
#include "tile.h"


/*
.. realated to TIFF
*/
typedef enum              /* tiff datatype {1, 2, .., 12}*/
{
  BYTE = 1,
  ASCII,
  SHORT,
  LONG,
  RATIONAL,
  SBYTE,
  UNDEFINED,
  SSHORT,
  SLONG,
  SRATIONAL,
  FLOAT,
  DOUBLE
} TIFFType;

typedef enum                          /* tiff tag types */
{
  TIFF_TAG_IMAGE_WIDTH        = 256,   /* SHORT or LONG */
  TIFF_TAG_IMAGE_LENGTH       = 257,   /* SHORT or LONG */
  TIFF_TAG_BITS_PER_SAMPLE    = 258,   /* SHORT */
  TIFF_TAG_COMPRESSION        = 259,   /* SHORT */
  TIFF_TAG_PHOTOMETRIC_INTERP = 262,   /* SHORT */

  TIFF_TAG_STRIP_OFFSETS      = 273,   /* SHORT or LONG */
  TIFF_TAG_SAMPLES_PER_PIXEL  = 277,   /* SHORT */
  TIFF_TAG_ROWS_PER_STRIP     = 278,   /* SHORT or LONG */
  TIFF_TAG_STRIP_BYTE_COUNTS  = 279,   /* SHORT or LONG */

  TIFF_TAG_PLANAR_CONFIG      = 284,   /* SHORT */
  TIFF_TAG_PREDICTOR          = 317,   /* SHORT */

  TIFF_TAG_TILE_WIDTH         = 322,   /* SHORT or LONG */
  TIFF_TAG_TILE_LENGTH        = 323,   /* SHORT or LONG */
  TIFF_TAG_TILE_OFFSETS       = 324,   /* LONG */
  TIFF_TAG_TILE_BYTE_COUNTS   = 325,   /* SHORT or LONG */

  TIFF_TAG_SAMPLE_FORMAT      = 339,   /* SHORT */

  /* geotiff specific */
  TIFF_TAG_GEO_PIXEL_SCALE    = 33550,         /* DOUBLE */
  TIFF_TAG_GEO_TIE_POINT      = 33922,         /* DOUBLE */
  TIFF_TAG_GEO_KEY_DIR        = 34735,         /* DOUBLE */ 
  TIFF_TAG_GEO_DOUBLE_PARAMS  = 34736,         /* DOUBLE */
  TIFF_TAG_GEO_ASCII_PARAMS   = 34737,          /* ASCII */
    
} TIFFTag;

typedef enum           /* image tile compression details */
{
  TIFF_COMP_NONE      = 1,
  TIFF_COMP_LZW       = 5,
  TIFF_COMP_PACKBITS  = 32773,
  TIFF_COMP_DEFLATE   = 8
} TIFFCompression;

typedef enum                /* Predictor for compression */
{
  TIFF_PREDICTOR_DEFAULT = 1,
  TIFF_PREDICTOR_HORIZONTAL_DIFFERENCING
} TIFFCompressionPredictor;


/*
.. Related to image file directory (IFD) and it's entries
*/
typedef struct Entry                        /* IFD entry */
{
  uint16_t tag, type; 
  uint32_t count, value;
} Entry;

typedef struct Ifd         /* IFD : image file directory */
{
  Entry * entries;
  int n;
  uint32_t next;
} Ifd;

Ifd ifd_read (uint32_t offset, char ** m)
{

  *m = jump (offset);
  uint16_t nentries = u16 (m);       /* number of entries */
  if (!nentries)
    return (Ifd) {0};
   
  Entry * e = malloc (nentries * sizeof (Entry)),
    * entries = e;

  printf ( "no: entries %u\n",  (unsigned) nentries);
  printf ( "     Tag"
           " Type"
           "    Count"
           " Value/Offset\n");

  for (uint16_t n = 0; n < nentries; ++n)
  {
    uint16_t tag   = u16 (m);
    uint16_t type  = u16 (m);
    uint32_t count = u32 (m);
    uint32_t value = u32 (m);

    printf ( "%6u " "%4u " "%8u " "%13u\n",
      tag,
      type,
      count,
      value
    );

    *e++ = (Entry) { tag, type, count, value };
  }

  return (Ifd) {entries, nentries, u32 (m)};
}

void ifd_free (Ifd i)
{
  if (i.entries == NULL)
    return;
  free (i.entries);
  i.entries = NULL;
  i.n = 0;
}

/*
.. Image details.
*/

enum {
  hh = 0,
  ww = 1,
} TIFFDim;


typedef struct Image
{
  uint32_t dim [2];  /* dimension of img in pixels*/
  uint32_t tdim [2]; /* dimension of tile in pixels */
  uint32_t n [2];    /* number of tiles in each dimension */
  Tile * tiles;      /* tile grid */

  struct {
    uint32_t samples; /* samples per pixel */
    uint32_t bits;    /* bits per sample    */
    uint32_t format;  /* data type : unsigned, float, double, etc */
  } pixel;

  uint32_t compression;
  uint32_t predictor;
} Image;

Image img_details (Ifd i)
{

  int nentries = i.n, is_tiles = 0, is_strips = 0, ntiles = 0;
  Entry * e = i.entries;
  Image img = {0};
  uint32_t bytes_at, offsets_at;
  int type;

  while (nentries--) {
    switch (e->tag) {

      /* Image dimension */
      case TIFF_TAG_IMAGE_WIDTH :
        img.dim [ww] = e->value;
        break;
      case TIFF_TAG_IMAGE_LENGTH :
        img.dim [hh] = e->value;
        break;

      /* Dimension of tile */
      case TIFF_TAG_TILE_WIDTH :
        img.tdim [ww] = e->value;
        is_tiles = 1;
        break;
      case TIFF_TAG_ROWS_PER_STRIP :
        img.tdim [ww] = 1;
        img.tdim [hh] = e->value;
        is_strips = 1;
        break;
      case TIFF_TAG_TILE_LENGTH :
        img.tdim [hh] = e->value;
        is_tiles = 1;
        break;

      /* Location of offset array & byte count array */
      case TIFF_TAG_STRIP_OFFSETS :
      case TIFF_TAG_TILE_OFFSETS :
        offsets_at = e->value;
        break;
      case TIFF_TAG_STRIP_BYTE_COUNTS :
      case TIFF_TAG_TILE_BYTE_COUNTS :
        bytes_at = e->value;
        ntiles = e->count;
        type = e->type == SHORT ? 0 : e->type == LONG ? 1 : -1;
        break;

      /* Compression details */
      case TIFF_TAG_COMPRESSION :
        img.compression = e->value;
        break;
      case TIFF_TAG_PREDICTOR :
        img.predictor = e->value;
        break;

      /* information on pixel data*/
      case TIFF_TAG_SAMPLE_FORMAT :
        img.pixel.format =  e->value;
        break;
      case TIFF_TAG_SAMPLES_PER_PIXEL :
        img.pixel.samples =  e->value;
        break;
      case TIFF_TAG_BITS_PER_SAMPLE :
        img.pixel.bits =  e->value;
        break;

      case TIFF_TAG_PLANAR_CONFIG :
        /*
        .. It matters if theres are more than one samples/pixel.
        .. Eg : if samples = {R, G, B} you can interleave them as
        .. RGBRGBRGB... or as RRRR..GGG...BBB....
        .. For Geotiff it is 1 sample (i.e {elevation}) per pixel.
        */
        break;

      /* geotiff specific */
      case TIFF_TAG_GEO_PIXEL_SCALE :
        break;
      case TIFF_TAG_GEO_TIE_POINT :
        break;
      case TIFF_TAG_GEO_KEY_DIR :
        break;
      case TIFF_TAG_GEO_DOUBLE_PARAMS :
        break;
      case TIFF_TAG_GEO_ASCII_PARAMS :
        break;

      /* unexpected */
      default :
    }
    e++;
  }

  if (is_strips && is_tiles)
    error ("image tile/strips redefined");
  if (img.tdim [hh] == 0 || img.tdim [ww] == 0 || 
      offsets_at == 0 || bytes_at == 0)
    error ("tile details not defined");
  if (type == -1)
    error ("expect only SHORT/LONG for byte length data type");
  if (img.pixel.bits == 0 || img.pixel.samples == 0 || img.pixel.format == 0)
    error ("pixel data information missing");

  img.n [hh] = (img.dim [hh] + img.tdim[hh] - 1) / img.tdim [hh];
  img.n [ww] = (img.dim [ww] + img.tdim[ww] - 1) / img.tdim [ww];

  printf ("image [h%6u x w%6u] pixels\n", img.dim [hh], img.dim [ww] );
  printf ("tile  [h%6u x w%6u] pixels\n", img.tdim [hh], img.tdim [ww] );
  printf ("grid  [ %6u x  %6u] tiles\n", img.n [hh], img.n [ww] );
  printf ("pixel [ %6u x  %6u] samples x bytes/sample \n", img.pixel.samples, img.pixel.bits >> 3);
  printf ("tile count %u (does it match %u ?)\n", ntiles, img.n[hh] * img.n [ww]);
  printf ("compression type %u\n", img.compression);
  printf ("location of {offsets %u, byte_lengths %u}\n", offsets_at, bytes_at);

  /* reading tile info */
  img.tiles = img_tiles (img.n [hh] * img.n [ww], offsets_at, bytes_at, type);
  if (img.tiles == NULL)
    error ("cannot create tile grid");

  return img;
}

void img_free (Image img)
{
  free (img.tiles);
  img.tiles = NULL;  
}


int main (int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf (stderr, "Usage: %s file.tif\n", argv[0]);
    return 1;
  }

  map (argv [1]);
  char * m = fp.data;

  if (m [0] == 'I' && m [1] == 'I') 
    printf ("Endian: Little (II)\n");
  else if (m[0] == 'M' && m[1] == 'M')
  {
    little_endian = 0;
    printf ("Endian: Big (MM)\n");
  }
  else
    error ("not a tiff file");
  m += 2;
  
  if (u16 (&m) != 42)
    error ("unexpected tiff magic number\n");

  Ifd i = ifd_read (u32 (&m), &m);
  if (i.entries == NULL)
    error ("ifd entries not found");

  if (i.next != 0)
    error ("expect only 1 ifd");
  
  Image img = img_details (i);

  size_t gridsize = img.tdim [0] * img.tdim [1] 
    * img.pixel.samples * (img.pixel.bits >> 3);
  unsigned char * buffer = malloc (gridsize);

  /*
  .. NOTES :
  .. 1. The edge tiles that may have smaller height and/or width than of the normal
  .. tiles, produce same lzw_decode () bytes compared to normal tiles. You have to
  .. disregard pixels outside  [ 0 : dim[hh]) x [0, dim [ww]).
  .. 2. Bytes (correspoonding to Pixels) are in row-major
  */
  Tile * tile = img.tiles;
  for (unsigned int i=0; i<img.n[0]; i++) 
    for (unsigned int j=0; j<img.n[1]; j++)
    {
      int err =
        lzw_decode (tile->start, tile->end - tile->start, buffer, gridsize);
      if (err)
        printf ("lzw_error @ tile [%u, %u] : err type %s", i, j, LZWD_ERROR(err));
    }

  /*
  .. A sample tile for visualization
  */
  do {
    /*
    FILE * fp = fopen ("sample.dat", "w");
    Tile t = img.t [10][30];
    int err = lzw_decode (t.start,t .end - t .start, buffer, gridsize);
    .. You still have to figure out
    .. 1. the byte order / unsigned (little/big endian)
    .. 2. data type is unsigned int [4] or float [4] or any other datatype
    .. 3. Whether any byte traversal order is used in compression (it's part of lzw routine)
    if (img.pixel.bits != 32) {
      error ("expected elevation data to be unsigned int [4 bytes]");
      break;
    }
    unsigned * 
    assert (err == 0);
    for (unsigned i=0; i<img.tdim [0]; ++i) {
      for (unsigned j=0; j<img.tdim [1]; ++j) {
        printf (
      }
    }
    */
  } while (0);


  free (buffer);
  img_free (img);
  ifd_free (i);
  error ("");
  return 0;
}
  
/*
.. 
~~~gnuplot
.. plot this
~~~
*/
