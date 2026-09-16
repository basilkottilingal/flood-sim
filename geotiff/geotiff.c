#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

/*
.. One of the limitation of this script is the reliance on POSIX
*/
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "filemap.h"
#include "lzw.h"
#include "geotiff.h"

#define error(e)                          filemap_close_all(e)
#define is_available(cur_,end_,reqd_)                        \
  if ( (size_t) ((end_) - (cur_)) < (size_t) (reqd_) )       \
    error ("tiff file : insufficient size")

static inline
const char * jump (const char * start, const char * end, size_t dest)
{
  is_available (start, end, dest);
  return start + dest;
}

static int little_endian = 1;

static inline
uint16_t u16 (const char ** m)
{
  uint8_t * b = (uint8_t *) *m;
  *m += 2;
  return little_endian ?
    ((b[1] << 8) | b[0]) :
    ((b[0] << 8) | b[1]);
}

static inline
uint32_t u32 (const char ** m)
{
  uint8_t * b = (uint8_t *) *m;
  *m += 4;
  return little_endian ?
    ((b[3] << 24) | (b[2] << 16) | (b[1] << 8) | b[0]) :
    ((b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3]);
}

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

typedef enum  /* data format for sample (for each pixel */
{
  T_UINT = 1, /* unsigned integer data                  */
  T_INT,      /* two’s complement signed integer data   */
  T_FLOAT,    /* IEEE floating point data [IEEE]        */
  T_UNDEF     /* undefined data format                  */
} TIFFFormat;

typedef enum                          /* tiff tag types */
{
  /* tiff tags */
  TIFF_TAG_IMAGE_WIDTH        = 256,   /* SHORT or LONG */
  TIFF_TAG_IMAGE_LENGTH       = 257,   /* SHORT or LONG */
  TIFF_TAG_BITS_PER_SAMPLE    = 258,   /* SHORT         */
  TIFF_TAG_COMPRESSION        = 259,   /* SHORT         */
  TIFF_TAG_PHOTOMETRIC_INTERP = 262,   /* SHORT         */

  TIFF_TAG_STRIP_OFFSETS      = 273,   /* SHORT or LONG */
  TIFF_TAG_SAMPLES_PER_PIXEL  = 277,   /* SHORT         */
  TIFF_TAG_ROWS_PER_STRIP     = 278,   /* SHORT or LONG */
  TIFF_TAG_STRIP_BYTE_COUNTS  = 279,   /* SHORT or LONG */

  TIFF_TAG_PLANAR_CONFIG      = 284,   /* SHORT         */
  TIFF_TAG_PREDICTOR          = 317,   /* SHORT         */

  TIFF_TAG_TILE_WIDTH         = 322,   /* SHORT or LONG */
  TIFF_TAG_TILE_LENGTH        = 323,   /* SHORT or LONG */
  TIFF_TAG_TILE_OFFSETS       = 324,   /* LONG          */
  TIFF_TAG_TILE_BYTE_COUNTS   = 325,   /* SHORT or LONG */

  TIFF_TAG_SAMPLE_FORMAT      = 339,   /* SHORT         */

  /* geotiff specific */
  TIFF_TAG_GEO_PIXEL_SCALE    = 33550, /* DOUBLE        */
  TIFF_TAG_GEO_TIE_POINT      = 33922, /* DOUBLE        */
  TIFF_TAG_GEO_KEY_DIR        = 34735, /* DOUBLE        */ 
  TIFF_TAG_GEO_DOUBLE_PARAMS  = 34736, /* DOUBLE        */
  TIFF_TAG_GEO_ASCII_PARAMS   = 34737, /* ASCII         */

  /* Not a (geo)tiff tag. Used as the end of tag array   */
  TIFF_TAG_NOT_A_TAG          = 0,
    
} TIFFTag;

typedef enum           /* image tile compression details */
{

  TIFF_COMP_NONE              = 1,
  TIFF_COMP_LZW               = 5,
  TIFF_COMP_PACKBITS          = 32773,
  TIFF_COMP_DEFLATE           = 8
} TIFFCompression;

typedef enum                /* Predictor for compression */
{
  TIFF_PREDICTOR_DEFAULT      = 1,
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

typedef struct
{
  uint32_t y, x;
} coord;

typedef struct Image
{
  coord dim;          /* dimension of img in pixels               */
  coord tdim;         /* dimension of tile in pixels              */
  coord n;            /* number of tiles in each dimension        */

  uint32_t offsets_at;
  uint32_t byte_counts_at;
  uint16_t byte_counts_type;

  struct
  {
    uint16_t pred;
    uint16_t type;
  } comp;             /* compression details & decompression tools */

  struct {
    uint16_t samples; /* samples per pixel                         */
    uint16_t bits;    /* bits per sample                           */
    uint16_t format;  /* data type : unsigned, float, double, etc  */
  } pixel;

} Image;

void write_decoded_pixels (void * decoded, coord tile, Image img)
{
  char * const db  = filemap_address (FILEMAP_PIXELS);
  assert (db != NULL);
  char * const dbend = db + filemap_size (FILEMAP_PIXELS);

  uint16_t
    nbits   = img.pixel.bits,
    samples = img.pixel.samples,
    hdiff   = img.comp.pred == TIFF_PREDICTOR_HORIZONTAL_DIFFERENCING,
    format  = img.pixel.format;

  if (samples != 1)
    error ("expects only one sample per pixel in geotiff");

  uint8_t * b = (uint8_t *) decoded;

  coord
    start = (coord) { tile.y * img.tdim.y    , tile.x * img.tdim.x     },
    end   = (coord) { (tile.y+1) * img.tdim.y, (tile.x+1) * img.tdim.x },
    lim   = (coord) { end.y > img.dim.y ? img.dim.y : end.y,
                                 end.x > img.dim.x ? img.dim.x : end.x };

  if (nbits == 32)
  {
    for (unsigned int h=start.y; h<lim.y; h++) {
      is_available (db, dbend, (h * img.dim.x + lim.x) * 4);
      char * data = & db [ (h * img.dim.x + start.x) * 4];
      uint32_t prev = 0;
      for (unsigned int w=start.x; w<lim.x; w++) {
        prev += little_endian ? 
          ((b[3] << 24) | (b[2] << 16) | (b[1] << 8) | b[0]) :
          ((b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3]);
        memcpy (data, &prev, 4);
        if (!hdiff)
          prev = 0;
        data += 4;
        b += 4;
      }
      for (unsigned int w=lim.x; w<end.x; w++)
        b += 4;
    }
    return;
  }

  if (nbits == 16)
  {
    for (unsigned int h=start.y; h<lim.y; h++) {
      is_available (db, dbend, (h * img.dim.x + lim.x) * 2);
      char * data = & db [ (h * img.dim.x + start.x) * 2];
      uint16_t prev = 0;
      for (unsigned int w=start.x; w<lim.x; w++) {
        prev += little_endian ? 
          ((b[1] << 8) | b[0]) :
          ((b[0] << 8) | b[1]);
        memcpy (data, &prev, 2);
        if (!hdiff)
          prev = 0;
        data += 2;
        b += 2;
      }
      for (unsigned int w=lim.x; w<end.x; w++)
        b += 2;
    }
    return;
  }

  error ("pixel nbits not supported");
}

Image img_details (Entry * entries)
{

  int is_tiles = 0, is_strips = 0, ntiles = 0;
  Image img = {0};

  while (1)
  {

    Entry entry = *entries++;
    if (entry.tag == TIFF_TAG_NOT_A_TAG)
      break;

    switch (entry.tag)
    {

      /* Image dimension */
      case TIFF_TAG_IMAGE_WIDTH :
        img.dim.x = entry.value;
        break;
      case TIFF_TAG_IMAGE_LENGTH :
        img.dim.y = entry.value;
        break;

      /* Dimension of tile */
      case TIFF_TAG_TILE_LENGTH :
        img.tdim.y = entry.value;
        is_tiles = 1;
        break;
      case TIFF_TAG_TILE_WIDTH :
        img.tdim.x = entry.value;
        is_tiles = 1;
        break;
      case TIFF_TAG_ROWS_PER_STRIP :
        img.tdim.y = entry.value;
        is_strips = 1;
        break;

      /* Location of offset array & byte count array */
      case TIFF_TAG_STRIP_OFFSETS :
      case TIFF_TAG_TILE_OFFSETS :
        img.offsets_at = entry.value;
        break;
      case TIFF_TAG_STRIP_BYTE_COUNTS :
      case TIFF_TAG_TILE_BYTE_COUNTS :
        ntiles           = entry.count;
        img.byte_counts_at   = entry.value;
        img.byte_counts_type = entry.type; 
        break;

      /* Compression details */
      case TIFF_TAG_COMPRESSION :
        img.comp.type = entry.value;
        break;
      case TIFF_TAG_PREDICTOR :
        img.comp.pred = entry.value;
        break;

      /* information on pixel data*/
      case TIFF_TAG_SAMPLE_FORMAT :
        img.pixel.format = entry.value;
        break;
      case TIFF_TAG_SAMPLES_PER_PIXEL :
        img.pixel.samples = entry.value;
        break;
      case TIFF_TAG_BITS_PER_SAMPLE :
        img.pixel.bits = entry.value;
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
  }


  /* some error/warning checks */
  if (is_strips && is_tiles)
    error ("image tile/strips redefined");
  if (is_strips)
    img.tdim.x = img.dim.x;
  if (img.tdim.x == 0 || img.tdim.y == 0 || 
      img.offsets_at == 0 || img.byte_counts_at == 0)
    error ("tile details not defined");
  if (! (img.byte_counts_type == SHORT || img.byte_counts_type == LONG) )
    error ("expect only SHORT/LONG for BYTE_COUNT data type");
  if (img.pixel.bits == 0 || img.pixel.samples == 0 || img.pixel.format == 0)
    error ("pixel data information missing");
  if ( ! (img.comp.type == TIFF_COMP_LZW) )
    error ("implementation error : only lzw compression expected");

  img.n.x = (img.dim.x + img.tdim.x - 1) / img.tdim.x;
  img.n.y = (img.dim.y + img.tdim.y - 1) / img.tdim.y;

  printf ("image [h%6u x w%6u] pixels\n", img.dim.y, img.dim.x );
  printf ("tile  [h%6u x w%6u] pixels\n", img.tdim.y, img.tdim.x );
  printf ("grid  [ %6u x  %6u] tiles\n", img.n.y, img.n.x );
  printf ("pixel [ %6u x  %6u] samples x bytes/sample \n",
    img.pixel.samples, img.pixel.bits >> 3);
  printf ("tile count %u (does it match %u ?)\n", ntiles, img.n.x * img.n.y);
  //printf ("compression type %u\n", img.comp.type);
  printf ("location of {offsets %u, byte_counts %u}\n",
    img.offsets_at, img.byte_counts_at);
  printf ("data format %u", img.pixel.format);

  return img;
}

Image ifd_read ()
{
  const char * const start = filemap_address (FILEMAP_TIFF), * m = start;
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  /*
  .. read tiff metadata. check at least 8 bytes available
  .. 2 : MM/II
  .. 2 : tiff magic number
  .. 4 : first ifd loc
  */
  is_available (m, end, 8);
  little_endian = 1;
  if (m[0] == 'M' && m[1] == 'M')
    little_endian = 0;
  else if ( !(m [0] == 'I' && m [1] == 'I') )
    error ("not a tiff file");
  m += 2;
  
  if (u16 (&m) != 42)
    error ("unexpected tiff magic number\n");

  /*
  .. We read only the first IFD and expect only one IFD.
  .. jump to the offset of the first ifd
  */
  uint32_t ifd_loc = u32 (&m);
  m = jump (start, end, ifd_loc);

  is_available (m, end, 2);
  uint16_t nentries = u16 (&m);
  if (!nentries)
    error ("ifd entries not found");

  printf (
    "no: entries %u\n"
    "%6s %4s %8s %16s\n",
    (unsigned) nentries,
    "Tag", "Type", "Count", "Value/Offset"
  );
   
  Entry * entries = malloc ((1 + nentries) * sizeof (Entry));
  is_available (m, end, nentries * 12);
  for (uint16_t n = 0; n < nentries; ++n)
  {
    uint16_t tag   = u16 (&m);
    uint16_t type  = u16 (&m);
    uint32_t count = u32 (&m);
    uint32_t value = u32 (&m);

    printf ( "%6u %4u %8u %16u\n",
      tag, type, count, value );

    entries [n] = (Entry) { tag, type, count, value };
  }
  entries [nentries] = (Entry) {.tag = TIFF_TAG_NOT_A_TAG};

  Image img = img_details (entries);

  free (entries);

  /* next ifd */
  is_available (m, end, 4);
  uint32_t next = u32 (&m);
  if (next != 0u)
    error ("expects only one IFD");

  return img;
}

void img_unwrap (Image img)
{
  /*
  .. Mapping decoded tile data onto the 2D array requires following information
  .. (1) pixel order : Row-major (default)
  .. (2) Byte Ordering (Endianness) 
  .. (3) Sample ordering : Planar (RRR..GGG..BBB..) / Chunked (RGBRGB...)
  .. (4) Sample datatype : 1 unsigned, 2 int, 3 float, 4 Custom
  .. (5) Sample data length in no: of bits.
  .. (6) Predictor (if used before encoding) : 1 not used, 2 for horizontal diff

  size_t bytes = img.pixel.bits >> 3;
  switch (bytes) {
    case 8 :
      break;
    case 4 :
      break;
    case 2 :
      break;
    case 1 :
      break;
    default :
      error ("sample size not a power of 2");
  }

  //fixme : for the moment we assume ...
  
    for (int i = tilex [0] * img.tdim [0]; i< (tilex [0] + 1) * img.tdim [0]; ++i) {
      //char * m = & db.start [tilex.x * im
      for (int j = tilex [1] * img.tdim [1]; i< (tilex [1] + 1) * img.tdim [1]; ++i) {
        if (
      }
    }
  */

  const char * const start  = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end = start + filemap_size (FILEMAP_TIFF);

  size_t
    gridsize = img.tdim.y * img.tdim.x 
               * img.pixel.samples * (img.pixel.bits >> 3),
    imgsize  = img.dim.y * img.dim.x
               * img.pixel.samples * (img.pixel.bits >> 3);

  if (imgsize > (1<<30))
    error ("tiff image is too big. reprogramme to load image partially");
  
  unsigned char * buffer = malloc (gridsize);

  filemap_pixels (imgsize);

  is_available (start, end, img.offsets_at + img.n.y * img.n.x * 4);
  is_available (start, end,
    img.byte_counts_at + img.n.y * img.n.x *
    (img.byte_counts_type == SHORT ? 2 : 4)
  );

  const char
    * offset_at      = start + img.offsets_at,
    * byte_counts_at = start + img.byte_counts_at;

  for (unsigned int i=0; i<img.n.y; i++) 
    for (unsigned int j=0; j<img.n.x; j++)
    {
      uint32_t
        offset      = u32 (&offset_at),
        byte_counts = img.byte_counts_type == SHORT ?
                      u16 (&byte_counts_at) : u32 (&byte_counts_at);
      is_available (start, end, offset + byte_counts);

      if (img.comp.type == TIFF_COMP_LZW)
      {
        if (lzw_decode (start + offset, byte_counts, buffer, gridsize))
          error ("decoding failed");
        //lzw_error (err);
      }
      write_decoded_pixels (buffer, (coord) {i, j}, img);
    }

  free (buffer);
}

void geotiff_map (const char * tiff)
{
  filemap_tiff (tiff);
  Image img = ifd_read ();
  img_unwrap (img);

  #if 1
  float * val = (float *) filemap_address (FILEMAP_PIXELS);
  assert (val);
  uint32_t W = img.dim.x;
  //17414 x w  9756
  FILE * sample = fopen ("SampleGrid.dat", "w"); 
  for (uint32_t i=15; i<527; ++i) {
    float * row = & val [i*W + 9000];
    for (int j=0; j<512; ++j)
      fprintf (sample, "%f ", row [j]);
    fprintf (sample, "\n");
  }
  fclose (sample);
  #endif
}

void geotiff_map_destroy ()
{
  filemap_close_all (NULL);
}

#undef error
#undef is_available
