#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

/*
.. One of the limitation of this script is the reliance on POSIX
*/
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "filemap.h"
#include "lzw.h"
#include "coordinates.h"
#include "tiff.h"
#include "geotiff.h"
#include "parse.h"

#define error(e)                          filemap_close_all(e)
#define is_available(cur_,end_,reqd_)                        \
  if ( (size_t) ((end_) - (cur_)) < (size_t) (reqd_) )       \
    error ("tiff file : insufficient size")

struct
{
  const char * address;
  const char * end;
  size_t size;
} map;

static inline
const char * jump (const char * start, const char * end, size_t dest)
{
  is_available (start, end, dest);
  return start + dest;
}

#define entry_size(entry) (entry->count * tiff_datasize (entry->type))

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
    double tiepoint [6]; /* tie point. pixel +  refcoord         */
    double scale    [3]; /* scale (degrees per pixel)            */
  } coordMap;

  struct
  {
    uint16_t pred;
    uint16_t type;
  } comp;             /* compression details & decompression tools */

  struct
  {
    uint16_t samples; /* samples per pixel                         */
    uint16_t bits;    /* bits per sample                           */
    uint16_t format;  /* data type : unsigned, float, double, etc  */
  } pixel;

} Image;

/*
static void geo_tag_entries (TIFFEntry * entries)
{
  TIFFEntry
    KeyDirectry  = {.tag = TIFF_TAG_NOT_A_TAG},
    DoubleParams = {.tag = TIFF_TAG_NOT_A_TAG},
    AsciiParams  = {.tag = TIFF_TAG_NOT_A_TAG};
}
*/

static void read_geo_key_dir (Image * img, TIFFEntry * entry)
{
  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == GeoKeyDirectoryTag);
  assert (entry->count % 4 == 0);
  assert (entry->type == SHORT);

  is_available (start, end, entry->value + entry->count * 2);

  const char * array = start + entry->value;
  uint16_t
    KeyDirectoryVersion = u16 (&array),
    KeyRevision         = u16 (&array),
    MinorRevision       = u16 (&array),
    NumberOfKeys        = u16 (&array);
  printf ("KeyDirectoryVersion %u, KeyRevision %u, MinorRevision %u, NumberOfKeys %u\n",
    KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys);

  for (int i=0; i<NumberOfKeys; ++i)
  {
    geotiff_key ( (GeoKey)
      {
        u16 (&array),
        u16 (&array),
        u16 (&array),
        u16 (&array)
      }
    );

  }
}
        
static
void read_geo_pixel_scale (Image * img, TIFFEntry * entry)
{
  /*
  .. scaling of longitude and latitude per pixel.
  */

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == ModelPixelScaleTag);
  assert (entry->count == 3);
  assert (entry->type == DOUBLE);
   
  is_available (start, end, entry->value + entry->count * 8);

  const char * array = start + entry->value;
  printf ("pixel scale\n");
  double scale [3] = {d64 (&array), d64 (&array), d64 (&array)};
  printf ("\t(%g, %g, %g) degrees per pixel\n", scale [0], scale [1], scale [2]);
  memcpy (img->coordMap.scale, scale, sizeof (scale));

}

static
void read_geo_tie_point (Image * img, TIFFEntry * entry)
{

  /*
  .. Tie point(s) is a tuple of 6 double numbers.
  .. (I, J, K) represent which pixel corresponds to reference coordinate's origin
  .. (X, Y, Z) represent longitude (-180 deg W, 180 deg E],
  .. latitude [-90 deg N, 90 deg N], and elevation of the reference pixel.
  .. NOTE :
  .. (a) the reference pixel (I, J, K) can be fractional. 
  .. (b) There can be multiple tie points. 

  */

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag == ModelTiepointTag);
  assert (entry->count % 6 == 0);
  assert (entry->type == DOUBLE);
   
  is_available (start, end, entry->value + entry->count * 8);

  const char * array = start + entry->value;
  printf ("tie point\n");

  double tiepoint [6] = 
    { 
      d64 (&array), d64 (&array), d64 (&array),
      d64 (&array), d64 (&array), d64 (&array)
    };

  printf ("\ttiepoint pixel [%g, %g, %g]\n", tiepoint [0], tiepoint [1], tiepoint [2]);
  printf ("\ttiepoint coord (%g, %g, %g)\n", tiepoint [3], tiepoint [4], tiepoint [5]);
  printf ("\tInference (longitude %g latitude %g)\n", tiepoint [3], tiepoint [4]);

  if (entry->count > 6)
    error ("warning : library not designed for multiple tie points");

  memcpy (img->coordMap.tiepoint, tiepoint, sizeof (tiepoint));
}

static
void read_geo_ascii_params (Image * img, TIFFEntry * entry)
{

  const char * const start = filemap_address (FILEMAP_TIFF);
  assert (start != NULL);
  const char * const end   = start + filemap_size (FILEMAP_TIFF);

  assert (entry->tag  == GeoAsciiParamsTag);
  assert (entry->type == ASCII);
   
  is_available (start, end, entry->value + entry->count);

  const char * params = start + entry->value;
  printf ("geo ascii params\n\t%s\n", params);

}

static
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

  const char * b = (const char *) decoded;

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
        prev += u32 (&b);
        memcpy (data, &prev, 4);
        if (!hdiff)
          prev = 0;
        data += 4;
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
        prev += u16 (&b);
        memcpy (data, &prev, 2);
        if (!hdiff)
          prev = 0;
        data += 2;
      }
      for (unsigned int w=lim.x; w<end.x; w++)
        b += 2;
    }
    return;
  }

  error ("pixel nbits not supported");
}

static
Image img_details (TIFFEntry * entries)
{

  int is_tiles = 0, is_strips = 0, ntiles = 0;
  int found_geo_tags = 0;
  Image img = {0};

  while (1)
  {

    TIFFEntry entry = *entries++;
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
      case GeoKeyDirectoryTag:
        read_geo_key_dir (&img, &entry);
        found_geo_tags |= 1;
        break;
      case GeoDoubleParamsTag:
        found_geo_tags |= 2;
        break;
      case GeoAsciiParamsTag:
        read_geo_ascii_params (&img, &entry);
        found_geo_tags |= 4;
        break;
      /*
      .. there are two models for raster->coord mapping
      .. (a) pixel scale + tie point
      .. (b) matrix transformation model
      */
      case ModelPixelScaleTag :
        read_geo_pixel_scale (&img, &entry);
        found_geo_tags |= 8;
        break;
      case ModelTiepointTag :
        read_geo_tie_point (&img, &entry);
        found_geo_tags |= 16;
        break;
      case ModelTransformationTag :
        error ("implemetation error : model transformation");
        found_geo_tags |= 32;
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
  if (img.tdim.x     == 0 || img.tdim.y         == 0 || 
      img.offsets_at == 0 || img.byte_counts_at == 0 )
    error ("tile details not defined");
  if (! (img.byte_counts_type == SHORT || img.byte_counts_type == LONG) )
    error ("expect only SHORT/LONG for BYTE_COUNT data type");
  if (img.pixel.bits == 0 || img.pixel.samples == 0 || img.pixel.format == 0)
    error ("pixel data information missing");
  if ( ! (img.comp.type == TIFF_COMP_LZW) )
    error ("implementation error : only lzw compression expected");
  if ( ! (found_geo_tags == (1|2|4|8|16) || found_geo_tags == (1|2|4|32) ) )
    error ("some geo tags missing");

  img.n.x = (img.dim.x + img.tdim.x - 1) / img.tdim.x;
  img.n.y = (img.dim.y + img.tdim.y - 1) / img.tdim.y;

  printf ("image [h%6u x w%6u] pixels\n", img.dim.y , img.dim.x  );
  printf ("tile  [h%6u x w%6u] pixels\n", img.tdim.y, img.tdim.x );
  printf ("grid  [ %6u x  %6u] tiles \n", img.n.y   , img.n.x    );
  printf ("pixel [ %6u x  %6u] samples x bytes/sample \n",
    img.pixel.samples, img.pixel.bits >> 3);
  printf ("tile count %u (does it match %u ?)\n", ntiles, img.n.x * img.n.y);
  printf ("compression type : %s\n", img.comp.type == TIFF_COMP_LZW ? 
    "lzw" : "unknown" );
  printf ("location of {offsets %u, byte_counts %u}\n",
    img.offsets_at, img.byte_counts_at);
  printf ("data format %u\n", img.pixel.format);

  printf ("coord bound\n");
  double
    * tiepixel = img.coordMap.tiepoint,
    * tiecoord = &img.coordMap.tiepoint [3],
    * scale    = img.coordMap.scale;
  printf ("(%#3.4gE %#3.4gN)     (%#3.4gE %#3.4gN)\n",
    - tiepixel[0] * scale[0] + tiecoord [0],
      tiepixel[1] * scale[1] + tiecoord [1],
     ((double) img.dim.x - tiepixel[0]) * scale[0] + tiecoord [0],
      tiepixel[1] * scale[1] + tiecoord [1]);

  int tieApprox [2] =
    {
      (int) (10 * tiepixel [0] / (double) img.dim.x),
      (int) ( 5 * tiepixel [1] / (double) img.dim.y),
    };
  for (int iy = 0; iy <=5; ++iy)
  {
    printf ("               ");
    for (int ix = 0; ix <=10; ++ix)
      printf ("%c",
        ix == tieApprox [0] && iy == tieApprox [1] ? '+' :
        ix == 0 || ix == 10 || iy == 0 || iy == 5  ? '.' : ' ');
    if (iy == 3)
      printf ("       +  tiepoint");
    printf ("\n");
  }
    
  printf ("(%#3.4gE %#3.4gN)     (%#3.4gE %#3.4gN)\n",
    -tiepixel[0] * scale[0] + tiecoord [0],
    (tiepixel[1] - (double) img.dim.y) * scale[1] + tiecoord [1],
    ((double) img.dim.x - tiepixel[0]) * scale[0] + tiecoord [0],
    (tiepixel[1] - (double) img.dim.y) * scale[1] + tiecoord [1]);

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
  if (m [0] == 'I' && m [1] == 'I')
    parser_endianness (LITTLE_ENDIAN);
  else if (m[0] == 'M' && m[1] == 'M')
    parser_endianness (BIG_ENDIAN);
  else
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
   
  TIFFEntry * entries = malloc ((1 + nentries) * sizeof (TIFFEntry));
  is_available (m, end, nentries * 12);

  for (uint16_t n = 0; n < nentries; ++n)
  {
    uint16_t tag   = u16 (&m);
    uint16_t type  = u16 (&m);
    uint32_t count = u32 (&m);
    uint32_t value = u32 (&m);

    printf ( "%6u %4u %8u %16u\n",
      tag, type, count, value );

    entries [n] = (TIFFEntry) { tag, type, count, value };
  }
  entries [nentries] = (TIFFEntry) {.tag = TIFF_TAG_NOT_A_TAG};

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

static int is_running = 0;
Image geotiff = (Image) {0};

#if 0
static
int geotiff_getval (double val [2][2], int i, int j)
{
  char * address = filemap_address (FILEMAP_PIXELS);
  assert (address != NULL);
  if (geotiff.pixel.format == SHORT)
  {
    uint16_t * elevation = (uint16_t *) address, cpy;
    for (int ii=0; ii<2; ++ii)
      for (int jj=0; jj<2; ++jj)
      {
        memcpy (& cpy, & elevation [geotiff.dim.x * (i + ii) + ( j + jj)], 2);
        val [ii][jj] = (double) cpy;
      }

    return 0;
  }
  return  1;
}
#endif

double geotiff_elevation (coord c)
{
  if (geotiff.pixel.bits == 32)
  {
    if (geotiff.pixel.format == T_UINT)
      return (double)
        ((uint32_t *) filemap_address (FILEMAP_PIXELS)) [c.y * geotiff.dim.x + c.x]; 
    if (geotiff.pixel.format == T_FLOAT)
      return (double)
        ((float *)    filemap_address (FILEMAP_PIXELS)) [c.y * geotiff.dim.x + c.x]; 
  }
  if (geotiff.pixel.bits == 16)
  {
    if (geotiff.pixel.format == T_UINT)
      return (double)
        ((uint16_t *) filemap_address (FILEMAP_PIXELS)) [c.y * geotiff.dim.x + c.x]; 
  }
  /* we have skipped T_INT as elevation are usually stored as float/uint */
  error ("pixel format not implemented");
  return NAN;
}

double geotiff_elevation_at (double geo_coord [])
{
  if (!is_running)
    return NAN;
  double
    * tie   = geotiff.coordMap.tiepoint,
    * scale = geotiff.coordMap.scale,
    x       =   (geo_coord [0] - tie [3]) / scale [0]  + tie [0],
    y       = - (geo_coord [1] - tie [4]) / scale [1]  + tie [1];
  int i = floor (x), j = floor (y);
  if ( i < 0 || i >= geotiff.dim.x || j < 0 || j >= geotiff.dim.y )
    return NAN;

  /* fixme : interpolate */
  return geotiff_elevation ( (coord) {j, i} );

  return NAN;
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

  is_running = 1;
  geotiff = img;
}

void geotiff_map_destroy ()
{
  filemap_close_all (NULL);
  is_running = 0;
  geotiff = (Image) {0};
}

#undef error
#undef is_available
