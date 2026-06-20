/*
.. Display Metadata of the GEOTIFF file (topography_ESA_Copernicus_30m_resolution.tif)
.. gcc -O2 -o run tiff_meta.c && ./run topography_ESA_Copernicus_30m_resolution.tif
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

/*
.. short and long 
*/ 

#define nan    UINT32_MAX
int little_endian = 1;

static uint16_t read_u16 (FILE *fp)
{
  uint8_t b[2];
  if (fread (b, 1, 2, fp) < 2)
    return UINT16_MAX;

  if (little_endian)
    return b[0] | (b[1] << 8);
  else
    return (b[0] << 8) | b[1];
}

static uint32_t read_u32 (FILE *fp)
{
  uint8_t b[4];
  if (fread (b, 1, 4, fp) < 4)
    return nan;

  if (little_endian)
    return b[0]
       | (b[1] << 8)
       | (b[2] << 16)
       | (b[3] << 24);
  else
    return (b[0] << 24)
       | (b[1] << 16)
       | (b[2] << 8)
       | b[3];
}

/*
.. tiff data types
*/
typedef enum {
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

/*
.. tiff tag types
*/
typedef enum {
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

  /*
  .. geotiff specific
  */
  TIFF_TAG_GEO_PIXEL_SCALE    = 33550, /* double [3] */
  TIFF_TAG_GEO_TIE_POINT      = 33922, /* double [6] */
  TIFF_TAG_GEO_KEY_DIR        = 34735, /* short [32] */ 
  TIFF_TAG_GEO_DOUBLE_PARAMS  = 34736, /* double [2] */
  TIFF_TAG_GEO_ASCII_PARAMS   = 34737, /* ascii [*] */
    
} TIFFTag;

typedef enum {
  TIFF_COMP_NONE      = 1,
  TIFF_COMP_LZW       = 5,
  TIFF_COMP_PACKBITS  = 32773,
  TIFF_COMP_DEFLATE   = 8
} TIFFCompression;

struct img {
  uint16_t tag, type;
  uint32_t count, value;
};

struct meta {
  uint32_t h, w,           /* height x width             */
    th, tw,                /* tile height x weight       */
    nt,                    /* no: of tiles               */
    offsets_loc,           /* offset of offsets of tiles */
    bytes_loc;             /* offset of byte counts of tiles */
  uint32_t compression;    /* compression type           */
  /* etc.. */
};

struct meta data;

int img_details (struct img * imgs, uint16_t n) {

  struct img * i = imgs;

  data.h = data.w = data.th = data.tw = data.nt = 
    data.compression = nan;

  while (n--) {
    switch (i->tag) {

  case TIFF_TAG_IMAGE_WIDTH           :   /* SHORT or LONG */
    data.w = i->value;
    break;
  case TIFF_TAG_IMAGE_LENGTH          :   /* SHORT or LONG */
    data.h = i->value;
    break;
  case TIFF_TAG_BITS_PER_SAMPLE       :   /* SHORT */
    break;
  case TIFF_TAG_COMPRESSION           :   /* SHORT */
    data.compression = i->value;
    break;
  case TIFF_TAG_PHOTOMETRIC_INTERP    :   /* SHORT */
    break;
  case TIFF_TAG_STRIP_OFFSETS         :   /* SHORT or LONG */
    break;
  case TIFF_TAG_SAMPLES_PER_PIXEL     :   /* SHORT */
    break;
  case TIFF_TAG_ROWS_PER_STRIP        :   /* SHORT or LONG */
    break;
  case TIFF_TAG_STRIP_BYTE_COUNTS     :   /* SHORT or LONG */
    break;
  case TIFF_TAG_PLANAR_CONFIG         :   /* SHORT */
    break;
  case TIFF_TAG_PREDICTOR             :   /* SHORT */
    break;
  case TIFF_TAG_TILE_WIDTH            :   /* SHORT or LONG */
    data.tw = i->value;
    break;
  case TIFF_TAG_TILE_LENGTH           :   /* SHORT or LONG */
    data.th = i->value;
    break;
  case TIFF_TAG_TILE_OFFSETS          :   /* LONG */
    data.nt = i->count;
    data.offsets_loc = i->value;
    break;
  case TIFF_TAG_TILE_BYTE_COUNTS      :   /* SHORT or LONG */
    data.bytes_loc = i->value;
    break;
  case TIFF_TAG_SAMPLE_FORMAT         :   /* SHORT */
    break;
  /*
  .. geotiff specific
  */
    break;
  case TIFF_TAG_GEO_PIXEL_SCALE       : /* DOUBLE */
    break;
  case TIFF_TAG_GEO_TIE_POINT         : /* DOUBLE */
    break;
  case TIFF_TAG_GEO_KEY_DIR           : /* SHORT */ 
    break;
  case TIFF_TAG_GEO_DOUBLE_PARAMS     : /* DOUBLE */
    break;
  case TIFF_TAG_GEO_ASCII_PARAMS      : /* ASCII */
    break;
  default :
    }
    i++;
  }

  if (data.w == nan || data.h == nan) {
    fprintf (stderr, "\nerror : img dimensions missing");
    return -1;
  }
  printf ("\nimage [h%u x w%u] pixels", data.h, data.w);
  if (data.nt == nan || data.tw == nan || data.th == nan) {
    fprintf (stderr, "\nwarning : not stored as tiles ??");
    return -2;
  }
  printf ("\ntile [h%u x w%u] pixels", data.th, data.tw);
  printf ("\ntile count %u (does it match %u ?)", data.nt,
    ((uint32_t) ((data.w + data.tw - 1) / data.tw)) *
    ((uint32_t) ((data.h + data.th - 1) / data.th))
    );
  printf ("\ncompression type %u", data.compression);
  printf ("\nlocation of {offsets %u, byte_lengths %u}", 
    data.offsets_loc, data.bytes_loc);

  return 0;
}




int main (int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf (stderr, "Usage: %s file.tif\n", argv[0]);
    return 1;
  }

  FILE *fp = fopen (argv[1], "rb");

  if (!fp)
  {
    perror ("fopen");
    return 1;
  }

  /*
  ..  Read endianness
  */
  char endian[2];
  if (fread (endian, 1, 2, fp) < 2)
    perror ("reading endian");

  if (endian[0] == 'I' && endian[1] == 'I')
  {
    printf ("Endian: Little (II)\n");
  }
  else if (endian[0] == 'M' && endian[1] == 'M')
  {
    little_endian = 0;
    printf ("Endian: Big (MM)\n");
  }
  else
  {
    printf ("Not a TIFF file\n");
    fclose (fp);
    return 1;
  }

  /*
  .. verify "tiff" magic number
  */
  uint16_t magic = read_u16 (fp);
  if (magic != 42)
  {
    printf ("Unexpected TIFF magic number\n");
  }

  uint32_t ifd_offset = read_u32 (fp);
  printf ("First IFD Offset: %u\n", ifd_offset);
  fseek (fp, ifd_offset, SEEK_SET);

  /*
  .. see number of image file directories (IFD)
  */
  uint16_t num_entries = read_u16 (fp);
  printf ("IFD Entries: %u\n\n", num_entries);
  struct img * imgs = malloc (num_entries * sizeof (struct img));

  struct img * i = imgs;
  for (uint16_t j = 0; j < num_entries; j++)
  {
    uint16_t tag   = read_u16 (fp);
    uint16_t type  = read_u16 (fp);
    uint32_t count = read_u32 (fp);
    uint32_t value = read_u32 (fp);

    printf (
      "Entry %3u: "
      "Tag %6u "
      "Type %2u "
      "Count %8u "
      "Value/Offset %10u\n",
      j,
      tag,
      type,
      count,
      value
    );

    *i = (struct img) {
      tag, type, count, value
    };
    i++;
  }

  img_details (imgs, num_entries);

  /*
  .. Next IFD entry's offset
  */
  assert (read_u32 (fp) == 0u);

  free (imgs);
  fclose (fp);
  return 0;
}
