/*
.. Convert geotiff file to 2d array of points (mmap() is used to avoid RAM overload)
.. gcc -O2 -o run geotiff.c && ./run topography_ESA_Copernicus_30m_resolution.tif
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

/*
.. One of the limitation of this script is the reliance on POSIX
*/
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

/*
.. mmap () virtual memory to file
*/
struct file {
  int fd;
  size_t size;
  char * data;
  char * end;
};

struct file fp = (struct file) {-1, 0, NULL, NULL};

void unmap () {
  if (fp.fd == -1)
    return;
  close (fd);
  if (fp.data != NULL)
    munmap (fp.data, fp.size);
}
  
#define error(e) do {                \
    if (errno)                       \
      perror (e);                    \
    else                             \
      fprintf (stderr, "%s", e);     \
    unmap ();                        \
    exit (-1);                       \
  } while (0)

/*
.. Map a file to a contiguous virtual memory
*/
int map (const char * f) {

  int fd = open (f, O_RDONLY);
  if (fd == -1)
    error ("open ()");
  fp.fd = fd;

  struct stat s;
  if (fstat (fd, &s) == -1)
    error ("fstat ()");
  fp.size = s.st_size;

  char * data = mmap (NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (data == NULL)
    error ("mmap ()");
  fp.data = data;
  fp.end = data + fp.size;

  return 0;
}

#define nan    UINT32_MAX
int little_endian = 1;

static uint16_t u16 (char ** m)
{
  if (fp.end - *m < 2)
    return UINT16_MAX;
  uint8_t * b = (uint8_t *) *m;
  *m += 2;
  return little_endian ?
    ((b[1] << 8) | b[0]) :
    ((b[0] << 8) | b[1]);
}

static uint32_t u32 (char ** m)
{
  if (fp.end - *m < 4)
    return nan;
  uint8_t * b = (uint8_t *) *m;
  *m += 4;
  return little_endian ?
    ((b[3] << 24) | (b[2] << 16) | (b[1] << 8) | b[0]) :
    ((b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3]);
}

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
} TIFFType;             /* tiff datatype {1, 2, .., 12}*/

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
  TIFF_TAG_GEO_PIXEL_SCALE    = 33550,         /* DOUBLE */
  TIFF_TAG_GEO_TIE_POINT      = 33922,         /* DOUBLE */
  TIFF_TAG_GEO_KEY_DIR        = 34735,         /* DOUBLE */ 
  TIFF_TAG_GEO_DOUBLE_PARAMS  = 34736,         /* DOUBLE */
  TIFF_TAG_GEO_ASCII_PARAMS   = 34737,          /* ASCII */
    
} TIFFTag;                             /* tiff tag types */

typedef enum {
  TIFF_COMP_NONE      = 1,
  TIFF_COMP_LZW       = 5,
  TIFF_COMP_PACKBITS  = 32773,
  TIFF_COMP_DEFLATE   = 8
} TIFFCompression;

struct meta {
  uint32_t h, w,           /* height x width             */
    th, tw,                /* tile height x weight       */
    nt,                    /* no: of tiles               */
    offsets_loc,           /* offset of offsets of tiles */
    bytes_loc;             /* offset of byte counts of tiles */
  uint32_t compression;    /* compression type           */
  /* etc.. */
};

struct entry {             /* IFD entry */
  uint16_t tag, type; 
  uint32_t count, value;
};

struct ifd {               /* IFD : image file directory */
  struct entry * entries;
  int n;
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

  map (argv [1]);
  char * m = fp.data;

  if (fp.end - m < 2)
    error ("not a tiff file");
  assert (fp.end - m >= 2);

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

  uint32_t ifd_offset = u32 (&m);     /* offset of first ifd */
  printf ("First IFD Offset: %u\n", ifd_offset);
  m = fp.data + ifd_offset;

  uint16_t nentries = u16 (&m);  /* number of entries in this ifd */ 
  printf ("IFD Entries: %u\n\n", nentries);
  struct img * imgs = malloc (nentries * sizeof (struct img));

  struct img * i = imgs;
  for (uint16_t j = 0; j < nentries; j++)
  {
    uint16_t tag   = u16 (&m);
    uint16_t type  = u16 (&m);
    uint32_t count = u32 (&m);
    uint32_t value = u32 (&m);

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

    *i++ = (struct entry) { tag, type, count, value };
  }

  img_details (imgs, nentries);

  /*
  .. Next IFD entry's offset
  */
  ifd_offset = u32 (fp);
  printf ("next IFD's offset %u", ifd_offset);

  free (imgs);
  fclose (fp);
  return 0;
}
