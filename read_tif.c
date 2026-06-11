#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

static uint16_t read_u16(FILE *fp, int little_endian)
{
  uint8_t b[2];
  fread(b, 1, 2, fp);

  if (little_endian)
    return b[0] | (b[1] << 8);
  else
    return (b[0] << 8) | b[1];
}

static uint32_t read_u32(FILE *fp, int little_endian)
{
  uint8_t b[4];
  fread(b, 1, 4, fp);

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

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf(stderr, "Usage: %s file.tif\n", argv[0]);
    return 1;
  }

  FILE *fp = fopen(argv[1], "rb");

  if (!fp)
  {
    perror("fopen");
    return 1;
  }

  char endian[2];
  fread(endian, 1, 2, fp);

  int little_endian;

  if (endian[0] == 'I' && endian[1] == 'I')
  {
    little_endian = 1;
    printf("Endian: Little (II)\n");
  }
  else if (endian[0] == 'M' && endian[1] == 'M')
  {
    little_endian = 0;
    printf("Endian: Big (MM)\n");
  }
  else
  {
    printf("Not a TIFF file\n");
    fclose(fp);
    return 1;
  }

  uint16_t magic = read_u16(fp, little_endian);

  printf("Magic Number: %u\n", magic);

  if (magic != 42)
  {
    printf("Unexpected TIFF magic number\n");
  }

  uint32_t ifd_offset = read_u32(fp, little_endian);

  printf("First IFD Offset: %u\n", ifd_offset);

  fseek(fp, ifd_offset, SEEK_SET);

  uint16_t num_entries = read_u16(fp, little_endian);

  printf("IFD Entries: %u\n\n", num_entries);

  for (uint16_t i = 0; i < num_entries; i++)
  {
    uint16_t tag  = read_u16(fp, little_endian);
    uint16_t type = read_u16(fp, little_endian);
    uint32_t count = read_u32(fp, little_endian);
    uint32_t value = read_u32(fp, little_endian);

    printf(
      "Entry %3u: "
      "Tag=%5u "
      "Type=%2u "
      "Count=%8u "
      "Value/Offset=%10u\n",
      i,
      tag,
      type,
      count,
      value
    );
  }

  fclose(fp);
  return 0;
}
