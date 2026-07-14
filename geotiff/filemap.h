#ifndef _GEOTIFF_FILEMAP_H_
#define _GEOTIFF_FILEMAP_H_

  #include <stdio.h>
  #include <stdint.h>
  #include <string.h>
  #include <stdlib.h>
  #include <errno.h>
  
  #include <fcntl.h>
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <unistd.h>

  #include "lzw-decode.h"
  
  typedef struct File
  {
    int fd;
    size_t size;
    char * data;
    char * end;
  } File;

  #define DATABASE "data.bin" 

  /* global var */
  int little_endian = 1;
  File fp = (File) {-1, 0, NULL, NULL};
  File db = (File) {-1, 0, NULL, NULL};
 
  /* exit with msg */ 
  #define error(e)                     \
    do                                 \
      {                                \
      if (errno)                       \
        perror (e);                    \
      else                             \
        fprintf (stderr, "%s", e);     \
      unmap ();                        \
      exit (-1);                       \
    } while (0)
  
  /* munmap memory mapped to file "fp' */
  void unmap ()
  {
    
    if (fp.fd != -1) {
      close (fp.fd);
      if (fp.data != NULL)
        munmap (fp.data, fp.size);
    }
    fp = (File) {-1, 0, NULL, NULL};

    if (db.fd != -1) {
      unlink (DATABASE);
      close (db.fd);
      if (db.data != NULL)
        munmap (db.data, db.size);
    }
    db = (File) {-1, 0, NULL, NULL};

  }

  char * jump (size_t loc)
  {
    if (fp.size < loc)
      error ("offset out of file size");
    return fp.data + loc;
  }

  uint16_t u16 (char ** m)
  {
    if (fp.end - *m < 2)
      return UINT16_MAX;
    uint8_t * b = (uint8_t *) *m;
    *m += 2;
    return little_endian ?
      ((b[1] << 8) | b[0]) :
      ((b[0] << 8) | b[1]);
  }

  uint32_t u32 (char ** m)
  {
    if (fp.end - *m < 4)
      return UINT32_MAX;
    uint8_t * b = (uint8_t *) *m;
    *m += 4;
    return little_endian ?
      ((b[3] << 24) | (b[2] << 16) | (b[1] << 8) | b[0]) :
      ((b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3]);
  }
  
  /* map file "f" to virtual memory using mmap */
  int tiff_map (const char * f)
  {
  
    int fd = open (f, O_RDONLY);
    if (fd == -1)
      error ("open () tiff file");
    fp.fd = fd;
  
    struct stat s;
    if (fstat (fd, &s) == -1)
      error ("fstat () tiff file");
    fp.size = s.st_size;
    if (fp.size < 8)
      error ("not a tiff file");
  
    char * data =
      mmap (NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (data == NULL)
      error ("mmap () tiff file");
  
    fp.data = data;
    fp.end = data + fp.size;

    return 0;
  }

  /* create a 2D persistent paageable array */
  int database (size_t size)
  {
    int fd = open (DATABASE, O_RDWR | O_CREAT, 0666);
    if (fd == -1)
      error ("open () database");

    db.fd = fd;
    db.size = size; 
    if (ftruncate (fd, (off_t) size) == -1)
      error ("ftruncate () database");
  
    char * data =
      mmap (NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (data == NULL)
      error ("mmap () database");
    db.data = data;
  
    return 0;
  }
  
#endif
