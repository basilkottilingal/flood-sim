#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "filemap.h"

typedef struct
{
  int fd;
  size_t size;
  char * data;
} Map;

#define unmapped   (Map) {-1, 0, NULL}
#define database   "data.bin"
#define error(e)   filemap_close_all (e)
#define unmap(f)                      \
  do                                  \
  {                                   \
    if (f.fd != -1)                   \
    {                                 \
      close (f.fd);                   \
      if (f.data != NULL)             \
        munmap (f.data, f.size);      \
      f = unmapped;                   \
    }                                 \
  } while (0)

static Map fp = unmapped;
static Map db = unmapped;

char * filemap_address (int m)
{
  return
    m == FILEMAP_TIFF   ? fp.fd < 0 ? NULL : fp.data : 
    m == FILEMAP_PIXELS ? db.fd < 0 ? NULL : db.data :
    NULL;
}

size_t filemap_size (int m)
{
  return
    m == FILEMAP_TIFF   ? fp.fd < 0 ? 0 : fp.size : 
    m == FILEMAP_PIXELS ? db.fd < 0 ? 0 : db.size :
    0;
}

/* munmap memory mapped to file "fp' and 'db'*/
void filemap_close_all (const char * error_if_any)
{
  unmap (fp);
  unmap (db);
  if (error_if_any != NULL)
  {
    if (errno)
      perror (error_if_any);
    else
      fprintf (stderr, "%s", error_if_any);
    fflush (stderr);
    exit (-1);
  }
}

/* map file "f" to virtual memory using mmap */
int filemap_tiff (const char * tiff)
{
  if (fp.fd != -1)
    error ("another tiff file mapped");

  int fd = open (tiff, O_RDONLY);
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

  return 0;
}

/* create a 2D persistent paageable array */
int filemap_pixels (size_t size)
{
  if (db.fd != -1)
    error ("database exists");

  int fd = open (database, O_RDWR | O_CREAT, 0666);
  if (fd == -1)
    error ("open () database");

  db.fd   = fd;
  db.size = size;
  if (ftruncate (fd, (off_t) size) == -1)
    error ("ftruncate () database");

  char * data =
    mmap (NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (data == NULL)
    error ("mmap () database");
  db.data = data;

  /* this will make sure the file is not saved */
  unlink (database);
  return 0;
}

#undef error
#undef database
#undef unmap
#undef unmapped
