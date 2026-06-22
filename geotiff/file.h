#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct File
{
  int fd;
  size_t size;
  char * data;
  char * end;
} File;

File fp = (File) {-1, 0, NULL, NULL};

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

void unmap ()
{
  if (fp.fd == -1)
    return;
  close (fp.fd);
  if (fp.data != NULL)
    munmap (fp.data, fp.size);
}

int map (const char * f)
{

  int fd = open (f, O_RDONLY);
  if (fd == -1)
    error ("open ()");
  fp.fd = fd;

  struct stat s;
  if (fstat (fd, &s) == -1)
    error ("fstat ()");
  fp.size = s.st_size;
  if (fp.size < 8)
    error ("not a tiff file");

  char * data =
    mmap (NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (data == NULL)
    error ("mmap ()");

  fp.data = data;
  fp.end = data + fp.size;

  return 0;
}

