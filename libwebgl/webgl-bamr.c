/*
.. A .bamr file as a blocks layer: see webgl-bamr.h. A reader of its own,
.. so that libwebgl needs no libb-amr: tests/test_bamr.c checks it against
.. libb-amr's writer, bit for bit.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "webgl-bamr.h"

enum { HEAD = 80 };

static uint32_t get32 (const unsigned char * p)
{
  return (uint32_t) p[0] | (uint32_t) p[1] << 8 | (uint32_t) p[2] << 16 |
    (uint32_t) p[3] << 24;
}

static uint64_t get64 (const unsigned char * p)
{
  return (uint64_t) get32 (p) | (uint64_t) get32 (p + 4) << 32;
}

static double get_double (const unsigned char * p)
{
  uint64_t v = get64 (p);
  double x;
  memcpy (&x, &v, sizeof (x));
  return x;
}

/* bytes at offset: 0, or -1 if the file has not got them */
static int read_at (FILE * fp, uint64_t size, uint64_t offset, void * buf,
                    uint64_t bytes)
{
  if (offset > size || bytes > size - offset)
    return -1;
  if (!bytes)
    return 0;
  if (fseeko (fp, (off_t) offset, SEEK_SET) ||
      fread (buf, 1, bytes, fp) != bytes)
    return -1;
  return 0;
}

int webgl_bamr (WebglScene * S, const char * path)
{
  if (!S || !path)
    return -1;
  FILE * fp = fopen (path, "rb");
  if (!fp)
    return -1;
  int layer = -1;
  char ** names = NULL;
  double * box = NULL;
  float * level = NULL, * v = NULL;
  unsigned char * raw = NULL;
  uint32_t n = 0;

  unsigned char h[HEAD];
  if (fseeko (fp, 0, SEEK_END))
    goto done;
  uint64_t size = (uint64_t) ftello (fp);
  if (read_at (fp, size, 0, h, HEAD) || memcmp (h, "BAMRGRID", 8) ||
      get32 (h + 8) != 1 || get32 (h + 12) != 0x01020304u ||
      get32 (h + 16) != 2 || get32 (h + 72) != 0)
    goto done;
  uint32_t m = get32 (h + 20), rs = get32 (h + 24);
  n = get32 (h + 28);
  uint64_t N = get64 (h + 32);
  double o[2] = { get_double (h + 40), get_double (h + 48) };
  double h0 = get_double (h + 64);
  if (m < 1 || m > 1024 || (rs != 4 && rs != 8) || n > 4096 || N < 1 ||
      N > (1ull << 32) || !(h0 > 0.))
    goto done;

  /* the names */
  names = calloc (n ? n : 1, sizeof (char *));
  if (!names)
    goto done;
  uint64_t at = HEAD;
  for (uint32_t k = 0; k < n; k++) {
    unsigned char b[4];
    if (read_at (fp, size, at, b, 4))
      goto done;
    uint32_t len = get32 (b);
    if (len < 1 || len > 4096 || !(names[k] = malloc (len + 1)) ||
        read_at (fp, size, at + 4, names[k], len))
      goto done;
    names[k][len] = '\0';
    at += 4 + len;
  }
  uint64_t H = (at + 7) / 8 * 8, mm = (uint64_t) m * m;
  if (size != H + 16 * N + (uint64_t) n * N * mm * rs)
    goto done;

  /* the leaves: boxes and levels */
  raw = malloc (16 * N > N * mm * rs ? 16 * N : N * mm * rs);
  box = malloc (4 * N * sizeof (double));
  level = malloc (N * sizeof (float));
  v = malloc (N * mm * sizeof (float));
  if (!raw || !box || !level || !v || read_at (fp, size, H, raw, 16 * N))
    goto done;
  for (uint64_t q = 0; q < N; q++) {
    int32_t l = (int32_t) get32 (raw + 16 * q);
    int32_t i = (int32_t) get32 (raw + 16 * q + 4);
    int32_t j = (int32_t) get32 (raw + 16 * q + 8);
    int32_t k = (int32_t) get32 (raw + 16 * q + 12);
    if (l < 0 || l > 30 || i < 0 || j < 0 || i >= 1 << l || j >= 1 << l ||
        k != 0)
      goto done;
    double w = h0 / (double) (1L << l);
    box[4 * q] = o[0] + i * w, box[4 * q + 1] = o[1] + j * w;
    box[4 * q + 2] = box[4 * q + 3] = w;
    level[q] = (float) l;
  }

  /* the layer, then a field a scalar: the file's cells are leaf by leaf,
     i fastest, as WEBGL_CELL */
  const char * base = strrchr (path, '/');
  layer = webgl_blocks (S, base ? base + 1 : path, (int) N, (int) m, box);
  if (layer < 0)
    goto done;
  for (uint32_t k = 0; k < n; k++) {
    if (read_at (fp, size, H + 16 * N + (uint64_t) k * N * mm * rs, raw,
                 N * mm * rs)) {
      layer = -1;              /* the size was checked: not expected */
      goto done;
    }
    for (uint64_t c = 0; c < N * mm; c++)
      if (rs == 8) {
        double x;
        memcpy (&x, raw + 8 * c, 8);
        v[c] = (float) x;
      }
      else
        memcpy (&v[c], raw + 4 * c, 4);
    webgl_field (S, layer, names[k], WEBGL_CELL, v, 0);
  }
  webgl_field (S, layer, "level", WEBGL_BLOCK, level, 1);

done:
  for (uint32_t k = 0; names && k < n; k++)
    free (names[k]);
  free (names), free (box), free (level), free (v), free (raw);
  fclose (fp);
  return layer;
}
