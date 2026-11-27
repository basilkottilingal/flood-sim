#ifndef _WEBGL_TEST_DECODE_H_
#define _WEBGL_TEST_DECODE_H_

  /*
  .. For the tests: a written HTML file read back, a blob found by its id
  .. and decoded from base64 by a decoder of its own (not the writer's).
  */

  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>

  static char * read_file (const char * path, size_t * n)
  {
    FILE * fp = fopen (path, "rb");
    if (!fp)
      return NULL;
    fseek (fp, 0, SEEK_END);
    long len = ftell (fp);
    fseek (fp, 0, SEEK_SET);
    char * s = malloc ((size_t) len + 1);
    if (!s || fread (s, 1, (size_t) len, fp) != (size_t) len) {
      fclose (fp);
      free (s);
      return NULL;
    }
    fclose (fp);
    s[len] = '\0';
    if (n)
      *n = (size_t) len;
    return s;
  }

  static int b64_value (int c)
  {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
  }

  /* base64 text of length n to bytes, *len their number; NULL if the text
     is not base64 */
  static unsigned char * b64_decode (const char * t, size_t n, size_t * len)
  {
    if (n % 4)
      return NULL;
    unsigned char * out = malloc (n / 4 * 3 + 1);
    size_t k = 0;
    for (size_t i = 0; i < n; i += 4) {
      int v[4], pad = 0;
      for (size_t j = 0; j < 4; j++) {
        if (t[i + j] == '=' && i + 4 == n && j >= 2) {
          v[j] = 0, pad++;
          continue;
        }
        v[j] = b64_value ((unsigned char) t[i + j]);
        if (v[j] < 0 || pad) {
          free (out);
          return NULL;
        }
      }
      unsigned w = (unsigned) v[0] << 18 | (unsigned) v[1] << 12 |
        (unsigned) v[2] << 6 | (unsigned) v[3];
      out[k++] = (unsigned char) (w >> 16);
      if (pad < 2)
        out[k++] = (unsigned char) (w >> 8);
      if (pad < 1)
        out[k++] = (unsigned char) w;
    }
    *len = k;
    return out;
  }

  /* the text between <script id="ID" ...> and </script>, its length in *n */
  static const char * script_text (const char * html, const char * id,
                                   size_t * n)
  {
    char open[64];
    snprintf (open, sizeof (open), "<script id=\"%s\"", id);
    const char * a = strstr (html, open);
    if (!a || !(a = strchr (a, '>')))
      return NULL;
    a++;
    const char * b = strstr (a, "</script>");
    if (!b)
      return NULL;
    *n = (size_t) (b - a);
    return a;
  }

  /* blob ID decoded, its bytes in *len */
  static unsigned char * blob_bytes (const char * html, const char * id,
                                     size_t * len)
  {
    size_t n;
    const char * t = script_text (html, id, &n);
    return t ? b64_decode (t, n, len) : NULL;
  }

#endif
