/*
.. LZW enoding is a lossless encoding.
.. Here only lzw decoding routine is defined that too with the tiff specs.
.. You may refer to the libtiff repo for industry grade encoding/decoding here.
.. https://gitlab.com/libtiff/libtiff/-/blob/master/libtiff/tif_lzw.c
*/

#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>

#define CLEAR_CODE 256
#define EOI_CODE   257
#define FIRST_CODE 258
#define MAX_CODE   4095

typedef struct
{
  unsigned char * start;
  size_t seek;
  size_t size;
} Bitfile ;

static void bitfile_set (Bitfile * b, void * start, size_t len)
{
  b->start = (unsigned char *) start;
  b->seek = 0;
  b->size = len << 3;
}

/*
.. NOTE : TIFF use MSB ordering. That is Most significant bit is
.. read first
*/

static inline
unsigned int get_code (Bitfile * b, size_t len)
{
  /* assert (len >= 9 && len <= 12); */
  if (b->size - b->seek < len)
    return (unsigned int) EOF;
 
  unsigned char * stream = b->start; 
  unsigned int code = 0;
  for (size_t i = b->seek; i < b->seek + len; ++i) {
    code <<= 1;
    if (stream [i >> 3] & (1u << (7 - (i & 7))))
      code |= 1;
  }
  b->seek += len;
  return code; 
}

#define LZWD_ERR_NO_EOI             -1
#define LZWD_ERR_WRONG_ENCODING     -2
#define LZWD_ERR_INCOMPLETE_READ    -3
#define LZWD_ERR_BUFFER_OOM         -4
#define LZWD_ERR_OUTPUT_MISMATCH    -5

#define LZWD_ERROR(e)                                            \
  (e == 0) ? "success ":                                         \
  (e == LZWD_ERR_NO_EOI) ? "No EOI found" :                      \
  (e == LZWD_ERR_INCOMPLETE_READ) ? "Encode size mismatch" :     \
  (e == LZWD_ERR_BUFFER_OOM) ? "Buffer Out Of Memory" :          \
  (e == LZWD_ERR_OUTPUT_MISMATCH) ? "Unexpected Output length" : \
  (e == LZWD_ERR_WRONG_ENCODING) ? "Wrong lzw tiff encoding" :   \
  "unknown error"

void lzw_error (int e)
{
  fprintf (stderr, "lzw_decode () : %s\n", LZWD_ERROR (e));
}

/*
.. lzw_decode () : LZW decompression tailor made for TIFF
.. Args :
.. (a) encoding : Encoded data
.. (b) len      : size of "encoding". Expects last code in the encoding as EOI
.. (c) buff     : buffer to store decoded bytes
.. (d) outlen   : expected decode output length. "buff" size should be >= outlen.
.. Return :
..  error of type int. Please see LZWD_ERROR()
..
.. Assumes
.. 1. endcoding [0:len) is readable.
.. 2. buff [] has capacity of atleast "outlen" bytes.
*/
int lzw_decode (void * encoding, size_t len, void * buff, size_t outlen)
{

  #define put(c) do                       \
  {                                       \
    if (outlen--)                         \
      *out++ = c;                         \
    else                                  \
      return LZWD_ERR_BUFFER_OOM;         \
  } while (0)

  static unsigned char stack [MAX_CODE + 1];
  static struct {
    unsigned char suffix_char [MAX_CODE - FIRST_CODE];
    unsigned int prefix_code [MAX_CODE - FIRST_CODE]; 
  } dictionary;

  unsigned char * out = (unsigned char *) buff; 

  unsigned int code, last, next = FIRST_CODE;
  int code_len = 9;

  Bitfile bitfile;
  bitfile_set (& bitfile, encoding, len);

  while ( (last = get_code (& bitfile, code_len)) == CLEAR_CODE ) {}
  if (last >= 256)
    return LZWD_ERR_WRONG_ENCODING;

  unsigned char firstchar = (unsigned char) last;
  put (firstchar);

  while (1)
  {
    code = get_code (&bitfile, code_len);
    
    /* EOF */
    if ((int) code ==  EOF)
      return LZWD_ERR_NO_EOI;

    /* Error */
    if (code > next)
      return LZWD_ERR_WRONG_ENCODING;

    /* End of encoded stream */
    if (code == EOI_CODE)
      return get_code (&bitfile, code_len) == EOF ?
        (outlen == 0 ? 0 : LZWD_ERR_OUTPUT_MISMATCH) : 
        LZWD_ERR_INCOMPLETE_READ;

    /* reset */
    if (code == CLEAR_CODE) {
      next = FIRST_CODE;
      code_len = 9;
      while ( (last = get_code (& bitfile, code_len)) == CLEAR_CODE ) {}
      if (last > CLEAR_CODE)
        switch (last) {
          case EOI_CODE :
            return get_code (&bitfile, code_len) == EOF ?
              (outlen == 0 ? 0 : LZWD_ERR_OUTPUT_MISMATCH) : 
              LZWD_ERR_INCOMPLETE_READ;
          case EOF :
            return LZWD_ERR_NO_EOI;
          default : 
            return LZWD_ERR_WRONG_ENCODING;
        }
      firstchar = (unsigned char) last;
      put (firstchar);
      continue;
    }

    /* using stack to avoid recursive function */  
    int depth = 0;                     
    unsigned int icode = code == next ? last : code; 
    if (code == next)
      stack [depth ++] = firstchar;
    while (icode >= FIRST_CODE) {
      stack [depth ++] = dictionary.suffix_char [icode - FIRST_CODE];
      icode = dictionary.prefix_code [icode - FIRST_CODE];
    }
    firstchar = stack [depth ++] = (unsigned char) icode;
    
    while (depth--) {
      put (stack [depth]);
    }

    if (next <= MAX_CODE)
    {
      dictionary.prefix_code[next - FIRST_CODE] = last;
      dictionary.suffix_char[next - FIRST_CODE] = firstchar;
      next++;

      if ((1u<<code_len) - 1 == next && code_len < 12)
        code_len++;
    }

    last = code;
  }

  #undef put

  return 0;
}

