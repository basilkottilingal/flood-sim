#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "scalar.h"

#define nn ((NSCALAR_MAX + 31) / 32)

/*
.. bit s of 'bits' is set while the id s is free. The bits from NSCALAR_MAX to
.. the end of the last word stay clear: those ids are never handed out.
*/
static uint32_t bits [nn];

# if (defined(__GNUC__) || defined(__clang__)) && (UINT_MAX == 0xFFFFFFFF)
#   if defined(__has_builtin)
#     if __has_builtin(__builtin_ctz)
#       define HAVE_CTZ
#     endif
#   endif
# endif

static inline unsigned ctz32 (uint32_t v)
{
# ifdef HAVE_CTZ
  return (unsigned) __builtin_ctz (v);
# else
  /*
  .. Portable fallback in case __builtin_ctz function not available.
  .. de Bruijn sequence lookup.
  */
  static const uint8_t tbl [32] =
  {
    0,  1,  28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
    31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
  };
  return tbl [(uint32_t) ((v & (0u - v)) * 0x077CB531u) >> 27];
# endif
}

/* the lowest free id, now taken, or -1 */
static inline int popbit ()
{
  for (unsigned i = 0; i < nn; i++)
  {
    uint32_t v = bits [i];
    if (v != 0)
    {
      bits [i] = v & (v - 1);
      return (int) (i * 32 + ctz32 (v));
    }
  }
  return -1;
}

static inline void pushbit (int b)
{
  bits [b/32] |= (uint32_t) 1 << (b & 31);
}

scalar scalar_new ()
{
  int s = popbit ();
  if (s < 0)
  {
    fprintf (stderr, "scalar_new: all %d scalars are in use (NSCALAR_MAX)\n",
             NSCALAR_MAX);
    abort ();
  }
  return (scalar) s;
}

int scalar_in_use (scalar s)
{
  return s >= 0 && s < NSCALAR_MAX
         && !((bits [s/32] >> (s & 31)) & 1);
}

void scalar_free (scalar s)
{
  /*
  .. an id not in use is left alone: freed twice, or never set, as when a goto
  .. jumps past the declaration of a local scalar.
  */
  if (scalar_in_use (s))
    pushbit ((int) s);
}

int scalar_count ()
{
  int nfree = 0;
  for (unsigned i = 0; i < nn; i++)
    for (uint32_t v = bits [i]; v; v &= v - 1)
      nfree++;
  return NSCALAR_MAX - nfree;
}

/*
.. cwrap puts a slot for each scalar at file scope (and each static scalar in
.. a function) in the section "scalars". The linker joins the slots of every
.. file between __start_scalars and __stop_scalars, which are weak for a
.. program with none. Priority 101 runs this before other constructors, so
.. every id is free before any code can ask for one.
*/
extern scalar * const __start_scalars [] __attribute__ ((weak));
extern scalar * const __stop_scalars [] __attribute__ ((weak));
int n_global_scalars = 0;

__attribute__ ((constructor (101))) static void number_global_scalars ()
{
  for (int s = 0; s < NSCALAR_MAX; s++)       /* all ids free */
    pushbit (s);
  for (scalar * const * p = __start_scalars; p < __stop_scalars; p++)
    **p = scalar_new ();
  n_global_scalars = (int) (__stop_scalars - __start_scalars);
}
