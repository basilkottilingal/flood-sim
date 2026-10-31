#include "tree.h"
#include "scalar.h"
#include "morton-curve.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>

static RAllocator * alloc_scalars = NULL;
static RAllocator * alloc_block = NULL;

#define w SCALAR_ARRAY_WIDTH
#define  nn    (((NSCALAR_MAX + 31) & ~31) / 32)
/* bits which are available for 'scalar'(s) */
static uint32_t bits [nn] = {0xFFFFFFFF};

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
  return (unsigned) __builtin_ctz(v);
# else
  /*
  .. Portable fallback in case __builtin_ctz function not available.
  .. de Bruijn sequence lookup.
  */
  static const uint8_t tbl[32] =
  {
    0,  1,  28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
    31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
  };
  return tbl [(uint32_t)((v & (0u - v)) * 0x077CB531u) >> 27];
# endif
}

static inline int popbit ()
{
  for (unsigned i = 0; i < nn; i++)
  {
    uint32_t v = bits [i];
    if (v != 0)
    {
      unsigned b = ctz32 (v);
      bits [i] = v & (v - 1);
      int at = (int) (i * 32 + b);
      return at >= NSCALAR_MAX ? -1 : at;
    }
  }
  return -1;
}

#define scalars() do {                        \
    for (unsigned i = 0; i < nn; i++)         \
    {                                         \
      uint32_t v = ~bits [i];                 \
      while (v) {                             \
        unsigned b = ctz32 (v);               \
        v &= (v - 1);                         \
        scalar s = (scalar) (i * 32 + b);     \
        (void) s;                             \
        /* do something here */
#define end()                                 \
      }                                       \
    }                                         \
  } while (0)                                 \

static inline void pushbit (int b)
{
  assert (b >=0 && b < NSCALAR_MAX);
  bits [b/32] |= (uint32_t) 1 << (b & 31);
}

scalar scalar_new ()
{
  int s = popbit (); 
  assert (s >= 0 && "hit nscalar limit");
  /* fixme : allocate */
  morton_start (0, 0, 0) {
    if (!sblocks [s/w])
      sblocks [s/w] = r_allocate (alloc_scalars);
    sblocks [s/w][s%w] = r_allocate (alloc_block);
  } morton_end ();
  return (scalar) s;
}

void scalar_free (scalar s)
{
  /* fixme : free allocated blocks */
  pushbit ((int) s);
}

vector vector_new ()
{
  vector v = (vector) {-1,-1,-1};
  scalar * V = &v.x;
  for (int d=0; d < dimension; ++d)
    V [d] = scalar_new ();
  return v;
}

void vector_free (vector v)
{
  scalar * V = &v.x;
  for (int d=0; d < dimension; ++d)
    scalar_free (V[d]);
}

void scalar_init (int blocksize)
{
  alloc_scalars = r_allocator (w * sizeof (void *));
  alloc_block   = r_allocator (blocksize * sizeof (double));
}

void scalar_free_all ()
{
  r_allocator_free (alloc_block);
  r_allocator_free (alloc_scalars);
  alloc_block  = NULL;
  alloc_scalars  = NULL;
}

void scalar_blocks (void **** sblocksp)
{
  assert (alloc_block && alloc_scalars);
  void *** sblocks = *sblocksp = r_allocate (alloc_scalars);
  scalars () {
    if (!sblocks [s/w])
      sblocks [s/w] = r_allocate (alloc_scalars);
    sblocks [s/w][s%w] = r_allocate (alloc_block);
  } end ();
}

void scalar_blocks_free (void *** sblocks)
{
  assert (sblocks && alloc_scalars && alloc_block);
  scalar max = 0;
  scalars () {
    assert (sblocks [s/w] && sblocks [s/w][s%w]);
    r_allocate_free (alloc_block, sblocks [s/w][s%w]);
    max = s; /* s is already sorted */
  } end ();

  int lim = ((int) max + w - 1) / w;
  for (int i=0; i<lim; ++i)
    if (sblocks [i/w])
      r_allocate_free (alloc_scalars, sblocks [i/w]);
  r_allocate_free (alloc_scalars, sblocks);
}

