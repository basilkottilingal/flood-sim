#ifndef _BITS_H_
#define _BITS_H_

  /*
  .. Bit tricks on unsigned integers of 8, 16, 32 and 64 bits.
  ..
  ..   bits_popcount_uN (x)   set bits
  ..   bits_ctz_uN (x)        trailing zeros: N if x is 0
  ..   bits_clz_uN (x)        leading zeros: N if x is 0
  ..   bits_lowest_uN (x)     index of the lowest set bit: -1 if x is 0
  ..   bits_highest_uN (x)    index of the highest set bit: -1 if x is 0
  ..   bits_log2_floor_uN (x) floor (log2 x): bits_highest; -1 if x is 0
  ..   bits_log2_ceil_uN (x)  ceil (log2 x): 0 for 1; -1 if x is 0
  ..   bits_is_pow2_uN (x)    1 if x is a power of 2 (not 0), else 0
  ..   bits_ceil_pow2_uN (x)  the smallest power of 2 >= x: 1 for 0, and 0
  ..                          if it does not fit in N bits
  ..   bits_pop_lowest_uN (&x) index of the lowest set bit, which is cleared
  ..                          from x: -1, and x left alone, if x is 0. For a
  ..                          loop over the set bits:
  ..                            for (int b; (b = bits_pop_lowest_u64 (&x)) >= 0; )
  ..   bits_parity_uN (x)     popcount mod 2
  ..
  .. N is 8, 16, 32 or 64; the uN types are uint8_t ... uint64_t, and the
  .. results int, except bits_ceil_pow2_uN, a uN. Every function is defined
  .. for every argument, 0 included: none has undefined behaviour.
  ..
  .. In C11 and later, the same names without the suffix choose by the type
  .. of x: bits_highest (x), bits_popcount (x), ... A uint8_t, uint16_t,
  .. uint32_t or uint64_t takes its own width; any other type is converted
  .. to uint64_t. An int constant is then 64 bits wide: bits_clz (1) is 63.
  .. bits_pop_lowest (&x) chooses by the type pointed to.
  ..
  .. How. With GCC or clang, the builtins (__builtin_ctz, _clz, _popcount
  .. and their l, ll forms), each the instruction where the target has one
  .. (tzcnt, lzcnt, popcnt with -mbmi, -mlzcnt, -mpopcnt or -march=native;
  .. bsf, bsr and a few instructions without). Otherwise, or with
  .. -DBITS_NO_BUILTINS:
  ..   ctz       x & -x isolates the lowest set bit, a power of 2; times a de
  ..             Bruijn constant, its top 5 (6) bits are unique for each power,
  ..             and a 32 (64) entry table gives the index;
  ..   clz       the bits below the highest set one are smeared to ones, then
  ..             x ^ (x >> 1) keeps only the highest: the same table;
  ..   popcount  the sum of bits in parallel (SWAR): pairs, nibbles, bytes,
  ..             then one multiply adds the bytes.
  .. 8 and 16 bits use the 32 bit forms.
  ..
  .. Builds with C99 or later (stdint.h), strict ISO or GNU, and with GNU
  .. C89 (-std=gnu89). bits.c holds the two de Bruijn tables, which the
  .. fallback reads: link it in either way.
  */

  #include <stdint.h>
  #include <limits.h>

  #if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define BITS_INLINE static inline
  #elif defined (__GNUC__)
    #define BITS_INLINE static __inline__
  #else
    #define BITS_INLINE static
  #endif

  /* the builtins: GCC 3.4 on, and clang, which defines __GNUC__ too */
  #if !defined (BITS_NO_BUILTINS) && defined (__GNUC__) && \
      (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
    #define BITS_HAVE_BUILTINS 1
  #else
    #define BITS_HAVE_BUILTINS 0
  #endif

  /*
  .. The de Bruijn sequences and their tables (bits.c): for a power of 2,
  .. p = 2^k, table [(p * B) >> (N - log2 N)] is k.
  */
  #define BITS_DEBRUIJN_32  0x077CB531u
  #define BITS_DEBRUIJN_64  UINT64_C (0x03F79D71B4CB0A89)
  extern const unsigned char bits_debruijn_table_32 [32];
  extern const unsigned char bits_debruijn_table_64 [64];

  /* ------------------------------------------------------------------ */
  /* 32 and 64 bits: the cores                                          */
  /* ------------------------------------------------------------------ */

  /* the index of p, a power of 2 */
  BITS_INLINE int bits_debruijn_index_32 (uint32_t p)
  {
    return bits_debruijn_table_32 [(uint32_t) (p * BITS_DEBRUIJN_32) >> 27];
  }

  BITS_INLINE int bits_debruijn_index_64 (uint64_t p)
  {
    return bits_debruijn_table_64 [(uint64_t) (p * BITS_DEBRUIJN_64) >> 58];
  }

  BITS_INLINE int bits_ctz_u32 (uint32_t x)
  {
    if (!x)
      return 32;
  #if BITS_HAVE_BUILTINS && UINT_MAX >= 0xFFFFFFFFu
    return __builtin_ctz ((unsigned int) x);
  #elif BITS_HAVE_BUILTINS
    return __builtin_ctzl ((unsigned long) x);
  #else
    return bits_debruijn_index_32 (x & (0u - x));
  #endif
  }

  BITS_INLINE int bits_ctz_u64 (uint64_t x)
  {
    if (!x)
      return 64;
  #if BITS_HAVE_BUILTINS && ULONG_MAX >= 0xFFFFFFFFFFFFFFFFu
    return __builtin_ctzl ((unsigned long) x);
  #elif BITS_HAVE_BUILTINS
    return __builtin_ctzll ((unsigned long long) x);
  #else
    return bits_debruijn_index_64 (x & (UINT64_C (0) - x));
  #endif
  }

  BITS_INLINE int bits_clz_u32 (uint32_t x)
  {
    if (!x)
      return 32;
  #if BITS_HAVE_BUILTINS && UINT_MAX == 0xFFFFFFFFu
    return __builtin_clz ((unsigned int) x);
  #elif BITS_HAVE_BUILTINS
    /* a wider unsigned long: its leading zeros, less the extra width */
    return __builtin_clzl ((unsigned long) x)
           - (int) (sizeof (unsigned long) * CHAR_BIT - 32);
  #else
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return 31 - bits_debruijn_index_32 (x ^ (x >> 1));
  #endif
  }

  BITS_INLINE int bits_clz_u64 (uint64_t x)
  {
    if (!x)
      return 64;
  #if BITS_HAVE_BUILTINS && ULONG_MAX == 0xFFFFFFFFFFFFFFFFu
    return __builtin_clzl ((unsigned long) x);
  #elif BITS_HAVE_BUILTINS
    return __builtin_clzll ((unsigned long long) x)
           - (int) (sizeof (unsigned long long) * CHAR_BIT - 64);
  #else
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return 63 - bits_debruijn_index_64 (x ^ (x >> 1));
  #endif
  }

  BITS_INLINE int bits_popcount_u32 (uint32_t x)
  {
  #if BITS_HAVE_BUILTINS && UINT_MAX >= 0xFFFFFFFFu
    return __builtin_popcount ((unsigned int) x);
  #elif BITS_HAVE_BUILTINS
    return __builtin_popcountl ((unsigned long) x);
  #else
    x = x - ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x + (x >> 4)) & 0x0F0F0F0Fu;
    return (int) ((uint32_t) (x * 0x01010101u) >> 24);
  #endif
  }

  BITS_INLINE int bits_popcount_u64 (uint64_t x)
  {
  #if BITS_HAVE_BUILTINS && ULONG_MAX >= 0xFFFFFFFFFFFFFFFFu
    return __builtin_popcountl ((unsigned long) x);
  #elif BITS_HAVE_BUILTINS
    return __builtin_popcountll ((unsigned long long) x);
  #else
    x = x - ((x >> 1) & UINT64_C (0x5555555555555555));
    x = (x & UINT64_C (0x3333333333333333))
        + ((x >> 2) & UINT64_C (0x3333333333333333));
    x = (x + (x >> 4)) & UINT64_C (0x0F0F0F0F0F0F0F0F);
    return (int) ((uint64_t) (x * UINT64_C (0x0101010101010101)) >> 56);
  #endif
  }

  /* ------------------------------------------------------------------ */
  /* 8 and 16 bits: from the 32 bit forms                               */
  /* ------------------------------------------------------------------ */

  BITS_INLINE int bits_ctz_u8 (uint8_t x)
  {
    return x ? bits_ctz_u32 (x) : 8;
  }

  BITS_INLINE int bits_ctz_u16 (uint16_t x)
  {
    return x ? bits_ctz_u32 (x) : 16;
  }

  BITS_INLINE int bits_clz_u8 (uint8_t x)
  {
    return bits_clz_u32 (x) - 24;
  }

  BITS_INLINE int bits_clz_u16 (uint16_t x)
  {
    return bits_clz_u32 (x) - 16;
  }

  BITS_INLINE int bits_popcount_u8 (uint8_t x)
  {
    return bits_popcount_u32 (x);
  }

  BITS_INLINE int bits_popcount_u16 (uint16_t x)
  {
    return bits_popcount_u32 (x);
  }

  /* ------------------------------------------------------------------ */
  /* the rest, the same for every width                                 */
  /* ------------------------------------------------------------------ */

  #define BITS_DEFINE_DERIVED(N)                                          \
                                                                          \
  BITS_INLINE int bits_lowest_u##N (uint##N##_t x)                        \
  {                                                                       \
    return x ? bits_ctz_u##N (x) : -1;                                    \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_highest_u##N (uint##N##_t x)                       \
  {                                                                       \
    return (N - 1) - bits_clz_u##N (x);           /* -1 for 0 */          \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_log2_floor_u##N (uint##N##_t x)                    \
  {                                                                       \
    return bits_highest_u##N (x);                                         \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_is_pow2_u##N (uint##N##_t x)                       \
  {                                                                       \
    return x != 0 && (x & (uint##N##_t) (x - 1u)) == 0;                   \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_log2_ceil_u##N (uint##N##_t x)                     \
  {                                                                       \
    if (!x)                                                               \
      return -1;                                                          \
    return bits_highest_u##N (x) + !bits_is_pow2_u##N (x);                \
  }                                                                       \
                                                                          \
  BITS_INLINE uint##N##_t bits_ceil_pow2_u##N (uint##N##_t x)             \
  {                                                                       \
    int k = bits_log2_ceil_u##N (x);              /* -1 for 0 */          \
    if (k <= 0)                                                           \
      return 1;                                                           \
    if (k >= N)                                                           \
      return 0;                                   /* does not fit */      \
    return (uint##N##_t) ((uint##N##_t) 1 << k);                          \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_pop_lowest_u##N (uint##N##_t * x)                  \
  {                                                                       \
    int b = bits_lowest_u##N (*x);                /* -1 for 0 */          \
    *x = (uint##N##_t) (*x & (uint##N##_t) (*x - 1u));   /* 0 stays 0 */  \
    return b;                                                             \
  }                                                                       \
                                                                          \
  BITS_INLINE int bits_parity_u##N (uint##N##_t x)                        \
  {                                                                       \
    return bits_popcount_u##N (x) & 1;                                    \
  }

  BITS_DEFINE_DERIVED (8)
  BITS_DEFINE_DERIVED (16)
  BITS_DEFINE_DERIVED (32)
  BITS_DEFINE_DERIVED (64)

  /* ------------------------------------------------------------------ */
  /* C11: by the type of the argument                                   */
  /* ------------------------------------------------------------------ */

  #if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    #define BITS_GENERIC(f, x)                                            \
      _Generic ((x),                                                      \
        uint8_t:  bits_##f##_u8,                                          \
        uint16_t: bits_##f##_u16,                                         \
        uint32_t: bits_##f##_u32,                                         \
        default:  bits_##f##_u64) (x)

    #define bits_popcount(x)    BITS_GENERIC (popcount, x)
    #define bits_ctz(x)         BITS_GENERIC (ctz, x)
    #define bits_clz(x)         BITS_GENERIC (clz, x)
    #define bits_lowest(x)      BITS_GENERIC (lowest, x)
    #define bits_highest(x)     BITS_GENERIC (highest, x)
    #define bits_log2_floor(x)  BITS_GENERIC (log2_floor, x)
    #define bits_log2_ceil(x)   BITS_GENERIC (log2_ceil, x)
    #define bits_is_pow2(x)     BITS_GENERIC (is_pow2, x)
    #define bits_ceil_pow2(x)   BITS_GENERIC (ceil_pow2, x)
    #define bits_parity(x)      BITS_GENERIC (parity, x)
    #define bits_pop_lowest(p)                                            \
      _Generic ((p),                                                      \
        uint8_t *:  bits_pop_lowest_u8,                                   \
        uint16_t *: bits_pop_lowest_u16,                                  \
        uint32_t *: bits_pop_lowest_u32,                                  \
        uint64_t *: bits_pop_lowest_u64) (p)
  #endif

#endif
