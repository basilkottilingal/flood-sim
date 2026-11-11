/* Scalars: ids from a fixed pool of NSCALAR_MAX.
 *
 * With -dimension=1, 2 or 3, cwrap includes this header in every file it
 * compiles (gcc -include), with CWRAP_DIMENSION defined as that dimension. It
 * links scalar.o into every program. Then
 *
 *   scalar g;                  at file scope, gets an id before main () runs
 *   scalar x;                  in a function, gets scalar_new (), and
 *                              scalar_free (x) as its scope ends
 *   scalar k = scalar_new ();  in a function, keeps its id until freed
 *
 * In 2 and 3 dimensions, vector and tensor are structs of scalars, and each
 * of their scalars is numbered or freed in the same way, in the order
 * x, y (, z), and for a tensor x.x, x.y, ... y.x, ...
 */
#ifndef CWRAP_SCALAR_H
#define CWRAP_SCALAR_H

  /* 64 scalars: a practical limit, at least for now */
  #define NSCALAR_MAX 64

  typedef int scalar;

  #if CWRAP_DIMENSION == 2
  typedef struct { scalar x, y; } vector;
  typedef struct { vector x, y; } tensor;
  #elif CWRAP_DIMENSION == 3
  typedef struct { scalar x, y, z; } vector;
  typedef struct { vector x, y, z; } tensor;
  #endif

  /* The scalars at file scope, which take ids 0 to n_global_scalars - 1 */
  extern int n_global_scalars;

  /* api */
  scalar scalar_new    (void);         /* the lowest free id; aborts if none */
  void   scalar_free   (scalar);       /* an id not in use is left alone;
                                          a global one stops the run: a
                                          global scalar is never freed */
  int    scalar_in_use (scalar);
  int    scalar_count  (void);         /* ids in use */

  /*
  .. If not NULL, scalar_free (s) calls scalar_free_hook (s) first, while s is
  .. still in use, so the hook can release whatever it keeps for s. It is not
  .. called for an id not in use, so freeing twice calls it once.
  */
  extern void (* scalar_free_hook) (scalar);

  /*
  .. If not NULL, scalar_new () calls scalar_new_hook (s) with the id it
  .. hands out, once s is in use, so the hook can set up what it keeps for s.
  .. Scalars at file scope get their ids before main () runs, before any hook
  .. can be set; foreach_scalar (s) reaches them.
  */
  extern void (* scalar_new_hook) (scalar);

  /*
  .. The free ids, for foreach_scalar (): bit b of word w is set while the id
  .. 64 w + b is free. Bits from NSCALAR_MAX on stay clear and are never
  .. handed out. Only scalar.c writes it.
  */
  #define __CWRAP_SCALAR_WORDS ((NSCALAR_MAX + 63) / 64)
  extern unsigned long long __cwrap_scalar_free_bits [__CWRAP_SCALAR_WORDS];

  /* the lowest set bit of v, which is not 0 */
  static inline int __cwrap_scalar_ctz (unsigned long long v)
  {
  #if defined (__GNUC__) || defined (__clang__)
    return __builtin_ctzll (v);
  #else
    /* de Bruijn lookup, on the lowest 64 bits */
    static const unsigned char tbl [64] =
    {
       0,  1, 48,  2, 57, 49, 28,  3, 61, 58, 50, 42, 38, 29, 17,  4,
      62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12,  5,
      63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
      46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19,  9, 13,  8,  7,  6
    };
    unsigned long long low = (v & (0ull - v)) & 0xFFFFFFFFFFFFFFFFull;
    return tbl [((low * 0x03F79D71B4CB0A89ull) & 0xFFFFFFFFFFFFFFFFull) >> 58];
  #endif
  }

  /* the ids in use in word w, as bits */
  static inline unsigned long long __cwrap_scalar_used (int w)
  {
    unsigned long long all = (w == __CWRAP_SCALAR_WORDS - 1 && NSCALAR_MAX % 64)
      ? (1ull << (NSCALAR_MAX % 64)) - 1 : 0xFFFFFFFFFFFFFFFFull;
    return ~__cwrap_scalar_free_bits [w] & all;
  }

  /*
  .. foreach_scalar (s) { ... } runs the body once for each id in use, lowest
  .. first, with s that id, jumping from one id in use to the next (bit pops),
  .. so it costs one step per scalar in use and one per 64 ids.
  .. s is an int, not a scalar, so cwrap neither numbers it nor objects to it.
  .. break leaves the loop and continue goes to the next id, as in any loop.
  .. Each 64 ids are read when the loop reaches them: freeing s in the body is
  .. fine; an id taken in the body may or may not be visited.
  ..
  .. The innermost for runs once and declares s: on break, __cwrap_brk stays
  .. 1, which stops the two outer loops.
  */
  #define foreach_scalar(s)                                               \
    for (int __cwrap_w = 0, __cwrap_brk = 0;                              \
         !__cwrap_brk && __cwrap_w < __CWRAP_SCALAR_WORDS; __cwrap_w++)   \
      for (unsigned long long __cwrap_v = __cwrap_scalar_used (__cwrap_w); \
           !__cwrap_brk && __cwrap_v; __cwrap_v &= __cwrap_v - 1)          \
        for (int s = 64 * __cwrap_w + __cwrap_scalar_ctz (__cwrap_v),     \
                 __cwrap_k = (__cwrap_brk = 1);                           \
             ((void) (s), __cwrap_k); __cwrap_k = __cwrap_brk = 0)

  /* In 2 and 3 dimensions, also (below):
  ..   vector vector_new (void);    void vector_free (vector);
  ..   tensor tensor_new (void);    void tensor_free (tensor);
  */

  /* For cwrap: __attribute__ ((cleanup)) passes the scalar's address. */
  static inline void __cwrap_scalar_free (scalar * s)
  {
    scalar_free (*s);
  }

  /*
  .. Vectors and tensors: each scalar gets its id in the order x, y (, z), and
  .. for a tensor x.x, x.y, ... y.x, ... They are freed in the reverse order.
  .. The locals here have initializers, so cwrap leaves them alone.
  */
  #if CWRAP_DIMENSION == 2
  static inline vector vector_new (void)
  {
    scalar x = scalar_new (), y = scalar_new ();
    return (vector) { x, y };
  }

  static inline void vector_free (vector v)
  {
    scalar_free (v.y);
    scalar_free (v.x);
  }

  static inline tensor tensor_new (void)
  {
    vector x = vector_new (), y = vector_new ();
    return (tensor) { x, y };
  }

  static inline void tensor_free (tensor t)
  {
    vector_free (t.y);
    vector_free (t.x);
  }
  #elif CWRAP_DIMENSION == 3
  static inline vector vector_new (void)
  {
    scalar x = scalar_new (), y = scalar_new (), z = scalar_new ();
    return (vector) { x, y, z };
  }

  static inline void vector_free (vector v)
  {
    scalar_free (v.z);
    scalar_free (v.y);
    scalar_free (v.x);
  }

  static inline tensor tensor_new (void)
  {
    vector x = vector_new (), y = vector_new (), z = vector_new ();
    return (tensor) { x, y, z };
  }

  static inline void tensor_free (tensor t)
  {
    vector_free (t.z);
    vector_free (t.y);
    vector_free (t.x);
  }
  #endif

  /* For cwrap: the cleanups of local vectors and tensors. */
  #if CWRAP_DIMENSION >= 2
  static inline void __cwrap_vector_free (vector * v)
  {
    vector_free (*v);
  }

  static inline void __cwrap_tensor_free (tensor * t)
  {
    tensor_free (*t);
  }
  #endif

#endif
