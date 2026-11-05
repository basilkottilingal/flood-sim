/* Caches: arrays that grow as needed, and rebuild themselves when marked dirty.
 *
 * Including this header turns on two pieces of syntax in cwrap, for the rest
 * of the file:
 *
 *   cache (T) v;         an array of T: v.a [0] to v.a [v.n - 1]
 *   foreach (x in v) ... each element in turn
 *
 * cwrap defines struct cache_T, once per file, as
 *
 *   struct cache_T {
 *     T *a;
 *     int n, cap, is_dirty;
 *     void (*update) (struct cache_T *self);
 *   };
 *
 * Every field starts at zero: an empty cache with no update. Caches at file
 * scope and static ones are zero already, and cwrap gives a local cache
 * `= { 0 }` and frees it as its scope ends. A cache anywhere else, in a struct
 * or from malloc (), needs cache_init () or calloc ().
 *
 * Before a foreach over a cache, cwrap adds
 *
 *   if (v.update && v.is_dirty) { v.is_dirty = 0; v.update (&v); }
 *
 * Setting update and is_dirty is up to you.
 *
 * Two caches must not share their array, so cwrap stops with an error at
 * `w = v;`, `cache (T) w = v;` and `return v;`, and warns when a cache is passed
 * to a function by value: pass &v instead. A local cache may still have a
 * brace initializer, as in `cache (int) v = { .update = rebuild };`. */
#ifndef CWRAP_CACHE_H
#define CWRAP_CACHE_H

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

  /* The first capacity: it doubles from there. */
  #define CACHE_MIN 8

  /* api: each takes the cache itself, not its address */
  #define cache_init(V)      ((void) memset (&(V), 0, sizeof (V)))
  #define cache_clear(V)     ((void) ((V).n = 0))  /* keeps the memory */
  #define cache_reserve(V, m)                                                \
    ({                                                                       \
      __typeof__ (V) * __cwrap_v = &(V);                                     \
      long __cwrap_m = (m);                                                  \
      if (__cwrap_m > __cwrap_v->cap)                                        \
        __cwrap_v->a = __cwrap_cache_grow (__cwrap_v->a, &__cwrap_v->cap,    \
                                           __cwrap_m, sizeof *__cwrap_v->a); \
      (void) 0;                                                              \
    })
  /* The value is read after the array grows, so it may be an element of V.
  .. It is the rest of the arguments, so that a compound literal with commas,
  .. as in cache_push (c, (point) { 1, 2 }), is still one value. */
  #define cache_push(V, ...)                                                 \
    ({                                                                       \
      __typeof__ (V) * __cwrap_v = &(V);                                     \
      if (__cwrap_v->n == __cwrap_v->cap)                                    \
        __cwrap_v->a = __cwrap_cache_grow (__cwrap_v->a, &__cwrap_v->cap,    \
                                           __cwrap_v->n + 1L,                \
                                           sizeof *__cwrap_v->a);            \
      __cwrap_v->a [__cwrap_v->n++] = (__VA_ARGS__);                         \
      (void) 0;                                                              \
    })
  #define cache_free(V)                                                      \
    ({                                                                       \
      __typeof__ (V) * __cwrap_v = &(V);                                     \
      free (__cwrap_v->a);                                                   \
      __cwrap_v->a = NULL;                                                   \
      __cwrap_v->n = __cwrap_v->cap = 0;                                     \
      (void) 0;                                                              \
    })

  /* For the macros: a with room for at least need elements of size bytes. */
  static inline void * __cwrap_cache_grow (void * a, int * cap, long need,
                                           size_t size)
  {
    long n = *cap > 0 ? *cap : CACHE_MIN;
    while (n < need && n <= INT_MAX / 2)
      n *= 2;
    if (n < need || (size_t) n > SIZE_MAX / size)
    {
      fprintf (stderr, "cache: %ld elements are too many\n", need);
      abort ();
    }
    a = realloc (a, (size_t) n * size);
    if (!a)
    {
      fprintf (stderr, "cache: out of memory for %ld elements\n", n);
      abort ();
    }
    *cap = (int) n;
    return a;
  }

  /* For cwrap: the cleanup of a local cache. a, the first field, is a pointer
   * in every struct cache_T. */
  static inline void __cwrap_cache_free (void * c)
  {
    void * a;
    memcpy (&a, c, sizeof a);
    free (a);
  }

#endif
