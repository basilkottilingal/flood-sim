/* Scalars: ids from a fixed pool of NSCALAR_MAX.
 *
 * cwrap includes this header in every file it compiles (gcc -include) and
 * links scalar.o into every program, so that
 *
 *   scalar g;                  at file scope, gets an id before main () runs
 *   scalar x;                  in a function, gets scalar_new (), and
 *                              scalar_free (x) as its scope ends
 *   scalar k = scalar_new ();  in a function, keeps its id until freed
 */
#ifndef CWRAP_SCALAR_H
#define CWRAP_SCALAR_H

  /* 64 scalars: a practical limit, at least for now */
  #define NSCALAR_MAX 64

  typedef int scalar;

  /* The scalars at file scope, which take ids 0 to n_global_scalars - 1 */
  extern int n_global_scalars;

  /* api */
  scalar scalar_new    (void);         /* the lowest free id; aborts if none */
  void   scalar_free   (scalar);       /* an id not in use is left alone */
  int    scalar_in_use (scalar);
  int    scalar_count  (void);         /* ids in use */

  /* For cwrap: __attribute__ ((cleanup)) passes the scalar's address. */
  static inline void __cwrap_scalar_free (scalar * s)
  {
    scalar_free (*s);
  }

#endif
