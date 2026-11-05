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
  void   scalar_free   (scalar);       /* an id not in use is left alone */
  int    scalar_in_use (scalar);
  int    scalar_count  (void);         /* ids in use */

  /* For cwrap: __attribute__ ((cleanup)) passes the scalar's address. */
  static inline void __cwrap_scalar_free (scalar * s)
  {
    scalar_free (*s);
  }

  /* The same for vectors and tensors: the last scalar first. */
  #if CWRAP_DIMENSION >= 2
  static inline void __cwrap_vector_free (vector * v)
  {
  #if CWRAP_DIMENSION == 3
    scalar_free (v->z);
  #endif
    scalar_free (v->y);
    scalar_free (v->x);
  }

  static inline void __cwrap_tensor_free (tensor * t)
  {
  #if CWRAP_DIMENSION == 3
    __cwrap_vector_free (&t->z);
  #endif
    __cwrap_vector_free (&t->y);
    __cwrap_vector_free (&t->x);
  }
  #endif

#endif
