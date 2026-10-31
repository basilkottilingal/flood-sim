#ifndef _BAMR_SCALAR_H
#define _BAMR_SCALAR_H_

  #include "liballoc/reuse-allocator.h"

  /* 64 scalars : a practical limit (atleast for now) */
  # define SCALAR_ARRAY_WIDTH   8
  # define NSCALAR_MAX         (SCALAR_ARRAY_WIDTH * SCALAR_ARRAY_WIDTH)

  typedef int scalar;
  typedef struct {scalar x, y, z;} vector;

  /* api */
  scalar scalar_new  ();
  void   scalar_free (scalar);
  vector vector_new  ();
  void   vector_free (vector);
  void   scalar_init ();
  void   scalar_free_all ();
  void   scalar_blocks (void **** sbp);
  void   scalar_blocks_free (void *** sb);

#endif
