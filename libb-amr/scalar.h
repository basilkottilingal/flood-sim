#ifndef _BAMR_SCALAR_H
#define _BAMR_SCALAR_H_


  #ifndef NSCALAR_MAX
  # define NSCALAR_MAX 128
  #endif

  typedef int scalar;
  typedef struct {scalar x, y, z;} vector;

  /* api */
  scalar scalar_new  ();
  void   scalar_free (scalar);
  vector vector_new  ();
  void   vector_free (vector);

#endif
