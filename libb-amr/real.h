#ifndef _BAMR_REAL_H_
#define _BAMR_REAL_H_

  /*
  .. The precision of field values: double, or float with -DBAMR_SINGLE.
  ..
  ..   Real          a field value, as blocks store it
  ..   RealSum       an accumulator: always double, as a long sum in float
  ..                 loses digits
  ..   R (x)         a constant as a Real: R (0.5). A plain 0.5 is a double,
  ..                 and would turn float arithmetic into double and back
  ..   REAL_EPSILON  the epsilon of Real, for tolerances: a solver cannot
  ..                 reach 1e-12 in float
  ..
  .. Only field values are Real. Geometry stays double: at level 20 with
  .. 16 cells a block, a cell is 2^-24 wide, the whole mantissa of a float,
  .. so positions and cell sizes in float would stop being distinct.
  ..
  .. <tgmath.h> makes sqrt (), sin (), fabs (), ... take the type of their
  .. argument: sqrtf () for a float, sqrt () for a double.
  ..
  .. The library and the program must agree: each file that includes this
  .. header refers to bamr_real_is_float or bamr_real_is_double, and the
  .. library defines only the one it was built with, so a program built with
  .. the other precision does not link.
  */

  #include <float.h>
  #include <tgmath.h>

  #ifdef BAMR_SINGLE
    typedef float  Real;
    #define REAL_EPSILON  FLT_EPSILON
    #define REAL_MAX      FLT_MAX
    #define BAMR_REAL_TAG bamr_real_is_float
  #else
    typedef double Real;
    #define REAL_EPSILON  DBL_EPSILON
    #define REAL_MAX      DBL_MAX
    #define BAMR_REAL_TAG bamr_real_is_double
  #endif

  typedef double RealSum;

  #define R(x)  ((Real) (x))

  /* the precision check: see above */
  extern const int BAMR_REAL_TAG;
  static const int * const bamr_real_check __attribute__ ((used)) =
    &BAMR_REAL_TAG;

#endif
