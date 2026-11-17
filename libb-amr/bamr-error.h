#ifndef _BAMR_ERROR_H_
#define _BAMR_ERROR_H_

  #include <stdio.h>
  #include <stdlib.h>

  /*
  .. BAMR_ABORT (what): the run stopped, after one line on stderr: the
  .. file, the line and what went wrong. For what must not happen, as an
  .. allocation that fails, where assert () would vanish with NDEBUG.
  */
  #define BAMR_ABORT(what)  bamr_abort (__FILE__, __LINE__, (what))

  static inline _Noreturn void bamr_abort (const char * file, int line,
                                           const char * what)
  {
    fflush (stdout);
    fprintf (stderr, "%s:%d: %s\n", file, line, what);
    abort ();
  }

#endif
