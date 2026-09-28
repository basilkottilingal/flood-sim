#ifndef _ALLOC_GRID_H_
#define _ALLOC_GRID_H_

  #include <stdlib.h>

  /* API(s) for 2D grid */
  void ** grid_general (size_t type_size, int width, int height, int nghosts);
  #define grid(type,w,h)    (type **) grid_general (sizeof (type),w,h,2)
  #define grid_free(g)      free ( &g[-2][-2] ), free ( &g[-2] )
  
#endif
