#ifndef _ALLOC_GRID_H_
#define _ALLOC_GRID_H_

  #include <stdlib.h>

  typedef struct
  {
    int width, height, nghosts;
    size_t size;
  } GridData;

  /* API(s) for 2D grid */
  void ** grid_general (size_t type_size, int width, int height, int nghosts);
  #define grid(type,w,h)    (type **) grid_general (sizeof (type),w,h,2)
  #define grid_free(g)      free (((char *)&g[-2][-2]) - sizeof (GridData)), free (&g[-2])
  #define grid_data(g)      *(GridData *) (((char *)&g[-2][-2]) - sizeof (GridData))
  #define grid_copy(type,from,to,w,h) \
      memcpy (&to[-2][-2], &from[-2][-2], (w+4)*(h+4)*sizeof(type))
  
#endif
