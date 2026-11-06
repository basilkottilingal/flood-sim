#ifndef _BAMR_POOL_H_
#define _BAMR_POOL_H_

  #include <stddef.h>

  /*
  .. A pool of objects of one size, such as the blocks of a scalar.
  ..
  .. Objects come from chunks mapped with mmap: memory is committed only as it
  .. is touched, and pool_destroy () gives every chunk back at once. A freed
  .. object goes on a free list and is the next one handed out.
  ..
  .. pool_alloc () always returns zeroed memory. Objects of POOL_LARGE bytes or
  .. more are rounded up to whole pages; pool_free () hands their pages back to
  .. the system (madvise), so a coarsened grid does not keep the memory of its
  .. old blocks, and a reused object is zero without a memset.
  ..
  .. Smaller objects are aligned to 16 bytes, or 64 (a cache line) from 64
  .. bytes up; large ones to a page.
  */

  #ifndef POOL_LARGE
  #define POOL_LARGE  (16 * 1024)
  #endif

  /* the bytes of a chunk, at least: one object if it is larger */
  #ifndef POOL_CHUNK
  #define POOL_CHUNK  (1 << 20)
  #endif

  typedef struct _Pool Pool;

  /* api */
  Pool * pool_new     (size_t size);       /* objects of 'size' bytes */
  void * pool_alloc   (Pool * p);          /* zeroed; aborts if out of memory */
  void   pool_free    (Pool * p, void * o);   /* NULL is left alone */
  void   pool_destroy (Pool * p);          /* frees every object; NULL is fine */
  size_t pool_used    (const Pool * p);    /* objects handed out, not freed */
  size_t pool_size    (const Pool * p);    /* bytes per object, as rounded */

#endif
