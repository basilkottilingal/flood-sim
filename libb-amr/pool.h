#ifndef _BAMR_POOL_H_
#define _BAMR_POOL_H_

  #include <stddef.h>

  /*
  .. A pool of objects of one size, such as the blocks of a scalar.
  ..
  .. Objects come from chunks mapped with mmap: memory is committed only as it
  .. is touched, and pool_destroy () gives every chunk back at once.
  ..
  .. Freed objects wait on one of two lists, and are handed out again, the
  .. most recently freed first, as that one is likely still in cache:
  ..
  ..   hot   still in memory, with whatever they held: reusing one costs no
  ..         page fault. Every small object goes here.
  ..   cold  large objects (POOL_LARGE bytes or more, whole pages) whose pages
  ..         went back to the system (madvise): they use no memory, and read
  ..         as zero, but touching them again costs a page fault per page.
  ..
  .. A large object freed goes on the hot list while it holds fewer than
  .. pool_set_hot_max () objects (no limit at first, as malloc keeps freed
  .. memory), and on the cold list otherwise. pool_trim () moves hot objects
  .. to the cold list, for after a coarsening that freed many blocks.
  ..
  .. Which list an allocation takes from:
  ..
  ..   pool_alloc ()         zeroed memory: a cold object (zero already), or
  ..                         a hot one cleared with memset, or a fresh one
  ..                         (zero already)
  ..   pool_alloc_uninit ()  any content, for an object written in full
  ..                         straight away: a hot one (no fault, no clearing),
  ..                         or a cold or fresh one
  ..
  .. So memory is zeroed at most once: by the kernel for fresh or cold pages,
  .. or by memset for a hot object, and not at all by pool_alloc_uninit ().
  ..
  .. Alignment: 16 bytes, 64 (a cache line) from 64 bytes up, a page for large
  .. objects. Not thread safe: one pool, one thread at a time.
  ..
  ..
  .. Releasing memory: _BAMR_MADV_DONTNEED_RECOMMENDED_
  .. ---------------------------------------------------
  ..
  .. 1 (MADV_DONTNEED): as above. A large object over the hot limit goes
  ..   cold at once: its pages go back to the system, and it reads as zero.
  ..   Its first page faults back in at once, as the free list's link is
  ..   written into it.
  ..
  .. 0 (MADV_FREE, Linux 4.5 on): a third list between the two,
  ..
  ..   muzzy  large objects over the hot limit, given back lazily: all but
  ..          the first page (which holds the list's link) are advised
  ..          MADV_FREE. The system takes those pages only if it runs short
  ..          of memory; until then a write keeps a page, with no fault.
  ..          Each page reads as what it held, or as zero if the system took
  ..          it, and may turn to zero until it is written.
  ..
  ..   pool_alloc ()         a cold object, or a hot or muzzy one cleared with
  ..                         memset, or a fresh one: never a muzzy one as it is
  ..   pool_alloc_uninit ()  a hot one, or a muzzy one (no fault unless the
  ..                         system took its pages), or a cold or fresh one
  ..   pool_set_hot_max ()   hot objects over the limit go muzzy
  ..   pool_trim ()          hot objects beyond keep, and every muzzy one, go
  ..                         cold: memory really given back, and out of the
  ..                         resident size (RSS), which counts muzzy pages
  ..                         until the system takes them
  ..
  ..   For blocks that are written in full before being read, as a refinement
  ..   writes its children: reuse costs no fault while memory is plentiful,
  ..   and the memory still goes back when it is not. Without MADV_FREE in
  ..   <sys/mman.h>, or if the kernel refuses it, objects go cold instead.
  ..
  .. The default is 1, the behaviour before the choice existed; build with
  .. -D_BAMR_MADV_DONTNEED_RECOMMENDED_=0 to try MADV_FREE.
  */

  #ifndef _BAMR_MADV_DONTNEED_RECOMMENDED_
  #define _BAMR_MADV_DONTNEED_RECOMMENDED_  1
  #endif

  #ifndef POOL_LARGE
  #define POOL_LARGE  (16 * 1024)
  #endif

  /* the bytes of a chunk, at least: one object if it is larger */
  #ifndef POOL_CHUNK
  #define POOL_CHUNK  (1 << 20)
  #endif

  /* no limit on the hot list */
  #define POOL_HOT_ALL  ((size_t) -1)

  typedef struct _Pool Pool;

  /* api */
  Pool * pool_new          (size_t size);  /* objects of 'size' bytes */
  void * pool_alloc        (Pool * p);     /* zeroed; aborts if out of memory */
  void * pool_alloc_uninit (Pool * p);     /* any content; aborts if out of memory */
  void   pool_free         (Pool * p, void * o);   /* NULL is left alone */
  void   pool_destroy      (Pool * p);     /* frees every object; NULL is fine */

  void   pool_set_hot_max  (Pool * p, size_t n);  /* large objects kept hot */
  void   pool_trim         (Pool * p, size_t keep);  /* hot beyond keep: cold */

  size_t pool_used         (const Pool * p);  /* handed out, not freed     */
  size_t pool_hot          (const Pool * p);  /* freed, still in memory    */
  size_t pool_cold         (const Pool * p);  /* freed, pages given back   */
  size_t pool_muzzy        (const Pool * p);  /* freed, pages given back
                                                 lazily: 0 with
                                                 MADV_DONTNEED            */
  size_t pool_size         (const Pool * p);  /* bytes per object, rounded */

#endif
