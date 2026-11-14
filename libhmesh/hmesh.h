#ifndef _HMESH_H_
#define _HMESH_H_

  /*
  .. A simplicial manifold of dimension m immersed in R^D, D the euclidean
  .. dimension (cwrap's -dimension=D), 0 <= m <= D <= 3: points (m = 0), a
  .. curve (1), a surface (2) or a volume (3), as a half-edge mesh generalised
  .. to every level of the complex.
  ..
  ..
  .. Levels and simplices
  .. --------------------
  ..
  .. Level k, 0 <= k <= m, holds the k-simplices: points, edges, triangles,
  .. tetrahedra. Every face of every simplex is stored, once, at its own
  .. level: a triangle's three edges are level 1 simplices, each with its
  .. own id, so that each level carries its own fields (scalars, vectors),
  .. values on points, on edges (fluxes), on triangles, on tetrahedra.
  ..
  .. A k-simplex s, k >= 1, is its k + 1 vertices v0 ... vk, in an order that
  .. is its orientation, and has k + 1 halves, one for each facet:
  ..
  ..   half h = (k + 1) s + i     the facet opposite vi, the (k-1)-simplex
  ..                              v0 ... vi ... vk without vi, as a face of s
  ..
  .. Each half has
  ..
  ..   hm_half_vertex ()  vi, the vertex opposite its facet
  ..   hm_half_facet ()   the facet's id, at level k - 1
  ..   hm_half_sign ()    +1 if the orientation s induces on the facet,
  ..                      (-1)^i (v0 ... vi ... vk without vi), is the
  ..                      facet's own, -1 otherwise: the incidence numbers of
  ..                      the boundary operator, d (v0, v1) = v1 - v0
  ..   hm_half_next ()    the next half of the same simplex, i + 1 mod k + 1;
  ..                      computed, not stored
  ..   hm_half_twin ()    the next half, of another k-simplex, on the same
  ..                      facet: the halves on a facet form a cycle
  ..
  .. so every simplex but a point has a twin and a next:
  ..
  ..   k  halves                           twin cycles about
  ..   1  2 half-points (edge ends)        a point: the edges at it
  ..   2  3 half-edges                     an edge: the triangles on it
  ..   3  4 half-faces                     a triangle: the tetrahedra on it
  ..
  .. At level m, inside a manifold, a facet has two halves, and twin is the
  .. half-edge structure's twin: twin (twin (h)) = h, the neighbour across
  .. the facet, with the opposite sign. On a boundary the cycle is h alone
  .. (twin (h) = h, hm_half_boundary ()); where more than two meet (a
  .. junction of three river branches, a fin) it holds them all. Below
  .. level m a cycle holds every coface: the edges at a point of a
  .. triangulated surface, in no particular order.
  ..
  .. For a triangle the halves are the half-edges of a half-edge mesh: half
  .. i goes from v(i+1) to v(i+2) (indices mod 3), hm_half_origin ().
  ..
  .. A (k-1)-simplex, k <= m, keeps one half of level k on it, hm_use ()
  .. (HM_NONE if no k-simplex has it as a facet), the entry to its cycle:
  .. hm_use () of a point, at level 1, starts the walk through its edges.
  ..
  ..
  .. Storage
  .. -------
  ..
  .. A structure of arrays per level: a k-simplex's vertices, twins and
  .. facets at [(k + 1) s + i], each field component an array of its own
  .. (the coordinates of the points too: x[0], x[1], x[2]). Ids are 32 bit
  .. (up to 2^31 - 1 simplices per level).
  ..
  .. A simplex is added and removed in O(1) wherever it is: its id comes from
  .. a stack of freed ids (the last freed first, as its memory is still in
  .. cache), or the end of the arrays, which grow by doubling. An id stays
  .. the same until hm_compact (), so that twins, facets, uses and ids held
  .. outside stay valid through any sequence of additions and removals (a
  .. flip, a split, a collapse). A bit per id says whether it is in use;
  .. hm_foreach () visits those, skipping the holes 64 at a time.
  ..
  .. hm_compact () closes the holes: the simplices of each level move down,
  .. in the order of their ids, to 0 ... n - 1, and every reference follows;
  .. it can give back each level's map from old ids to new, for the links
  .. held outside the mesh (to blocks of the tree, to ranks).
  ..
  ..
  .. Adding and removing
  .. -------------------
  ..
  ..   hm_vertex_add (M, x)          a point at x (D coordinates)
  ..   hm_add (M, k, v)              the k-simplex of vertices v [0 ... k], in
  ..                                 that orientation, its faces too, found or
  ..                                 made, each new face in the orientation
  ..                                 its first coface induces; the one
  ..                                 already there if any (whatever the order
  ..                                 of v, its own orientation kept)
  ..   hm_find (M, k, v)             the k-simplex of vertices v, any order,
  ..                                 or HM_NONE: through the cycle of the
  ..                                 facet without v[0], no hash table
  ..   hm_remove (M, k, s)           s alone: refused (-1) while a (k+1)-
  ..                                 simplex has it as a face
  ..   hm_remove_closure (M, k, s)   s, and every face of it no other simplex
  ..                                 has left, points included
  ..
  .. Fields: hm_field_new (M, k, name, ncomp, esize) gives an id, for ncomp
  .. arrays of elements of esize bytes, one element per k-simplex,
  .. hm_field (M, f, c) the array of component c. hm_scalar_new () and
  .. hm_vector_new () make HmReal ones (double, or float with -DBAMR_SINGLE
  .. as libb-amr's Real), of 1 and D components. A new simplex's values are
  .. zero. A field's arrays move as its level grows: take the pointer again
  .. after an addition.
  ..
  .. One manifold for the program: manifold (dim) makes it as hmesh, and
  .. manifold_free () frees it; the functions take the mesh, for tests and
  .. for later couplings of several.
  ..
  .. Not thread safe.
  */

  #include <stdint.h>
  #include <stddef.h>
  #include <string.h>

  #include "libutils/bits.h"

  /* the euclidean dimension: cwrap's, or -DHM_DIMENSION=D without cwrap */
  #ifdef CWRAP_DIMENSION
    #define HM_DIMENSION CWRAP_DIMENSION
  #elif !defined (HM_DIMENSION)
    #error "build with cwrap -dimension=D, or -DHM_DIMENSION=D"
  #endif
  #if HM_DIMENSION < 1 || HM_DIMENSION > 3
    #error "HM_DIMENSION must be 1, 2 or 3"
  #endif

  #ifdef BAMR_SINGLE
    typedef float  HmReal;
  #else
    typedef double HmReal;
  #endif

  #define HM_LEVELS  4               /* levels 0 ... 3 */
  #define HM_NONE    UINT32_MAX      /* no simplex, no half */
  #define HM_FLIP    0x80000000u     /* a facet's sign bit, in Level.facet */
  #define HM_ID_MAX  0x7fffffffu     /* ids below */

  typedef struct {
    uint32_t   end;       /* ids < end have been used */
    uint32_t   cap;       /* ids < cap have memory */
    uint32_t   n;         /* ids in use */
    uint64_t * alive;     /* a bit per id: in use */
    uint32_t * freed;     /* ids < end not in use, the last freed on top */
    uint32_t   nfreed, capfreed;
    /* k >= 1: (k + 1) per simplex */
    uint32_t * vert;      /* the vertex opposite each half's facet */
    uint32_t * twin;      /* the next half on the same facet */
    uint32_t * facet;     /* the facet's id, HM_FLIP if its sign is -1 */
    /* k < m: one per simplex */
    uint32_t * use;       /* a half of level k + 1 on it, or HM_NONE */
  } HmLevel;

  typedef struct {
    int      level;       /* -1: a free field id */
    int      ncomp;
    int      real;        /* HmReal values (hm_scalar_new, hm_vector_new):
                             interpolated where hm_clip () cuts an edge */
    size_t   esize;
    char   * name;
    void  ** data;        /* ncomp arrays of level's cap elements */
  } HmField;

  typedef struct {
    int       dim;                     /* m */
    HmLevel   level[HM_LEVELS];
    double  * x[HM_DIMENSION];         /* the points' coordinates */
    HmField * field;
    int       nfield, capfield;
    /* the box the points must stay in: hm_set_domain () */
    int       domain;
    double    origin[HM_DIMENSION], size;
    /* changes with the points: added, removed, moved, renumbered, or the
       box changed; a cache of the points compares it */
    unsigned  generation;
    /* a cache another layer keeps with the mesh (hmesh-tree.h), freed with
       it by ext_free () */
    void    * ext;
    void   (* ext_free) (void *);
  } Hmesh;

  /* the library and the program must agree on D: see real.h in libb-amr */
  #define HM_TAG_(d)   hmesh_dimension_is_##d
  #define HM_TAG(d)    HM_TAG_(d)
  extern const int HM_TAG (HM_DIMENSION);
  static const int * const hmesh_dimension_check __attribute__ ((used)) =
    &HM_TAG (HM_DIMENSION);

  /* the mesh */
  Hmesh * hmesh_new (int dim);       /* NULL unless 0 <= dim <= D */
  void    hmesh_free (Hmesh * M);

  /* the program's manifold */
  extern Hmesh * hmesh;
  int  manifold (int dim);           /* 0, or -1 (dim out of range) */
  void manifold_free (void);

  /* simplices */
  uint32_t hm_vertex_add (Hmesh * M, const double * x);
  uint32_t hm_add (Hmesh * M, int k, const uint32_t * v);
  uint32_t hm_find (const Hmesh * M, int k, const uint32_t * v);
  int      hm_remove (Hmesh * M, int k, uint32_t s);
  int      hm_remove_closure (Hmesh * M, int k, uint32_t s);

  /*
  .. The box the points stay in: [origin, origin + size] along each axis,
  .. closed. hm_set_domain () refuses (-1, the box unchanged) a size not
  .. positive, or a box that would leave a point in use outside. Then
  .. hm_vertex_add () and hm_vertex_move () refuse a point outside (or NaN),
  .. and hm_check () finds one written outside through hm_x ().
  ..
  .. hm_vertex_move (M, v, x): v moved to x; -1 if v is not in use or x is
  .. outside. hm_touch (M): the points were written through hm_x (), for
  .. the caches of them (M->generation).
  */
  int  hm_set_domain (Hmesh * M, const double * origin, double size);
  int  hm_vertex_move (Hmesh * M, uint32_t v, const double * x);
  void hm_touch (Hmesh * M);

  /*
  .. A manifold crossing the box (hmesh-clip.c)
  .. ------------------------------------------
  ..
  .. hm_clip (M, origin, size): the mesh cut to the box, which then becomes
  .. its domain (hm_set_domain ()). Build the mesh first with its points
  .. anywhere (no domain, or after hm_set_domain () was never called), then
  .. clip; for points that moved out, write them through hm_x () and clip
  .. again. -1 if size is not positive.
  ..
  .. The box is cut by one plane at a time, x = x0, x = x0 + h0, y = y0, ...
  .. (Sutherland and Hodgman's reentrant clipping, on the complex):
  ..
  ..   - a point within HM_CLIP_SNAP size of the plane is moved onto it;
  ..     the others are in or out
  ..   - an edge from in to out gets one new point, where it crosses: the
  ..     edge is stored once, so every simplex on it shares the point and
  ..     the cut mesh stays conforming
  ..   - a maximal simplex (one of level m, or a lower one that is the face
  ..     of none) with points in and out is replaced by the simplices of its
  ..     part inside; one with none in goes; the faces nothing uses then go
  ..     too (hm_remove_closure ())
  ..   - the part inside is triangulated from its points of lowest id: a
  ..     polygon by the fan from its lowest point, a polyhedron by the cone
  ..     from its lowest point to the fans of its faces it is not on. A face
  ..     two tetrahedra share is cut the same on both sides, as both choose
  ..     by the same ids (Dompierre et al. 1999): no point added inside
  ..   - each new simplex in the orientation of the one it comes from: the
  ..     sign of the determinant of its points' barycentric coordinates in
  ..     that simplex
  ..
  .. Fields: a new point's HmReal values are interpolated along its edge
  .. (any other field takes the nearer end's), a new maximal simplex takes the
  .. values of the one it comes from; the other new faces are zero.
  ..
  .. Afterwards the halves on the box are boundary halves, twin (h) = h, next
  .. as ever; hm_half_on_domain () tells them from the manifold's own
  .. boundary by the faces of the box their facet lies on:
  ..
  ..   hm_domain_faces (M, k, s)    the faces of the box every point of s is
  ..                                on: bit 2 d the lower face along d, bit
  ..                                2 d + 1 the upper; 0 without a box
  ..   hm_half_on_domain (M, k, h)  hm_domain_faces () of h's facet if h is
  ..                                a boundary half, else 0
  ..
  .. The points are on the faces exactly (the coordinate is the plane's), so
  .. these are exact.
  */
  #define HM_CLIP_SNAP 1e-10
  int hm_clip (Hmesh * M, const double * origin, double size);
  int hm_domain_faces (const Hmesh * M, int k, uint32_t s);
  int hm_half_on_domain (const Hmesh * M, int k, uint32_t h);

  /*
  .. hm_compact (M, remap): ids 0 ... n - 1 on every level, in the order of
  .. the old ones. With remap not NULL, remap[k] is set to a malloc'ed map
  .. from each old id (< the old hm_end ()) to the new one, HM_NONE for an id
  .. not in use, which the caller frees.
  ..
  .. hm_compact_order (M, order, remap): the same, level k in the order
  .. order[k], the hm_count (M, k) ids in use, each once (order or order[k]
  .. NULL: the order of the ids). -1, and nothing changed, if an order[k] is
  .. not that.
  */
  void hm_compact (Hmesh * M, uint32_t ** remap);
  int  hm_compact_order (Hmesh * M, uint32_t * const * order,
                         uint32_t ** remap);

  /* fields */
  int    hm_field_new (Hmesh * M, int k, const char * name, int ncomp,
                       size_t esize);
  void   hm_field_free (Hmesh * M, int f);
  int    hm_scalar_new (Hmesh * M, int k, const char * name);
  int    hm_vector_new (Hmesh * M, int k, const char * name);

  /* geometry: k-dimensional measure (length, area, volume), centroid */
  double hm_measure (const Hmesh * M, int k, uint32_t s);
  double hm_volume_signed (const Hmesh * M, uint32_t s);   /* k = m = D */
  void   hm_centroid (const Hmesh * M, int k, uint32_t s, double * c);

  /*
  .. hm_check (M): every invariant of the structure, 0 if they hold, else -1
  .. and the first one broken on stderr. For tests: it is O(size).
  */
  int hm_check (const Hmesh * M);

  /* -------- inline accessors -------- */

  static inline uint32_t hm_count (const Hmesh * M, int k)
  { return M->level[k].n; }

  static inline uint32_t hm_end (const Hmesh * M, int k)
  { return M->level[k].end; }

  /* no holes: ids 0 ... hm_count () - 1, a plain loop (for SIMD) */
  static inline int hm_dense (const Hmesh * M, int k)
  { return M->level[k].n == M->level[k].end; }

  static inline int hm_alive (const Hmesh * M, int k, uint32_t s)
  {
    const HmLevel * l = &M->level[k];
    return s < l->end && (l->alive[s >> 6] >> (s & 63) & 1);
  }

  /* the first id in use from s on, or HM_NONE */
  static inline uint32_t hm_alive_from (const Hmesh * M, int k, uint32_t s)
  {
    const HmLevel * l = &M->level[k];
    if (s >= l->end)
      return HM_NONE;
    uint32_t w = s >> 6, nw = (l->end + 63) >> 6;
    uint64_t b = l->alive[w] & (~(uint64_t) 0 << (s & 63));
    while (!b) {
      if (++w >= nw)
        return HM_NONE;
      b = l->alive[w];
    }
    return (w << 6) + (uint32_t) bits_ctz_u64 (b);
  }

  /* every k-simplex in use, in the order of ids */
  #define hm_foreach(M, k, s)                                     \
    for (uint32_t s = hm_alive_from (M, k, 0); s != HM_NONE;      \
         s = hm_alive_from (M, k, s + 1))

  static inline uint32_t hm_vertex (const Hmesh * M, int k, uint32_t s, int i)
  { return k ? M->level[k].vert[(uint32_t) (k + 1) * s + (uint32_t) i] : s; }

  static inline uint32_t hm_half (int k, uint32_t s, int i)
  { return (uint32_t) (k + 1) * s + (uint32_t) i; }

  static inline uint32_t hm_half_simplex (int k, uint32_t h)
  { return h / (uint32_t) (k + 1); }

  static inline int hm_half_index (int k, uint32_t h)
  { return (int) (h % (uint32_t) (k + 1)); }

  static inline uint32_t hm_half_next (int k, uint32_t h)
  {
    uint32_t n = (uint32_t) (k + 1), i = h % n;
    return i + 1 == n ? h - i : h + 1;
  }

  static inline uint32_t hm_half_twin (const Hmesh * M, int k, uint32_t h)
  { return M->level[k].twin[h]; }

  static inline int hm_half_boundary (const Hmesh * M, int k, uint32_t h)
  { return M->level[k].twin[h] == h; }

  static inline uint32_t hm_half_vertex (const Hmesh * M, int k, uint32_t h)
  { return M->level[k].vert[h]; }

  static inline uint32_t hm_half_facet (const Hmesh * M, int k, uint32_t h)
  { return M->level[k].facet[h] & ~HM_FLIP; }

  static inline int hm_half_sign (const Hmesh * M, int k, uint32_t h)
  { return M->level[k].facet[h] & HM_FLIP ? -1 : 1; }

  /* a half-edge's origin, k = 2: v(i+1) */
  static inline uint32_t hm_half_origin (const Hmesh * M, uint32_t h)
  { return M->level[2].vert[hm_half_next (2, h)]; }

  static inline uint32_t hm_use (const Hmesh * M, int k, uint32_t s)
  { return M->level[k].use[s]; }

  static inline void * hm_field (const Hmesh * M, int f, int c)
  { return M->field[f].data[c]; }

  static inline HmReal * hm_val (const Hmesh * M, int f, int c)
  { return (HmReal *) M->field[f].data[c]; }

  static inline double * hm_x (const Hmesh * M, int c)
  { return M->x[c]; }

  /* x in the box, or no box; a NaN is not in it */
  static inline int hm_inside (const Hmesh * M, const double * x)
  {
    if (!M->domain)
      return 1;
    for (int d = 0; d < HM_DIMENSION; d++)
      if (!(x[d] >= M->origin[d] && x[d] <= M->origin[d] + M->size))
        return 0;
    return 1;
  }

#endif
