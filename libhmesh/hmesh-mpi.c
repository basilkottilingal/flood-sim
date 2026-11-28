/*
.. The partitioned mesh, stage 1: the data. See hmesh-mpi.h, and
.. doc/hmesh.tex, "Partitioned hmesh: design". Built with cwrap on
.. libb-amr's tree, -DEULERIAN_LAGRANGIAN; nothing without it.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "tree-mpi.h"
#include "block-ghosts.h"
#include "hmesh-mpi.h"
#include "hmesh-vof.h"
#ifdef _MPI
#include <mpi.h>
#endif

#ifdef EULERIAN_LAGRANGIAN

#define D HM_DIMENSION

struct HmPart {
  Hmesh  * M;
  int      rank, npe;
  int      gid;                    /* field: level 0, uint64_t */
  int      owner[HM_LEVELS];       /* fields: int32_t a simplex */
  int      nf[HM_LEVELS];          /* fields, k < m: uint32_t */
  uint32_t nowned[HM_LEVELS];
  uint64_t counter;                /* the next point made here */
  uint64_t * K;                    /* the partition's first keys at the
                                      last refresh: a change of them is a
                                      repartition */
  /* global id -> local id: open addressing, UINT64_MAX empty */
  uint64_t * hkey;
  uint32_t * hval;
  uint32_t   hcap;
};

static void save_K (HmPart * P);

static void * grab (size_t n, size_t e)
{
  void * p = calloc (n ? n : 1, e);
  if (!p) {
    fprintf (stderr, "hmesh-mpi: out of memory\n");
    abort ();
  }
  return p;
}

static uint64_t * gids (const HmPart * P)
{
  return hm_field (P->M, P->gid, 0);
}

static int32_t * owners (const HmPart * P, int k)
{
  return hm_field (P->M, P->owner[k], 0);
}

static uint32_t * cofaces (const HmPart * P, int k)
{
  return hm_field (P->M, P->nf[k], 0);
}

static void point_of (const Hmesh * M, uint32_t v, double * x)
{
  for (int d = 0; d < D; d++)
    x[d] = M->x[d][v];
}

/* -------- the tree: indices, keys, the owners' ranges -------- */

static Index index_of (int level, const long * i)
{
  Index p;
  memset (&p, 0, sizeof (p));
  p.level = level;
  p.i = (int) i[0];
#if D >= 2
  p.j = (int) i[1];
#endif
#if D >= 3
  p.k = (int) i[2];
#endif
  return p;
}

static int owner_of_key (const TreePartition * T, uint64_t key)
{
  int lo = 0, hi = T->npe;             /* key_first [lo] <= key < [hi] */
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (T->key_first[mid] <= key)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

/* the key of the node of the deepest level at x */
static uint64_t point_key (const double * x)
{
  int L = tree_get_depth_max ();
  double h0 = tree_get_size ();
  long n = 1L << L, i[3] = { 0, 0, 0 };
  for (int d = 0; d < D; d++) {
    double t = floor ((x[d] - tree_get_origin (d)) / h0 * (double) n);
    i[d] = t < 0. ? 0 : t >= (double) n ? n - 1 : (long) t;
  }
  return tree_morton_key (index_of (L, i));
}

int hm_part_point_owner (const HmPart * P, const double * x)
{
  const TreePartition * T = tree_get_partition ();
  if (!T)
    return P->rank;
  return owner_of_key (T, point_key (x));
}

/* the box [lo, hi] meets the node (level, i) and its key range is
   partly or wholly this process's */
static int meets (int rank, const TreePartition * T, int level,
                  const long * i, const double * lo, const double * hi)
{
  double h0 = tree_get_size (), w = h0 / (double) (1L << level);
  for (int d = 0; d < D; d++) {
    double a = tree_get_origin (d) + (double) i[d] * w;
    if (hi[d] < a || lo[d] > a + w)
      return 0;
  }
  uint64_t s = tree_morton_key (index_of (level, i));
  uint64_t span = (uint64_t) 1 << (D * (BAMR_MAX_DEPTH - level));
  uint64_t a = T->key_first[rank], b = T->key_first[rank + 1];
  if (s >= a && (b == UINT64_MAX || s + span <= b))
    return 1;                          /* wholly this process's */
  if (s + span <= a || s >= b)
    return 0;                          /* wholly another's */
  if (level >= tree_get_depth_max ())
    return s >= a && s < b;
  for (int c = 0; c < 1 << D; c++) {
    long j[3] = { 0, 0, 0 };
    for (int d = 0; d < D; d++)
      j[d] = 2 * i[d] + (c >> d & 1);
    if (meets (rank, T, level + 1, j, lo, hi))
      return 1;
  }
  return 0;
}

/* the box meets one of process rank's leaves */
static int box_meets_rank (const HmPart * P, int rank, const double * lo,
                           const double * hi)
{
  const TreePartition * T = tree_get_partition ();
  if (!T) {                            /* every leaf this process's */
    if (rank != P->rank)
      return 0;
    double h0 = tree_get_size ();
    for (int d = 0; d < D; d++) {
      double o = tree_get_origin (d);
      if (hi[d] < o || lo[d] > o + h0)
        return 0;
    }
    return 1;
  }
  long root[3] = { 0, 0, 0 };
  return meets (rank, T, 0, root, lo, hi);
}

int hm_part_box_meets (const HmPart * P, const double * lo,
                       const double * hi)
{
  return box_meets_rank (P, P->rank, lo, hi);
}

/* -------- global ids -------- */

static uint32_t hslot (const HmPart * P, uint64_t g)
{
  return (uint32_t) ((g * 0x9E3779B97F4A7C15ull) >> 32) & (P->hcap - 1);
}

static void hash_build (HmPart * P)
{
  Hmesh * M = P->M;
  uint32_t n = hm_count (M, 0), cap = 16;
  while (cap < 2 * n)
    cap *= 2;
  free (P->hkey), free (P->hval);
  P->hkey = grab (cap, sizeof (uint64_t));
  P->hval = grab (cap, sizeof (uint32_t));
  P->hcap = cap;
  memset (P->hkey, 0xff, cap * sizeof (uint64_t));
  const uint64_t * g = gids (P);
  hm_foreach (M, 0, v) {
    uint32_t a = hslot (P, g[v]);
    while (P->hkey[a] != UINT64_MAX)
      a = (a + 1) & (P->hcap - 1);
    P->hkey[a] = g[v], P->hval[a] = v;
  }
}

uint32_t hm_part_local (const HmPart * P, uint64_t g)
{
  if (!P->hcap || g == UINT64_MAX)
    return HM_NONE;
  uint32_t a = hslot (P, g);
  while (P->hkey[a] != UINT64_MAX) {
    if (P->hkey[a] == g)
      return P->hval[a];
    a = (a + 1) & (P->hcap - 1);
  }
  return HM_NONE;
}

uint64_t hm_part_gid (const HmPart * P, uint32_t v)
{
  return gids (P)[v];
}

int hm_part_owner (const HmPart * P, int k, uint32_t s)
{
  return owners (P, k)[s];
}

uint32_t hm_part_owned (const HmPart * P, int k)
{
  return P->nowned[k];
}

uint32_t hm_part_cofaces (const HmPart * P, int k, uint32_t f)
{
  if (k < 0 || k >= P->M->dim)
    return 0;
  return cofaces (P, k)[f];
}

/* the cofaces of f held here: the length of its cycle at level k + 1 */
static uint32_t held_cofaces (const Hmesh * M, int k, uint32_t f)
{
  uint32_t h = hm_use (M, k, f), n = 0;
  if (h == HM_NONE)
    return 0;
  uint32_t g = h;
  do {
    n++;
    g = hm_half_twin (M, k + 1, g);
  } while (g != h);
  return n;
}

int hm_part_half_cut (const HmPart * P, int k, uint32_t h)
{
  if (k < 1 || k > P->M->dim)
    return 0;
  uint32_t f = hm_half_facet (P->M, k, h);
  return held_cofaces (P->M, k - 1, f) < cofaces (P, k - 1)[f];
}

/* -------- the partition -------- */

HmPart * hm_part_new (Hmesh * M)
{
  if (!M || !M->domain || !M->cell_size)
    return NULL;
  HmPart * P = grab (1, sizeof (HmPart));
  P->M = M;
  P->rank = tree_rank ();
  P->npe = tree_npe ();
  P->gid = hm_field_new (M, 0, "gid", 1, sizeof (uint64_t));
  for (int k = 0; k < HM_LEVELS; k++)
    P->owner[k] = P->nf[k] = -1;
  for (int k = 0; k <= M->dim; k++) {
    P->owner[k] = hm_field_new (M, k, "owner", 1, sizeof (int32_t));
    if (k < M->dim)
      P->nf[k] = hm_field_new (M, k, "cofaces", 1, sizeof (uint32_t));
  }
  uint64_t * g = gids (P);
  hm_foreach (M, 0, v)
    g[v] = v;                          /* r0 = 0: made alike everywhere */
  P->counter = P->rank == 0 ? hm_end (M, 0) : 0;
  for (int k = 0; k <= M->dim; k++)
    P->nowned[k] = hm_end (M, k);
  hash_build (P);
  return P;
}

void hm_part_free (HmPart * P)
{
  if (!P)
    return;
  Hmesh * M = P->M;
  for (int k = 0; k <= M->dim; k++) {
    hm_field_free (M, P->owner[k]);
    if (P->nf[k] >= 0)
      hm_field_free (M, P->nf[k]);
  }
  hm_field_free (M, P->gid);
  free (P->hkey), free (P->hval), free (P->K);
  free (P);
}

/* owners by the min rule, from the points' */
static void set_owners (HmPart * P)
{
  Hmesh * M = P->M;
  int32_t * o0 = owners (P, 0);
  hm_foreach (M, 0, v) {
    double x[3];
    point_of (M, v, x);
    o0[v] = hm_part_point_owner (P, x);
  }
  for (int k = 1; k <= M->dim; k++) {
    int32_t * o = owners (P, k);
    o0 = owners (P, 0);
    hm_foreach (M, k, s) {
      int32_t r = INT32_MAX;
      for (int i = 0; i <= k; i++) {
        int32_t q = o0[hm_vertex (M, k, s, i)];
        r = q < r ? q : r;
      }
      o[s] = r;
    }
  }
}

static int compact_owned (HmPart * P);

/* s in A_r: owned, or its box grown by eps meets one of r's leaves; *far
   counts owned ones with a point in no leaf here */
static int reached (HmPart * P, int k, uint32_t s, int * far)
{
  Hmesh * M = P->M;
  int own = owners (P, k)[s] == P->rank;
  double lo[3], hi[3], hmax = 0.;
  for (int d = 0; d < D; d++)
    lo[d] = HUGE_VAL, hi[d] = -HUGE_VAL;
  for (int i = 0; i <= k; i++) {
    double x[3];
    point_of (M, hm_vertex (M, k, s, i), x);
    double h = hm_tree_cell_size (M, x);
    if (h <= 0.) {                     /* in no leaf here */
      if (own)
        (*far)++;
      return own;
    }
    hmax = fmax (hmax, h);
    for (int d = 0; d < D; d++)
      lo[d] = fmin (lo[d], x[d]), hi[d] = fmax (hi[d], x[d]);
  }
  if (own)
    return 1;
  for (int d = 0; d < D; d++)
    lo[d] -= 0.5 * hmax, hi[d] += 0.5 * hmax;
  return hm_part_box_meets (P, lo, hi);
}

int hm_part_distribute (HmPart * P)
{
  Hmesh * M = P->M;
  int m = M->dim, far = 0;
  set_owners (P);
  /* the cofaces in the whole mesh */
  for (int k = 0; k < m; k++) {
    uint32_t * nf = cofaces (P, k);
    hm_foreach (M, k, f)
      nf[f] = held_cofaces (M, k, f);
  }
  /* A_r, then its closure */
  char * held[HM_LEVELS] = { NULL };
  for (int k = 0; k <= m; k++) {
    held[k] = grab (hm_end (M, k), 1);
    hm_foreach (M, k, s)
      held[k][s] = (char) reached (P, k, s, &far);
  }
  for (int k = m; k >= 1; k--)
    hm_foreach (M, k, s)
      if (held[k][s])
        for (int i = 0; i <= k; i++)
          held[k - 1][hm_half_facet (M, k, hm_half (k, s, i))] = 1;
  /* the rest removed, the top first */
  int bad = 0;
  for (int k = m; k >= 0; k--)
    hm_foreach (M, k, s)
      if (!held[k][s] && hm_remove (M, k, s))
        bad++;
  for (int k = 0; k <= m; k++)
    free (held[k]);
  if (compact_owned (P))
    bad++;
  save_K (P);
  if (bad)
    fprintf (stderr, "hmesh-mpi: %d removals or compactions refused\n",
             bad);
  return far || bad ? -1 : 0;
}

/* every level renumbered, owned first; the map made again */
static int compact_owned (HmPart * P)
{
  Hmesh * M = P->M;
  int m = M->dim, bad = 0;
  uint32_t * order[HM_LEVELS] = { NULL };
  for (int k = 0; k <= m; k++) {
    order[k] = grab (hm_count (M, k), sizeof (uint32_t));
    const int32_t * o = owners (P, k);
    uint32_t n = 0;
    hm_foreach (M, k, s)
      if (o[s] == P->rank)
        order[k][n++] = s;
    P->nowned[k] = n;
    hm_foreach (M, k, s)
      if (o[s] != P->rank)
        order[k][n++] = s;
  }
  if (hm_compact_order (M, order, NULL))
    bad++;
  for (int k = 0; k <= m; k++)
    free (order[k]);
  hash_build (P);
  return bad;
}


/* -------- stage 2: a step -------- */

typedef struct {
  unsigned char * p;
  size_t n, cap, at;               /* at: the reading position */
} Buf;

static void put (Buf * b, const void * x, size_t n)
{
  if (b->n + n > b->cap) {
    size_t cap = b->cap ? 2 * b->cap : 4096;
    while (cap < b->n + n)
      cap *= 2;
    unsigned char * q = realloc (b->p, cap);
    if (!q) {
      fprintf (stderr, "hmesh-mpi: out of memory\n");
      abort ();
    }
    b->p = q, b->cap = cap;
  }
  memcpy (b->p + b->n, x, n);
  b->n += n;
}

static void get (Buf * b, void * x, size_t n)
{
  memcpy (x, b->p + b->at, n);
  b->at += n;
}

/* the peers: the processes this one exchanges with (libb-amr's plan) */
static int peers_of (const int ** peers)
{
#ifdef _MPI
  const TreeMpiPlan * plan = tree_mpi_plan ();
  if (plan) {
    *peers = plan->peers;
    return plan->npeers;
  }
#endif
  *peers = NULL;
  return 0;
}

/* send[i] to peers[i], recv[i] from it: sizes first, then the bytes */
static void exchange (int np, const int * peers, Buf * send, Buf * recv,
                      int tag)
{
#ifdef _MPI
  uint64_t * ns = grab ((size_t) np + 1, sizeof (uint64_t));
  uint64_t * nr = grab ((size_t) np + 1, sizeof (uint64_t));
  MPI_Request * rq = grab (2 * (size_t) np + 1, sizeof (MPI_Request));
  for (int i = 0; i < np; i++) {
    ns[i] = send[i].n;
    MPI_Irecv (&nr[i], 8, MPI_BYTE, peers[i], tag, MPI_COMM_WORLD, &rq[i]);
    MPI_Isend (&ns[i], 8, MPI_BYTE, peers[i], tag, MPI_COMM_WORLD,
               &rq[np + i]);
  }
  MPI_Waitall (2 * np, rq, MPI_STATUSES_IGNORE);
  int nq = 0;
  for (int i = 0; i < np; i++) {
    recv[i].n = recv[i].cap = nr[i], recv[i].at = 0;
    recv[i].p = grab (nr[i] + 1, 1);
    if (nr[i])
      MPI_Irecv (recv[i].p, (int) nr[i], MPI_BYTE, peers[i], tag + 1,
                 MPI_COMM_WORLD, &rq[nq++]);
    if (ns[i])
      MPI_Isend (send[i].p, (int) ns[i], MPI_BYTE, peers[i], tag + 1,
                 MPI_COMM_WORLD, &rq[nq++]);
  }
  MPI_Waitall (nq, rq, MPI_STATUSES_IGNORE);
  free (ns), free (nr), free (rq);
#else
  (void) peers, (void) send, (void) recv, (void) tag, (void) np;
#endif
}

static int reduce_min (int x)
{
#ifdef _MPI
  int y = x;
  MPI_Allreduce (&x, &y, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return y;
#else
  return x;
#endif
}

static long reduce_sum (long x)
{
#ifdef _MPI
  long y = x;
  MPI_Allreduce (&x, &y, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  return y;
#else
  return x;
#endif
}

/* the user's fields on level k: every field but the partition's */
static int is_ours (const HmPart * P, int f)
{
  if (f == P->gid)
    return 1;
  for (int k = 0; k < HM_LEVELS; k++)
    if (f == P->owner[k] || f == P->nf[k])
      return 1;
  return 0;
}

static size_t field_bytes (const HmPart * P, int k)
{
  const Hmesh * M = P->M;
  size_t n = 0;
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k && !is_ours (P, f))
      n += M->field[f].esize * (size_t) M->field[f].ncomp;
  return n;
}

static void put_fields (const HmPart * P, int k, uint32_t s, Buf * b)
{
  const Hmesh * M = P->M;
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k && !is_ours (P, f))
      for (int c = 0; c < M->field[f].ncomp; c++)
        put (b, (char *) M->field[f].data[c] + s * M->field[f].esize,
             M->field[f].esize);
}

static void get_fields (const HmPart * P, int k, uint32_t s, Buf * b)
{
  const Hmesh * M = P->M;
  for (int f = 0; f < M->nfield; f++)
    if (M->field[f].level == k && !is_ours (P, f))
      for (int c = 0; c < M->field[f].ncomp; c++)
        get (b, (char *) M->field[f].data[c] + s * M->field[f].esize,
             M->field[f].esize);
}

/* the box of simplex s grown by half the largest cell at its points */
static int grown_box (const HmPart * P, int k, uint32_t s, double * lo,
                      double * hi)
{
  Hmesh * M = P->M;
  double hmax = 0.;
  for (int d = 0; d < D; d++)
    lo[d] = HUGE_VAL, hi[d] = -HUGE_VAL;
  for (int i = 0; i <= k; i++) {
    double x[3];
    point_of (M, hm_vertex (M, k, s, i), x);
    double h = hm_tree_cell_size (M, x);
    if (h <= 0.)
      return -1;
    hmax = fmax (hmax, h);
    for (int d = 0; d < D; d++)
      lo[d] = fmin (lo[d], x[d]), hi[d] = fmax (hi[d], x[d]);
  }
  for (int d = 0; d < D; d++)
    lo[d] -= 0.5 * hmax, hi[d] += 0.5 * hmax;
  return 0;
}

/* -------- entries: a simplex as its points' global ids -------- */

/* k, the global ids in the stored order, the position (a point), n_f
   (below the top), the fields */
static void put_entry (const HmPart * P, Buf * b, int k, uint32_t s)
{
  Hmesh * M = P->M;
  put (b, &k, sizeof (int));
  for (int j = 0; j <= k; j++)
    put (b, &gids (P)[hm_vertex (M, k, s, j)], sizeof (uint64_t));
  if (k == 0)
    for (int d = 0; d < D; d++)
      put (b, &M->x[d][s], sizeof (double));
  if (k < M->dim)
    put (b, &cofaces (P, k)[s], sizeof (uint32_t));
  put_fields (P, k, s, b);
}

static size_t entry_bytes (const HmPart * P, int k)
{
  return sizeof (int) + (size_t) (k + 1) * sizeof (uint64_t) +
    (k == 0 ? D * sizeof (double) : 0) +
    (k < P->M->dim ? sizeof (uint32_t) : 0) + field_bytes (P, k);
}

/* s and every face of it */
static void put_closure (const HmPart * P, Buf * b, int k, uint32_t s)
{
  Hmesh * M = P->M;
  put_entry (P, b, k, s);
  if (k == 0)
    return;
  uint32_t v[HM_LEVELS];
  for (int j = 0; j <= k; j++)
    v[j] = hm_vertex (M, k, s, j);
  /* every face: the subsets of the points, by bit masks */
  for (unsigned mask = 1; mask < (1u << (k + 1)) - 1; mask++) {
    uint32_t u[HM_LEVELS];
    int n = 0;
    for (int j = 0; j <= k; j++)
      if (mask >> j & 1)
        u[n++] = v[j];
    uint32_t f = n == 1 ? u[0] : hm_find (M, n - 1, u);
    if (f != HM_NONE)
      put_entry (P, b, n - 1, f);
  }
}

typedef struct {
  uint32_t * id;
  uint32_t   n, cap;
} Ids;

static void ids_push (Ids * L, uint32_t s)
{
  if (L->n == L->cap) {
    L->cap = L->cap ? 2 * L->cap : 1024;
    L->id = realloc (L->id, L->cap * sizeof (uint32_t));
    if (!L->id) {
      fprintf (stderr, "hmesh-mpi: out of memory\n");
      abort ();
    }
  }
  L->id[L->n++] = s;
}

/* every entry of the buffers, a level at a time (faces before cofaces,
   in the sender's orientation); seen[k] the local ids met. The number of
   entries that could not be placed */
static int insert_entries (HmPart * P, int nb, Buf * recv, Ids * seen)
{
  Hmesh * M = P->M;
  int m = M->dim, bad = 0;
  M->no_edge_check = 1;              /* edges stretched by a move */
  for (int k = 0; k <= m; k++)
    for (int i = 0; i < nb; i++) {
      Buf * b = &recv[i];
      b->at = 0;
      while (b->at < b->n) {
        int kk;
        get (b, &kk, sizeof (int));
        if (kk != k) {
          b->at += entry_bytes (P, kk) - sizeof (int);
          continue;
        }
        uint64_t g[HM_LEVELS];
        uint32_t v[HM_LEVELS], s, nf = 0;
        for (int j = 0; j <= k; j++)
          get (b, &g[j], sizeof (uint64_t));
        if (k == 0) {
          double x[3] = { 0, 0, 0 };
          for (int d = 0; d < D; d++)
            get (b, &x[d], sizeof (double));
          s = hm_part_local (P, g[0]);
          if (s == HM_NONE) {
            s = hm_vertex_add (M, x);
            if (s == HM_NONE) {
              bad++;
              b->at += entry_bytes (P, 0) - sizeof (int) -
                sizeof (uint64_t) - D * sizeof (double);
              continue;
            }
            gids (P)[s] = g[0];
            if (2 * (hm_count (M, 0) + 1) > P->hcap)
              hash_build (P);
            else {
              uint32_t h = hslot (P, g[0]);
              while (P->hkey[h] != UINT64_MAX)
                h = (h + 1) & (P->hcap - 1);
              P->hkey[h] = g[0], P->hval[h] = s;
            }
          }
          else
            for (int d = 0; d < D; d++)
              M->x[d][s] = x[d];
        }
        else {
          int ok = 1;
          for (int j = 0; j <= k; j++)
            ok &= (v[j] = hm_part_local (P, g[j])) != HM_NONE;
          s = ok ? hm_find (M, k, v) : HM_NONE;
          if (ok && s == HM_NONE)
            s = hm_add (M, k, v);
          if (s == HM_NONE) {
            bad++;
            b->at += entry_bytes (P, k) - sizeof (int) -
              (size_t) (k + 1) * sizeof (uint64_t);
            continue;
          }
        }
        if (k < m) {
          get (b, &nf, sizeof (uint32_t));
          cofaces (P, k)[s] = nf;
        }
        get_fields (P, k, s, b);
        if (seen)
          ids_push (&seen[k], s);
      }
    }
  M->no_edge_check = 0;
  hm_touch (M);
  return bad;
}

/* S5 and S6: the refresh. Who must hold each simplex: the peers whose
   leaves its grown box meets, or that own it (A_r), closed under faces
   from the top down; each owner sends what it owns; the receivers add,
   and drop what they neither own nor were sent; owned first again */
static int refresh (HmPart * P)
{
  Hmesh * M = P->M;
  int m = M->dim, bad = 0;
  const int * peers;
  int np = peers_of (&peers);
  if (np > 64) {
    fprintf (stderr, "hmesh-mpi: more than 64 peers\n");
    return 1;
  }
  uint64_t * mask[HM_LEVELS] = { NULL };
  for (int k = m; k >= 0; k--) {
    mask[k] = grab (hm_end (M, k), sizeof (uint64_t));
    if (!np)
      continue;
    const int32_t * o = owners (P, k);
    hm_foreach (M, k, s) {
      double lo[3], hi[3];
      int boxed = grown_box (P, k, s, lo, hi) == 0;
      for (int i = 0; i < np; i++)
        if (o[s] == peers[i] ||
            (boxed && box_meets_rank (P, peers[i], lo, hi)))
          mask[k][s] |= (uint64_t) 1 << i;
    }
  }
  for (int k = m; k >= 1; k--)
    hm_foreach (M, k, s)
      for (int i = 0; i <= k; i++)
        mask[k - 1][hm_half_facet (M, k, hm_half (k, s, i))] |= mask[k][s];
  Buf * send = grab ((size_t) np + 1, sizeof (Buf));
  Buf * recv = grab ((size_t) np + 1, sizeof (Buf));
  for (int k = 0; k <= m; k++) {
    const int32_t * o = owners (P, k);
    hm_foreach (M, k, s) {
      if (o[s] != P->rank)
        continue;
      if (k < m)
        cofaces (P, k)[s] = held_cofaces (M, k, s);   /* holds its star */
      for (int i = 0; i < np; i++)
        if (mask[k][s] >> i & 1)
          put_entry (P, &send[i], k, s);
    }
  }
  for (int k = 0; k <= m; k++)
    free (mask[k]);
  exchange (np, peers, send, recv, 300);
  Ids seen[HM_LEVELS];
  memset (seen, 0, sizeof (seen));
  bad += insert_entries (P, np, recv, seen);
  for (int i = 0; i < np; i++)
    free (send[i].p), free (recv[i].p);
  free (send), free (recv);
  set_owners (P);
  for (int k = m; k >= 0; k--) {
    char * keep = grab (hm_end (M, k), 1);
    for (uint32_t a = 0; a < seen[k].n; a++)
      keep[seen[k].id[a]] = 1;
    const int32_t * o = owners (P, k);
    hm_foreach (M, k, s)
      if (o[s] != P->rank && !keep[s] && hm_remove (M, k, s))
        bad++;
    free (keep);
    free (seen[k].id);
  }
  bad += compact_owned (P);
  return bad;
}

/* -------- stage 3: the tree changed -------- */

/* the partition's first keys, as at the last refresh */
static void save_K (HmPart * P)
{
  const TreePartition * T = tree_get_partition ();
  free (P->K);
  P->K = NULL;
  if (!T)
    return;
  P->K = grab ((size_t) T->npe + 1, sizeof (uint64_t));
  memcpy (P->K, T->key_first, ((size_t) T->npe + 1) * sizeof (uint64_t));
}

static int K_changed (const HmPart * P)
{
  const TreePartition * T = tree_get_partition ();
  if (!T || !P->K)
    return T != NULL || P->K != NULL;
  return memcmp (P->K, T->key_first, ((size_t) T->npe + 1) *
                 sizeof (uint64_t)) != 0;
}

/* the owner of x under the saved keys */
static int old_owner (const HmPart * P, const double * x)
{
  if (!P->K)
    return P->rank;
  TreePartition T;
  T.npe = P->npe;
  T.key_first = P->K;
  return owner_of_key (&T, point_key (x));
}

/* buffers to any process: counts by MPI_Alltoall, then point to point
   where there is something; send and recv have npe entries */
static void exchange_all (HmPart * P, Buf * send, Buf * recv, int tag)
{
#ifdef _MPI
  int npe = P->npe;
  long long * ns = grab ((size_t) npe, sizeof (long long));
  long long * nr = grab ((size_t) npe, sizeof (long long));
  for (int r = 0; r < npe; r++)
    ns[r] = (long long) send[r].n;
  MPI_Alltoall (ns, 1, MPI_LONG_LONG, nr, 1, MPI_LONG_LONG,
                MPI_COMM_WORLD);
  MPI_Request * rq = grab (2 * (size_t) npe, sizeof (MPI_Request));
  int nq = 0;
  for (int r = 0; r < npe; r++) {
    recv[r].n = recv[r].cap = (size_t) nr[r], recv[r].at = 0;
    recv[r].p = grab ((size_t) nr[r] + 1, 1);
    if (r == P->rank) {
      if (send[r].n)
        memcpy (recv[r].p, send[r].p, send[r].n);
      continue;
    }
    if (nr[r])
      MPI_Irecv (recv[r].p, (int) nr[r], MPI_BYTE, r, tag, MPI_COMM_WORLD,
                 &rq[nq++]);
    if (ns[r])
      MPI_Isend (send[r].p, (int) ns[r], MPI_BYTE, r, tag, MPI_COMM_WORLD,
                 &rq[nq++]);
  }
  MPI_Waitall (nq, rq, MPI_STATUSES_IGNORE);
  free (ns), free (nr), free (rq);
#else
  recv[0] = send[0];
  send[0].p = NULL;
  (void) tag;
#endif
}

int hm_part_tree_changed (HmPart * P)
{
  Hmesh * M = P->M;
  int m = M->dim, bad = 0;
  if (K_changed (P)) {
    /* R1: the old owners of the points, then the new ones */
    int32_t * was = grab (hm_end (M, 0), sizeof (int32_t));
    hm_foreach (M, 0, v) {
      double x[3];
      point_of (M, v, x);
      was[v] = old_owner (P, x);
    }
    set_owners (P);
    /* R2: each simplex, with its closure, from the old owner of its
       carrier (its first point owned by the new owner) to the new owner */
    int npe = P->npe;
    Buf * send = grab ((size_t) npe, sizeof (Buf));
    Buf * recv = grab ((size_t) npe, sizeof (Buf));
    const int32_t * o0 = owners (P, 0);
    for (int k = 0; k <= m; k++) {
      const int32_t * o = owners (P, k);
      hm_foreach (M, k, s) {
        int q = o[s];
        if (q == P->rank)
          continue;
        int j = 0;
        while (j < k && o0[hm_vertex (M, k, s, j)] != q)
          j++;
        if (was[hm_vertex (M, k, s, j)] == P->rank)
          put_closure (P, &send[q], k, s);
      }
    }
    free (was);
    exchange_all (P, send, recv, 400);
    bad += insert_entries (P, npe, recv, NULL);
    for (int r = 0; r < npe; r++) {
      free (send[r].p), free (recv[r].p);
      memset (&send[r], 0, sizeof (Buf)), memset (&recv[r], 0, sizeof (Buf));
    }
    free (send), free (recv);
    set_owners (P);
    /* R3: each maximal simplex owned here, with its closure, to the new
       peers its grown box reaches: then every owner holds its stars */
    const int * peers;
    int np = peers_of (&peers);
    send = grab ((size_t) np + 1, sizeof (Buf));
    recv = grab ((size_t) np + 1, sizeof (Buf));
    for (int k = 0; k <= m; k++) {
      const int32_t * o = owners (P, k);
      hm_foreach (M, k, s) {
        if (o[s] != P->rank || (k < m && hm_use (M, k, s) != HM_NONE))
          continue;
        double lo[3], hi[3];
        if (grown_box (P, k, s, lo, hi))
          continue;
        for (int i = 0; i < np; i++)
          if (box_meets_rank (P, peers[i], lo, hi))
            put_closure (P, &send[i], k, s);
      }
    }
    exchange (np, peers, send, recv, 500);
    bad += insert_entries (P, np, recv, NULL);
    for (int i = 0; i < np; i++)
      free (send[i].p), free (recv[i].p);
    free (send), free (recv);
    set_owners (P);
  }
  /* R4, or the refresh after an adaptation (K kept) */
  bad += refresh (P);
  save_K (P);
  return reduce_sum (bad) ? -3 : 0;
}

int hm_part_advect (HmPart * P, face_vector uf, scalar c, double dt,
                    face_vector uf_prev, double dt_prev)
{
  Hmesh * M = P->M;
  int m = M->dim;
  const int * peers;
  int np = peers_of (&peers);
  if (np > 64) {
    fprintf (stderr, "hmesh-mpi: more than 64 peers\n");
    return -2;
  }

  /* S1: the owned points' displacements; all move, or none */
  uint32_t end0 = hm_end (M, 0);
  char * which = grab (end0, 1);
  const int32_t * o0 = owners (P, 0);
  hm_foreach (M, 0, v)
    which[v] = o0[v] == P->rank;
  HmAdvect * A = hm_advect_begin (M, uf, dt, uf_prev, dt_prev, which);
  int st = reduce_min (hm_advect_status (A));
  if (st) {
    hm_advect_free (A);
    free (which);
    return st;
  }
  /* S2's lists, from the leaves before the move */
  const HmTreeCache * TC = hm_tree_cache (M);
  const TreeCaches * tc = tree_caches ();
  Buf * sendlist = grab ((size_t) np + 1, sizeof (Buf));
  hm_foreach (M, 0, v) {
    if (!which[v] || !np)
      continue;
    Index p = tc->leaves[TC->leaf_of[v]].p;
    double w = tree_get_size () / (double) (1L << p.level);
    double lo[3], hi[3];
    long ijk[3] = { p.i, 0, 0 };
#if D >= 2
    ijk[1] = p.j;
#endif
#if D >= 3
    ijk[2] = p.k;
#endif
    for (int d = 0; d < D; d++)
      lo[d] = tree_get_origin (d) + (double) ijk[d] * w, hi[d] = lo[d] + w;
    for (int i = 0; i < np; i++)
      if (box_meets_rank (P, peers[i], lo, hi))
        put (&sendlist[i], &v, sizeof (uint32_t));
  }
  hm_advect_move (A);

  /* S2: the new positions to the peers, exact */
  Buf * send = grab ((size_t) np + 1, sizeof (Buf));
  Buf * recv = grab ((size_t) np + 1, sizeof (Buf));
  for (int i = 0; i < np; i++) {
    for (size_t a = 0; a < sendlist[i].n; a += sizeof (uint32_t)) {
      uint32_t v;
      memcpy (&v, sendlist[i].p + a, sizeof (v));
      uint64_t g = gids (P)[v];
      put (&send[i], &g, sizeof (g));
      for (int d = 0; d < D; d++)
        put (&send[i], &M->x[d][v], sizeof (double));
    }
    free (sendlist[i].p);
  }
  free (sendlist);
  exchange (np, peers, send, recv, 100);
  for (int i = 0; i < np; i++) {
    while (recv[i].at < recv[i].n) {
      uint64_t g;
      double x[3];
      get (&recv[i], &g, sizeof (g));
      for (int d = 0; d < D; d++)
        get (&recv[i], &x[d], sizeof (double));
      uint32_t v = hm_part_local (P, g);
      if (v != HM_NONE)
        for (int d = 0; d < D; d++)
          M->x[d][v] = x[d];
    }
    free (send[i].p), free (recv[i].p);
    memset (&send[i], 0, sizeof (Buf)), memset (&recv[i], 0, sizeof (Buf));
  }
  hm_touch (M);
  free (which);

  /* S3: owners by the min rule, the old ones kept for S4 */
  int32_t * old[HM_LEVELS] = { NULL };
  for (int k = 0; k <= m; k++) {
    old[k] = grab (hm_end (M, k), sizeof (int32_t));
    memcpy (old[k], owners (P, k), hm_end (M, k) * sizeof (int32_t));
  }
  set_owners (P);

  /* S4: the fields of the simplices changing hands, to their new owners
     (who hold them already) */
  int bad = 0;
  for (int k = 0; k <= m; k++) {
    const int32_t * o = owners (P, k);
    hm_foreach (M, k, s) {
      if (old[k][s] != P->rank || o[s] == P->rank)
        continue;
      int i = 0;
      while (i < np && peers[i] != o[s])
        i++;
      if (i == np) {
        bad++;
        continue;
      }
      put (&send[i], &k, sizeof (int));
      for (int j = 0; j <= k; j++) {
        uint64_t g = gids (P)[hm_vertex (M, k, s, j)];
        put (&send[i], &g, sizeof (g));
      }
      put_fields (P, k, s, &send[i]);
    }
  }
  for (int k = 0; k <= m; k++)
    free (old[k]);
  exchange (np, peers, send, recv, 200);
  for (int i = 0; i < np; i++) {
    while (recv[i].at < recv[i].n) {
      int k;
      uint32_t v[HM_LEVELS];
      get (&recv[i], &k, sizeof (int));
      int ok = 1;
      for (int j = 0; j <= k; j++) {
        uint64_t g;
        get (&recv[i], &g, sizeof (g));
        v[j] = hm_part_local (P, g);
        ok &= v[j] != HM_NONE;
      }
      uint32_t s = ok ? hm_find (M, k, v) : HM_NONE;
      if (s == HM_NONE) {
        bad++;
        break;                         /* the rest cannot be read */
      }
      get_fields (P, k, s, &recv[i]);
    }
    free (send[i].p), free (recv[i].p);
    memset (&send[i], 0, sizeof (Buf)), memset (&recv[i], 0, sizeof (Buf));
  }

  free (send), free (recv);

  /* S5 and S6: the refresh, owned first; then c from the mesh as it is
     now */
  bad += refresh (P);
  save_K (P);
  int longer = hm_advect_end (A, c);
  (void) longer;
  long own_long = 0;
  if (m >= 1) {
    const int32_t * o = owners (P, 1);
    hm_foreach (M, 1, e) {
      if (o[e] != P->rank)
        continue;
      double a[3], b[3];
      point_of (M, hm_vertex (M, 1, e, 0), a);
      point_of (M, hm_vertex (M, 1, e, 1), b);
      own_long += !hm_edge_fits (M, a, b);
    }
  }
  if (reduce_sum (bad)) {
    if (bad)
      fprintf (stderr, "hmesh-mpi: rank %d: %d messages not applied\n",
               P->rank, bad);
    return -3;
  }
  return (int) reduce_sum (own_long);
}

#endif
