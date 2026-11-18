/*
 * mpi_topology.c — discover the topology of an MPI job on a cluster/supercomputer.
 *
 * What it reports
 *   1. Per rank : host, node id, local rank, current CPU, socket (package),
 *         NUMA domain, core, L3 domain, CPU affinity mask.
 *   2. Per node : sockets, NUMA domains, cores, hardware threads, memory,
 *         GPUs (from PCI), high-speed NICs (InfiniBand / Slingshot).
 *   3. Node 0   : NUMA distance matrix (from the kernel/ACPI SLIT table).
 *   4. Network  : node-to-node MPI ping-pong latency matrix + histogram.
 *         Latency tiers reveal the fabric structure (same switch,
 *         same group/rack, across the spine/global links).
 *
 * Build
 *   mpicc -O2 -o mpi_topology mpi_topology.c    # sysfs only (Linux)
 *   mpicc -O2 -DHAVE_HWLOC -o mpi_topology mpi_topology.c $ (pkg-config --cflags --libs hwloc)  # with hwloc >= 2.0
 *
 * Run (examples)
 *   srun -N 8 --ntasks-per-node=4 --cpu-bind=cores ./mpi_topology
 *   mpirun -np 64 --map-by ppr:8:node --bind-to core ./mpi_topology --iters 2000
 *
 * Options
 *   --iters N    ping-pong iterations per node pair (default 1000)
 *   --no-latency   skip the network latency measurement
 *   --out PREFIX   CSV file prefix (default "topo") -> PREFIX_ranks.csv,
 *          PREFIX_nodes.csv, PREFIX_latency.csv
 */
#define _GNU_SOURCE
#include <mpi.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <glob.h>
#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

#define HOSTLEN   64
#define STRLEN   160
#define MAX_NUMA  64

typedef struct {
  int  rank, node_id, local_rank, local_size;
  int  cpu, socket, numa, core, l3;
  char host[HOSTLEN];
  char affinity[STRLEN];
} rank_info_t;

typedef struct {
  int  node_id, nranks;
  int  sockets, numa, cores, pus, l3s, gpus, nnics;
  double mem_gb;
  char   host[HOSTLEN];
  char   nics[STRLEN];
} node_info_t;

/* ---------------------------------------------------------------- helpers */

static int read_int_file (const char *path, int dflt)
{
  FILE *f = fopen (path, "r");
  int v;
  if (!f) return dflt;
  if (fscanf (f, "%i", &v) != 1) v = dflt;   /* %i accepts 0x... too */
  fclose (f);
  return v;
}

static void append (char *dst, size_t len, const char *src)
{
  size_t used = strlen (dst);
  if (used + 1 < len) snprintf (dst + used, len - used, "%s", src);
}

/* Compress a cpu_set_t into "0-3,8,10-11". */
static void cpuset_to_string (const cpu_set_t *set, char *buf, size_t len)
{
  char tmp[32];
  int start = -1;
  buf[0] = '\0';
  for (int c = 0; c <= CPU_SETSIZE; c++) {
    int in = (c < CPU_SETSIZE) && CPU_ISSET (c, set);
    if (in && start < 0) start = c;
    if (!in && start >= 0) {
      if (buf[0]) append (buf, len, ",");
      if (c - 1 == start) snprintf (tmp, sizeof tmp, "%d", start);
      else        snprintf (tmp, sizeof tmp, "%d-%d", start, c - 1);
      append (buf, len, tmp);
      start = -1;
    }
  }
}

/* NUMA node of a CPU from sysfs: /sys/devices/system/cpu/cpuN/nodeK */
static int sysfs_cpu_numa (int cpu)
{
  char pat[96];
  glob_t g;
  int node = -1;
  snprintf (pat, sizeof pat, "/sys/devices/system/cpu/cpu%d/node[0-9]*", cpu);
  if (glob (pat, 0, NULL, &g) == 0 && g.gl_pathc > 0) {
    const char *p = strrchr (g.gl_pathv[0], 'e');   /* ".../nodeK" */
    if (p) node = atoi (p + 1);
  }
  globfree (&g);
  return node;
}

static int sysfs_count (const char *pattern)
{
  glob_t g;
  int n = 0;
  if (glob (pattern, 0, NULL, &g) == 0) n = (int)g.gl_pathc;
  globfree (&g);
  return n;
}

/* GPUs = PCI display/3D controllers from NVIDIA, AMD or Intel
 * (skips BMC VGA chips like ASPEED/Matrox). Works with or without hwloc. */
static int sysfs_count_gpus (void)
{
  DIR *d = opendir ("/sys/bus/pci/devices");
  struct dirent *e;
  char path[512];
  int n = 0;
  if (!d) return -1;
  while ( (e = readdir (d))) {
    if (e->d_name[0] == '.') continue;
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/class", e->d_name);
    int cls = read_int_file (path, 0);
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/vendor", e->d_name);
    int ven = read_int_file (path, 0);
    if ( (cls >> 16) == 0x03 && (ven == 0x10de || ven == 0x1002 || ven == 0x8086))
      n++;
  }
  closedir (d);
  return n;
}

/* High-speed NICs: RDMA devices (mlx5_0, cxi0, hfi1_0, ...) and
 * Slingshot/IB netdevs (hsn0, ib0). */
static int sysfs_list_nics (char *buf, size_t len)
{
  const char *dirs[] = { "/sys/class/infiniband", "/sys/class/cxi", "/sys/class/net" };
  int n = 0;
  buf[0] = '\0';
  for (size_t i = 0; i < sizeof dirs / sizeof *dirs; i++) {
    DIR *d = opendir (dirs[i]);
    struct dirent *e;
    if (!d) continue;
    while ( (e = readdir (d))) {
      if (e->d_name[0] == '.') continue;
      if (i == 2 && strncmp (e->d_name, "hsn", 3) && strncmp (e->d_name, "ib", 2))
        continue;   /* /sys/class/net: only hsnX / ibX */
      if (buf[0]) append (buf, len, " ");
      append (buf, len, e->d_name);
      n++;
    }
    closedir (d);
  }
  return n;
}

/* ------------------------------------------------- per-rank / per-node */

static void fill_rank_info (rank_info_t *ri, void *topo_ptr)
{
  cpu_set_t set;
  CPU_ZERO (&set);
  if (sched_getaffinity (0, sizeof set, &set) == 0)
    cpuset_to_string (&set, ri->affinity, sizeof ri->affinity);
  else
    snprintf (ri->affinity, sizeof ri->affinity, "?");

  ri->cpu = sched_getcpu ();
  ri->socket = ri->numa = ri->core = ri->l3 = -1;
  if (ri->cpu < 0) return;

#ifdef HAVE_HWLOC
  hwloc_topology_t topo = * (hwloc_topology_t *)topo_ptr;
  hwloc_obj_t pu = hwloc_get_pu_obj_by_os_index (topo, (unsigned)ri->cpu);
  if (pu) {
    hwloc_obj_t o;
    if ( (o = hwloc_get_ancestor_obj_by_type (topo, HWLOC_OBJ_PACKAGE, pu)))  ri->socket = (int)o->os_index;
    if ( (o = hwloc_get_ancestor_obj_by_type (topo, HWLOC_OBJ_CORE, pu)))   ri->core   = (int)o->logical_index;
    if ( (o = hwloc_get_ancestor_obj_by_type (topo, HWLOC_OBJ_L3CACHE, pu)))  ri->l3   = (int)o->logical_index;
    if (pu->nodeset) ri->numa = hwloc_bitmap_first (pu->nodeset);
    return;
  }
#else
  (void)topo_ptr;
#endif
  char path[128];
  snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", ri->cpu);
  ri->socket = read_int_file (path, -1);
  snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/topology/core_id", ri->cpu);
  ri->core = read_int_file (path, -1);      /* unique only within a socket */
  snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/cache/index3/id", ri->cpu);
  ri->l3 = read_int_file (path, -1);
  ri->numa = sysfs_cpu_numa (ri->cpu);
}

static void fill_node_info (node_info_t *ni, void *topo_ptr)
{
#ifdef HAVE_HWLOC
  hwloc_topology_t topo = * (hwloc_topology_t *)topo_ptr;
  ni->sockets = hwloc_get_nbobjs_by_type (topo, HWLOC_OBJ_PACKAGE);
  ni->numa  = hwloc_get_nbobjs_by_type (topo, HWLOC_OBJ_NUMANODE);
  ni->cores   = hwloc_get_nbobjs_by_type (topo, HWLOC_OBJ_CORE);
  ni->pus   = hwloc_get_nbobjs_by_type (topo, HWLOC_OBJ_PU);
  ni->l3s   = hwloc_get_nbobjs_by_type (topo, HWLOC_OBJ_L3CACHE);
  ni->mem_gb  = (double)hwloc_get_root_obj (topo)->total_memory / (1024.0 * 1024 * 1024);
#else
  (void)topo_ptr;
  /* Count distinct packages and (package, core) pairs from sysfs. */
  int ncpu = (int)sysconf (_SC_NPROCESSORS_CONF);
  int *pk = malloc (sizeof (int) * ncpu), *co = malloc (sizeof (int) * ncpu);
  int npk = 0, nco = 0;
  char path[128];
  for (int c = 0; c < ncpu; c++) {
    snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/online", c);
    if (read_int_file (path, 1) == 0) continue;      /* offline */
    snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", c);
    int p = read_int_file (path, -1);
    snprintf (path, sizeof path, "/sys/devices/system/cpu/cpu%d/topology/core_id", c);
    int k = read_int_file (path, -1);
    if (p < 0) continue;
    int seen = 0;
    for (int i = 0; i < npk; i++) if (pk[i] == p) { seen = 1; break; }
    if (!seen) pk[npk++] = p;
    seen = 0;
    for (int i = 0; i < nco; i++) if (co[i] == (p << 16 | k)) { seen = 1; break; }
    if (!seen) co[nco++] = p << 16 | k;
  }
  free (pk); free (co);
  ni->sockets = npk;
  ni->cores   = nco;
  ni->pus   = (int)sysconf (_SC_NPROCESSORS_ONLN);
  ni->numa  = sysfs_count ("/sys/devices/system/node/node[0-9]*");
  ni->l3s   = -1;
  ni->mem_gb  = (double)sysconf (_SC_PHYS_PAGES) * sysconf (_SC_PAGESIZE) / (1024.0 * 1024 * 1024);
#endif
  ni->gpus  = sysfs_count_gpus ();
  ni->nnics = sysfs_list_nics (ni->nics, sizeof ni->nics);
}

static void print_numa_distances (void)
{
  char path[96], line[1024];
  int n = sysfs_count ("/sys/devices/system/node/node[0-9]*");
  if (n <= 0) return;
  printf ("\n=== NUMA distance matrix on node 0 (10 = local) ===\n");
  for (int k = 0, shown = 0; shown < n && k < 4 * MAX_NUMA; k++) {
    snprintf (path, sizeof path, "/sys/devices/system/node/node%d/distance", k);
    FILE *f = fopen (path, "r");
    if (!f) continue;            /* node ids can have gaps */
    if (fgets (line, sizeof line, f)) {
      line[strcspn (line, "\n")] = '\0';
      printf ("  numa%-3d: %s\n", k, line);
    }
    fclose (f);
    shown++;
  }
}

/* ----------------------------------------------------- network latency */

/* Round-robin tournament: n/2 disjoint pairs ping-pong at once,
 * n-1 rounds in total, so the cost grows with n rather than n^2.
 * Result (one-way latency in microseconds) is reduced to leader 0. */
static void measure_latency (MPI_Comm leaders, int iters, double *mat)
{
  int me, n;
  MPI_Comm_rank (leaders, &me);
  MPI_Comm_size (leaders, &n);
  int m = (n % 2) ? n + 1 : n;         /* add a dummy if odd */
  double *local = calloc ( (size_t)n * n, sizeof (double));
  char buf[8] = {0};
  int warm = iters / 10 + 10;

  for (int r = 0; r < m - 1; r++) {
    int partner;
    if (me == m - 1)    partner = r;
    else if (me == r)   partner = m - 1;
    else          partner = ( (2 * r - me) % (m - 1) + (m - 1)) % (m - 1);

    if (partner < n) {
      int first = me < partner;
      double t0 = 0.0;
      for (int k = 0; k < warm + iters; k++) {
        if (k == warm) t0 = MPI_Wtime ();
        if (first) {
          MPI_Send (buf, 8, MPI_CHAR, partner, 0, leaders);
          MPI_Recv (buf, 8, MPI_CHAR, partner, 0, leaders, MPI_STATUS_IGNORE);
        } else {
          MPI_Recv (buf, 8, MPI_CHAR, partner, 0, leaders, MPI_STATUS_IGNORE);
          MPI_Send (buf, 8, MPI_CHAR, partner, 0, leaders);
        }
      }
      if (first) {
        double us = (MPI_Wtime () - t0) / iters / 2.0 * 1e6;
        local[me * n + partner] = local[partner * n + me] = us;
      }
    }
    MPI_Barrier (leaders);
  }
  MPI_Reduce (local, mat, n * n, MPI_DOUBLE, MPI_MAX, 0, leaders);
  free (local);
}

static int cmp_double (const void *a, const void *b)
{
  double x = * (const double *)a, y = * (const double *)b;
  return (x > y) - (x < y);
}

static void report_latency (const double *mat, const node_info_t *nodes, int n, const char *prefix)
{
  int np = n * (n - 1) / 2;
  if (np == 0) { printf ("\n (only one node: no inter-node latency to measure)\n"); return; }

  double *v = malloc (sizeof (double) * np);
  for (int i = 0, k = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) v[k++] = mat[i * n + j];
  qsort (v, np, sizeof (double), cmp_double);
  double lo = v[0], hi = v[np - 1];

  printf ("\n=== Inter-node one-way latency, 8-byte messages (us) ===\n");
  printf ("  pairs=%d  min=%.2f  median=%.2f  max=%.2f\n", np, lo, v[np / 2], hi);

  /* Histogram: separate peaks usually match fabric levels. */
  int bins = 10, hist[10] = {0};
  double w = (hi - lo) / bins;
  for (int k = 0; k < np; k++) {
    int b = (w > 0) ? (int) ( (v[k] - lo) / w) : 0;
    if (b >= bins) b = bins - 1;
    hist[b]++;
  }
  printf ("  histogram (separate peaks ~ switch / group / global hops):\n");
  for (int b = 0; b < bins; b++) {
    int bar = (int) (50.0 * hist[b] / np + 0.5);
    printf ("   %7.2f-%7.2f | %5d ", lo + b * w, lo + (b + 1) * w, hist[b]);
    for (int k = 0; k < bar; k++) putchar ('#');
    putchar ('\n');
  }
  free (v);

  if (n <= 24) {
    printf ("\n  matrix (rows/cols = node id):\n    ");
    for (int j = 0; j < n; j++) printf ("%7d", j);
    putchar ('\n');
    for (int i = 0; i < n; i++) {
      printf ("  %4d  ", i);
      for (int j = 0; j < n; j++)
        if (i == j) printf ("%7s", "-");
        else    printf ("%7.2f", mat[i * n + j]);
      printf ("   %s\n", nodes[i].host);
    }
  } else {
    printf ("  (%d nodes: full matrix written to CSV only)\n", n);
  }

  char fn[256];
  snprintf (fn, sizeof fn, "%s_latency.csv", prefix);
  FILE *f = fopen (fn, "w");
  if (f) {
    fprintf (f, "node");
    for (int j = 0; j < n; j++) fprintf (f, ",%s", nodes[j].host);
    fputc ('\n', f);
    for (int i = 0; i < n; i++) {
      fprintf (f, "%s", nodes[i].host);
      for (int j = 0; j < n; j++) fprintf (f, ",%.3f", mat[i * n + j]);
      fputc ('\n', f);
    }
    fclose (f);
    printf ("  -> %s\n", fn);
  }
}

/* ------------------------------------------------------------------ main */

int main (int argc, char **argv)
{
  int iters = 1000, do_latency = 1;
  const char *prefix = "topo";

  MPI_Init (&argc, &argv);
  for (int i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "--iters") && i + 1 < argc)  iters = atoi (argv[++i]);
    else if (!strcmp (argv[i], "--no-latency"))      do_latency = 0;
    else if (!strcmp (argv[i], "--out") && i + 1 < argc) prefix = argv[++i];
  }
  if (iters < 1) iters = 1;

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  /* Ranks sharing memory = same node. */
  MPI_Comm nodecomm, leaders;
  MPI_Comm_split_type (MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &nodecomm);
  int lrank, lsize;
  MPI_Comm_rank (nodecomm, &lrank);
  MPI_Comm_size (nodecomm, &lsize);

  /* One leader per node; node id = rank among leaders. */
  MPI_Comm_split (MPI_COMM_WORLD, lrank == 0 ? 0 : MPI_UNDEFINED, rank, &leaders);
  int node_id = 0, nnodes = 0;
  if (lrank == 0) {
    MPI_Comm_rank (leaders, &node_id);
    MPI_Comm_size (leaders, &nnodes);
  }
  MPI_Bcast (&node_id, 1, MPI_INT, 0, nodecomm);
  MPI_Bcast (&nnodes,  1, MPI_INT, 0, MPI_COMM_WORLD);

  void *topo_ptr = NULL;
#ifdef HAVE_HWLOC
  hwloc_topology_t topo;
  hwloc_topology_init (&topo);
  hwloc_topology_load (topo);
  topo_ptr = &topo;
#endif

  rank_info_t me;
  memset (&me, 0, sizeof me);
  me.rank = rank;
  me.node_id = node_id;
  me.local_rank = lrank;
  me.local_size = lsize;
  int hl;
  char pname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name (pname, &hl);
  snprintf (me.host, sizeof me.host, "%.*s", HOSTLEN - 1, pname);
  fill_rank_info (&me, topo_ptr);

  rank_info_t *all = NULL;
  if (rank == 0) all = malloc (sizeof (rank_info_t) * size);
  MPI_Gather (&me, sizeof me, MPI_BYTE, all, sizeof me, MPI_BYTE, 0, MPI_COMM_WORLD);

  node_info_t *nodes = NULL;
  double *lat = NULL;
  if (lrank == 0) {
    node_info_t ni;
    memset (&ni, 0, sizeof ni);
    ni.node_id = node_id;
    ni.nranks  = lsize;
    snprintf (ni.host, sizeof ni.host, "%s", me.host);
    fill_node_info (&ni, topo_ptr);
    if (rank == 0) nodes = malloc (sizeof (node_info_t) * nnodes);
    MPI_Gather (&ni, sizeof ni, MPI_BYTE, nodes, sizeof ni, MPI_BYTE, 0, leaders);
  }

  /* ---- report: ranks ---- */
  if (rank == 0) {
    printf ("=== MPI job: %d ranks on %d nodes (%s) ===\n", size, nnodes,
#ifdef HAVE_HWLOC
         "hwloc"
#else
         "sysfs"
#endif
    );
    int show = size <= 256 ? size : 64;
    printf ("\n%5s %5s %-20s %5s %5s %6s %5s %5s %4s  %s\n",
         "rank", "node", "host", "lrank", "cpu", "socket", "numa", "core", "L3", "affinity");
    for (int i = 0; i < show; i++) {
      rank_info_t *r = &all[i];
      printf ("%5d %5d %-20.20s %5d %5d %6d %5d %5d %4d  %s\n", r->rank, r->node_id, r->host,
           r->local_rank, r->cpu, r->socket, r->numa, r->core, r->l3, r->affinity);
    }
    if (show < size) printf ("  ... (%d more ranks in CSV)\n", size - show);

    /* Warn about ranks sharing a core: usually a binding mistake. */
    int overlaps = 0;
    for (int i = 0; i < size; i++)
      for (int j = i + 1; j < size; j++)
        if (all[i].node_id == all[j].node_id && all[i].core >= 0 &&
          all[i].core == all[j].core && all[i].socket == all[j].socket)
          overlaps++;
    if (overlaps)
      printf ("\n  WARNING: %d rank pairs currently run on the same core "
           " (check --cpu-bind / --bind-to).\n", overlaps);

    char fn[256];
    snprintf (fn, sizeof fn, "%s_ranks.csv", prefix);
    FILE *f = fopen (fn, "w");
    if (f) {
      fprintf (f, "rank,node,host,local_rank,cpu,socket,numa,core,l3,affinity\n");
      for (int i = 0; i < size; i++) {
        rank_info_t *r = &all[i];
        fprintf (f, "%d,%d,%s,%d,%d,%d,%d,%d,%d,\"%s\"\n", r->rank, r->node_id, r->host,
            r->local_rank, r->cpu, r->socket, r->numa, r->core, r->l3, r->affinity);
      }
      fclose (f);
      printf ("\n  -> %s\n", fn);
    }

    /* ---- report: nodes ---- */
    printf ("\n=== Nodes ===\n%5s %-20s %6s %7s %5s %6s %5s %4s %8s %5s  %s\n",
         "node", "host", "ranks", "sockets", "numa", "cores", "pus", "L3", "mem_GB", "gpus", "nics");
    for (int i = 0; i < nnodes; i++) {
      node_info_t *n = &nodes[i];
      printf ("%5d %-20.20s %6d %7d %5d %6d %5d %4d %8.1f %5d  %s\n", n->node_id, n->host,
           n->nranks, n->sockets, n->numa, n->cores, n->pus, n->l3s, n->mem_gb, n->gpus,
           n->nnics ? n->nics : "-");
    }
    snprintf (fn, sizeof fn, "%s_nodes.csv", prefix);
    f = fopen (fn, "w");
    if (f) {
      fprintf (f, "node,host,ranks,sockets,numa,cores,pus,l3,mem_gb,gpus,nics\n");
      for (int i = 0; i < nnodes; i++) {
        node_info_t *n = &nodes[i];
        fprintf (f, "%d,%s,%d,%d,%d,%d,%d,%d,%.2f,%d,\"%s\"\n", n->node_id, n->host, n->nranks,
            n->sockets, n->numa, n->cores, n->pus, n->l3s, n->mem_gb, n->gpus, n->nics);
      }
      fclose (f);
      printf ("  -> %s\n", fn);
    }

    print_numa_distances ();
    fflush (stdout);
  }

  /* ---- network ---- */
  if (do_latency && lrank == 0) {
    if (rank == 0) lat = calloc ( (size_t)nnodes * nnodes, sizeof (double));
    measure_latency (leaders, iters, lat);
    if (rank == 0) report_latency (lat, nodes, nnodes, prefix);
  }

  if (leaders != MPI_COMM_NULL) MPI_Comm_free (&leaders);
  MPI_Comm_free (&nodecomm);
#ifdef HAVE_HWLOC
  hwloc_topology_destroy (topo);
#endif
  free (all); free (nodes); free (lat);
  MPI_Finalize ();
  return 0;
}
