/* mpi-topology.h: the part of a machine an MPI job was given. Without
.. _MPI, nothing. */

/* cpu_set_t, CPU_ISSET () and sched_getcpu () are GNU's: the macro must
.. come before the first system header, which fixes the features */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#ifdef _MPI

#include <dirent.h>
#include <fcntl.h>
#include <glob.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef HAVE_HWLOC
# include <hwloc.h>
#endif

#include "mpi-topology.h"

#define SYS_CPU   "/sys/devices/system/cpu"
#define SYS_NODE  "/sys/devices/system/node"
#define GIB       (1024. * 1024. * 1024.)

/*
.. Sets of CPUs and of NUMA domains: bit masks of NBITS bits in unsigned
.. longs, which MPI_BOR can join (the union over a node's ranks). 1024 as
.. glibc's CPU_SETSIZE: no node has more CPUs or NUMA domains yet.
*/
#define NBITS     1024
#define WORD      (8 * (int) sizeof (unsigned long))
#define NWORDS    (NBITS / WORD)

typedef struct { unsigned long w [NWORDS]; } Mask;

static int mask_has (const Mask * m, int i)
{
  return i >= 0 && i < NBITS && (m->w [i / WORD] >> (i % WORD)) & 1;
}

static void mask_add (Mask * m, int i)
{
  if (i >= 0 && i < NBITS)
    m->w [i / WORD] |= 1UL << (i % WORD);
}

static int mask_count (const Mask * m)
{
  int n = 0;
  for (int i = 0; i < NBITS; i++)
    n += mask_has (m, i);
  return n;
}

/* the first member, or -1 */
static int mask_first (const Mask * m)
{
  for (int i = 0; i < NBITS; i++)
    if (mask_has (m, i))
      return i;
  return -1;
}

/* ---------------------------------------------------------------------- */
/* text and files                                                         */
/* ---------------------------------------------------------------------- */

/* zeroed memory, or the job aborted */
static void * alloc (size_t n)
{
  void * p = calloc (1, n ? n : 1);
  if (!p)
  {
    fprintf (stderr, "mpi-topology: out of memory\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  return p;
}

/* src after the text in dst, cut at len */
static void append (char * dst, size_t len, const char * src)
{
  size_t used = strlen (dst);
  if (used + 1 < len)
    snprintf (dst + used, len - used, "%s", src);
}

/* a mask as a list "0-3,8,10-11", cut at len ("-" if empty) */
static void mask_list (const Mask * m, char * buf, size_t len)
{
  char run [32];
  int start = -1;
  buf [0] = '\0';
  /* one past the last bit, so that a run reaching it is closed */
  for (int i = 0; i <= NBITS; i++)
  {
    int in = mask_has (m, i);
    if (in && start < 0)
      start = i;
    else if (!in && start >= 0)
    {
      if (i - 1 == start)
        snprintf (run, sizeof run, "%s%d", buf [0] ? "," : "", start);
      else
        snprintf (run, sizeof run, "%s%d-%d", buf [0] ? "," : "", start,
                  i - 1);
      append (buf, len, run);
      start = -1;
    }
  }
  if (!buf [0])
    snprintf (buf, len, "-");
}

/* a list "0-3,8,10-11" (the kernel's form) into m; its members */
static int mask_parse (const char * s, Mask * m)
{
  memset (m, 0, sizeof * m);
  while (*s && *s != '\n')
  {
    char * end;
    long lo = strtol (s, &end, 10), hi = lo;
    if (end == s)
      break;
    s = end;
    if (*s == '-')
    {
      hi = strtol (s + 1, &end, 10);
      s = end;
    }
    for (long i = lo; i <= hi && i < NBITS; i++)
      mask_add (m, (int) i);
    if (*s == ',')
      s++;
  }
  return mask_count (m);
}

/* the first line of a file, its newline dropped; 0 if none */
static int read_line (const char * path, char * buf, size_t len)
{
  FILE * f = fopen (path, "r");
  if (!f)
    return 0;
  int ok = fgets (buf, (int) len, f) != NULL;
  fclose (f);
  if (ok)
    buf [strcspn (buf, "\n")] = '\0';
  return ok;
}

/*
.. the number in a file, or dflt: hexadecimal with 0x (PCI class and
.. vendor), else decimal. Not strtol's base 0, which would read a leading
.. 0 as octal
*/
static long read_number (const char * path, long dflt)
{
  char line [64];
  if (!read_line (path, line, sizeof line))
    return dflt;
  char * end;
  int hex = line [0] == '0' && (line [1] == 'x' || line [1] == 'X');
  long x = strtol (line, &end, hex ? 16 : 10);
  return end != line ? x : dflt;
}

/* the number in SYS_CPU/cpuN/what, or -1 */
static int cpu_number (int cpu, const char * what)
{
  char path [160];
  snprintf (path, sizeof path, SYS_CPU "/cpu%d/%s", cpu, what);
  return (int) read_number (path, -1);
}

#ifndef HAVE_HWLOC
/* the paths matching a pattern (mpi_topology_node () without hwloc) */
static int count_paths (const char * pattern)
{
  glob_t g;
  int n = 0;
  if (glob (pattern, 0, NULL, &g) == 0)
    n = (int) g.gl_pathc;
  globfree (&g);
  return n;
}
#endif

/* the NUMA domain of a CPU: sysfs links SYS_CPU/cpuN/nodeK, or -1 */
static int cpu_numa (int cpu)
{
  char pattern [96];
  glob_t g;
  int node = -1;
  if (cpu < 0)
    return -1;
  snprintf (pattern, sizeof pattern, SYS_CPU "/cpu%d/node[0-9]*", cpu);
  if (glob (pattern, 0, NULL, &g) == 0 && g.gl_pathc > 0)
    node = atoi (strrchr (g.gl_pathv [0], '/') + 1 + 4);    /* "nodeK" */
  globfree (&g);
  return node;
}

/* the memory nodes this process may allocate on (its cpuset cgroup), as
.. the kernel lists them in /proc/self/status; "?" if not found */
static void mems_allowed (char * buf, size_t len)
{
  FILE * f = fopen ("/proc/self/status", "r");
  char line [512];
  snprintf (buf, len, "?");
  if (!f)
    return;
  while (fgets (line, sizeof line, f))
    if (!strncmp (line, "Mems_allowed_list:", 18))
    {
      const char * s = line + 18;
      while (*s == ' ' || *s == '\t')
        s++;
      int n = (int) strcspn (s, "\n");
      snprintf (buf, len, "%.*s", n < (int) len ? n : (int) len - 1, s);
      break;
    }
  fclose (f);
}

void mpi_topology_cpuset (const void * set, char * buf, size_t len)
{
  const cpu_set_t * s = set;
  Mask m = { { 0 } };
  for (int c = 0; c < CPU_SETSIZE && c < NBITS; c++)
    if (CPU_ISSET (c, s))
      mask_add (&m, c);
  mask_list (&m, buf, len);
}

/* the CPUs this process may run on */
static Mask affinity_mask (void)
{
  Mask m = { { 0 } };
  cpu_set_t set;
  CPU_ZERO (&set);
  if (sched_getaffinity (0, sizeof set, &set) == 0)
    for (int c = 0; c < CPU_SETSIZE && c < NBITS; c++)
      if (CPU_ISSET (c, &set))
        mask_add (&m, c);
  return m;
}

/* ---------------------------------------------------------------------- */
/* GPUs                                                                   */
/* ---------------------------------------------------------------------- */

/* a GPU of the node, from sysfs and the driver */
typedef struct
{
  char pci [16];                      /* 0000:01:00.0                   */
  int vendor, device, numa;
  char model [64];
  char uuid [48];                     /* NVIDIA's GPU-..., else empty   */
  char file [32];                     /* the device file a job opens    */
  char local [MPI_TOPOLOGY_TEXT];     /* the CPUs close to it           */
} Gpu;

static int pci_cmp (const void * a, const void * b)
{
  return strcmp (((const Gpu *) a)->pci, ((const Gpu *) b)->pci);
}

/* the value after "key:" in the lines of a file, spaces trimmed */
static int read_field (const char * path, const char * key, char * buf,
                       size_t len)
{
  FILE * f = fopen (path, "r");
  char line [256];
  size_t k = strlen (key);
  int found = 0;
  if (!f)
    return 0;
  while (!found && fgets (line, sizeof line, f))
    if (!strncmp (line, key, k) && line [k] == ':')
    {
      const char * s = line + k + 1;
      while (*s == ' ' || *s == '\t')
        s++;
      int n = (int) strcspn (s, "\n");
      snprintf (buf, len, "%.*s", n < (int) len ? n : (int) len - 1, s);
      found = 1;
    }
  fclose (f);
  return found;
}

/*
.. is a PCI device of this class and vendor a GPU? A display controller
.. (class 0x03xxxx) or a processing accelerator (0x12xxxx: AMD's MI300X
.. lists as one) of NVIDIA, AMD or Intel; not a server's management VGA
.. chip (ASPEED, Matrox), nor another vendor's accelerator (Habana Gaudi,
.. Huawei Ascend), whose device files differ
*/
static int is_gpu (long cls, long vendor)
{
  return ((cls >> 16) == 0x03 || (cls >> 16) == 0x12)
         && (vendor == 0x10de || vendor == 0x1002 || vendor == 0x8086);
}

/*
.. The node's GPUs, at most max, in PCI order (is_gpu ()). NVIDIA's
.. driver names the model, the UUID and the minor of /dev/nvidiaN in
.. /proc/driver/nvidia/gpus/PCI/information; for the others the model is
.. AMD's product_name if there, else vendor:device, and the device file
.. the render node of the PCI device's drm directory.
*/
static int node_gpus (Gpu * g, int max)
{
  DIR * d = opendir ("/sys/bus/pci/devices");
  if (!d)
    return 0;
  struct dirent * e;
  char path [512], text [128];
  int n = 0;
  while ((e = readdir (d)) && n < max)
  {
    const char * pci = e->d_name;
    if (pci [0] == '.' || strlen (pci) >= sizeof g->pci)
      continue;
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/class", pci);
    long cls = read_number (path, 0);
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/vendor", pci);
    long vendor = read_number (path, 0);
    if (!is_gpu (cls, vendor))
      continue;
    Gpu * x = &g [n++];
    memset (x, 0, sizeof * x);
    snprintf (x->pci, sizeof x->pci, "%.*s", (int) sizeof x->pci - 1, pci);
    x->vendor = (int) vendor;
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/device", pci);
    x->device = (int) read_number (path, 0);
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/numa_node", pci);
    x->numa = (int) read_number (path, -1);    /* -1: one domain, or not
                                                  known to the firmware */
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/local_cpulist",
              pci);
    if (!read_line (path, x->local, sizeof x->local))
      snprintf (x->local, sizeof x->local, "-");
    snprintf (x->model, sizeof x->model, "%04x:%04x", x->vendor, x->device);
    snprintf (path, sizeof path, "/proc/driver/nvidia/gpus/%s/information",
              pci);
    if (vendor == 0x10de && read_field (path, "Model", text, sizeof text))
    {
      snprintf (x->model, sizeof x->model, "%.*s", (int) sizeof x->model - 1,
                text);
      if (read_field (path, "GPU UUID", text, sizeof text))
        snprintf (x->uuid, sizeof x->uuid, "%.*s", (int) sizeof x->uuid - 1,
                  text);
      if (read_field (path, "Device Minor", text, sizeof text))
        snprintf (x->file, sizeof x->file, "/dev/nvidia%d", atoi (text));
      continue;
    }
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/product_name",
              pci);
    if (read_line (path, text, sizeof text) && text [0])
      snprintf (x->model, sizeof x->model, "%.*s", (int) sizeof x->model - 1,
                text);
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/drm", pci);
    DIR * dd = opendir (path);
    struct dirent * de;
    while (dd && (de = readdir (dd)))
      if (!strncmp (de->d_name, "renderD", 7))
      {
        snprintf (x->file, sizeof x->file, "/dev/dri/%.20s", de->d_name);
        break;
      }
    if (dd)
      closedir (dd);
  }
  closedir (d);
  qsort (g, (size_t) n, sizeof (Gpu), pci_cmp);
  return n;
}

/* can this process open a device file? */
static int opens (const char * file)
{
  int fd = open (file, O_RDWR | O_CLOEXEC);
  if (fd < 0)
    return 0;
  close (fd);
  return 1;
}

/*
.. can this process compute on the GPU? Its device file opens, as the
.. cgroup allows. An AMD GPU needs /dev/kfd too, the one file of ROCm's
.. compute driver for all AMD GPUs: a job given the render node but not
.. /dev/kfd could draw, not compute, and the GPU is not counted
*/
static int gpu_opens (const Gpu * g)
{
  return g->file [0] && opens (g->file)
         && (g->vendor != 0x1002 || opens ("/dev/kfd"));
}

/* the *_VISIBLE_DEVICES set, as NAME=value;... ("-" if none) */
static void gpu_environment (char * buf, size_t len)
{
  const char * names [] = { "CUDA_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES",
                            "HIP_VISIBLE_DEVICES", "ZE_AFFINITY_MASK",
                            "GPU_DEVICE_ORDINAL", "SLURM_STEP_GPUS" };
  buf [0] = '\0';
  for (int i = 0; i < (int) (sizeof names / sizeof * names); i++)
  {
    const char * v = getenv (names [i]);
    if (!v)
      continue;
    if (buf [0])
      append (buf, len, ";");
    append (buf, len, names [i]);
    append (buf, len, "=");
    append (buf, len, v);
  }
  if (!buf [0])
    snprintf (buf, len, "-");
}

/* a set of GPU indices as a list "0-1,3" */
static void gpu_list (unsigned long bits, char * buf, size_t len)
{
  Mask m = { { 0 } };
  m.w [0] = bits;
  mask_list (&m, buf, len);
}

/* ---------------------------------------------------------------------- */
/* a rank                                                                 */
/* ---------------------------------------------------------------------- */

void mpi_topology_rank (MpiTopologyRank * r, void * topo)
{
  Mask a = affinity_mask ();
  mask_list (&a, r->affinity, sizeof r->affinity);
  mems_allowed (r->mems, sizeof r->mems);
  r->cpu = sched_getcpu ();
  /* the GPUs it can open, and those whose close CPUs hold its CPU */
  Gpu g [MPI_TOPOLOGY_GPUS];
  int ng = node_gpus (g, MPI_TOPOLOGY_GPUS);
  r->gpus = r->near = 0;
  for (int i = 0; i < ng; i++)
  {
    Mask local;
    if (gpu_opens (&g [i]))
      r->gpus |= 1UL << i;
    if (mask_parse (g [i].local, &local) && mask_has (&local, r->cpu))
      r->near |= 1UL << i;
  }
  gpu_environment (r->gpu_env, sizeof r->gpu_env);
  r->socket = r->numa = r->core = r->l3 = -1;
  if (r->cpu < 0)
    return;
#ifdef HAVE_HWLOC
  hwloc_topology_t t = * (hwloc_topology_t *) topo;
  hwloc_obj_t pu = hwloc_get_pu_obj_by_os_index (t, (unsigned) r->cpu);
  if (pu)
  {
    hwloc_obj_t o;
    if ((o = hwloc_get_ancestor_obj_by_type (t, HWLOC_OBJ_PACKAGE, pu)))
      r->socket = (int) o->os_index;
    if ((o = hwloc_get_ancestor_obj_by_type (t, HWLOC_OBJ_CORE, pu)))
      r->core = (int) o->logical_index;
    if ((o = hwloc_get_ancestor_obj_by_type (t, HWLOC_OBJ_L3CACHE, pu)))
      r->l3 = (int) o->logical_index;
    if (pu->nodeset)
      r->numa = hwloc_bitmap_first (pu->nodeset);
    return;
  }
  /* a CPU hwloc does not know (restricted topology): sysfs below */
#else
  (void) topo;
#endif
  r->socket = cpu_number (r->cpu, "topology/physical_package_id");
  r->core = cpu_number (r->cpu, "topology/core_id");
  r->l3 = cpu_number (r->cpu, "cache/index3/id");
  r->numa = cpu_numa (r->cpu);
}

/* ---------------------------------------------------------------------- */
/* a node's allocation and its caches                                     */
/* ---------------------------------------------------------------------- */

/* what a node gave the job: the unions over its ranks */
typedef struct
{
  int node, nranks, ncpus, ncores, ngpus, nodegpus;
  char gpus [MPI_TOPOLOGY_TEXT];       /* the GPUs (indices on the node) */
  char host [MPI_TOPOLOGY_HOST];
  char cpus [MPI_TOPOLOGY_TEXT];       /* the CPUs                       */
  char numa [MPI_TOPOLOGY_TEXT];       /* the NUMA domains of those CPUs */
  char mems [MPI_TOPOLOGY_TEXT];       /* the memory nodes allowed       */
} Allocation;

/*
.. A cache instance serving an allocated CPU: its level, type (d data, i
.. instruction, u unified), size, the CPUs sharing it and those of them
.. allocated, and the NUMA domain of the first of those
*/
typedef struct
{
  int node, numa, level, size_kb, nshared, nused;
  char type;
  char shared [MPI_TOPOLOGY_TEXT];
  char used [MPI_TOPOLOGY_TEXT];
} Cache;

/* a node's caches in short, to compare nodes: instances and kB a level */
enum { LEVELS = 5 };
typedef struct
{
  int count [LEVELS];
  long kb [LEVELS];
} CacheSummary;

static int int_cmp (const void * a, const void * b)
{
  int x = * (const int *) a, y = * (const int *) b;
  return (x > y) - (x < y);
}

/* the distinct values of v [0..n), sorted in place */
static int distinct_ints (int * v, int n)
{
  qsort (v, (size_t) n, sizeof (int), int_cmp);
  int d = 0;
  for (int i = 0; i < n; i++)
    d += !i || v [i] != v [i - 1];
  return d;
}

/* "48K", "1280K", "32M": kB, or -1 */
static int size_kb (const char * s)
{
  char * end;
  long v = strtol (s, &end, 10);
  if (end == s)
    return -1;
  return (int) (*end == 'M' ? v * 1024 : *end == 'G' ? v * 1024 * 1024 : v);
}

/* the allocation of this node from the unions of its CPUs, memory nodes
.. and GPUs */
static void allocation_of (const Mask * cpus, const Mask * mems,
                           unsigned long gpus, Allocation * a)
{
  Gpu g [MPI_TOPOLOGY_GPUS];
  a->nodegpus = node_gpus (g, MPI_TOPOLOGY_GPUS);
  a->ngpus = __builtin_popcountl (gpus);
  gpu_list (gpus, a->gpus, sizeof a->gpus);
  Mask numa = { { 0 } };
  int n = mask_count (cpus), k = 0;
  int * cores = alloc ((size_t) n * sizeof (int));
  for (int c = 0; c < NBITS; c++)
    if (mask_has (cpus, c))
    {
      mask_add (&numa, cpu_numa (c));
      /* core_id is unique within a package only: the pair, packed */
      int p = cpu_number (c, "topology/physical_package_id");
      int id = cpu_number (c, "topology/core_id");
      if (p >= 0 && id >= 0 && id <= 0xffff)
        cores [k++] = p << 16 | id;
    }
  a->ncpus = n;
  a->ncores = k ? distinct_ints (cores, k) : -1;
  free (cores);
  mask_list (cpus, a->cpus, sizeof a->cpus);
  mask_list (&numa, a->numa, sizeof a->numa);
  mask_list (mems, a->mems, sizeof a->mems);
}

/*
.. The cache instances serving the allocated CPUs, from sysfs:
.. SYS_CPU/cpuN/cache/indexK/{level,type,size,shared_cpu_list} for each
.. allocated CPU N. An instance is known by its level, type and the CPUs
.. sharing it, so that the CPUs it serves list it once. Returns the
.. instances, *n of them (none where sysfs has no cache directory).
*/
static Cache * caches_of (const Mask * cpus, int node, int * n)
{
  int cap = 64, k = 0;
  Cache * c = alloc ((size_t) cap * sizeof (Cache));
  char path [192], text [MPI_TOPOLOGY_TEXT];
  for (int cpu = 0; cpu < NBITS; cpu++)
  {
    if (!mask_has (cpus, cpu))
      continue;
    for (int index = 0; ; index++)
    {
      char base [128];
      snprintf (base, sizeof base, SYS_CPU "/cpu%d/cache/index%d", cpu,
                index);
      snprintf (path, sizeof path, "%s/level", base);
      int level = (int) read_number (path, -1);
      if (level < 0)
        break;                        /* no more caches for this CPU */
      Cache e = { .node = node, .level = level, .type = 'u' };
      snprintf (path, sizeof path, "%s/type", base);
      if (read_line (path, text, sizeof text))
        e.type = text [0] == 'D' ? 'd' : text [0] == 'I' ? 'i' : 'u';
      snprintf (path, sizeof path, "%s/size", base);
      e.size_kb = read_line (path, text, sizeof text) ? size_kb (text) : -1;
      snprintf (path, sizeof path, "%s/shared_cpu_list", base);
      if (!read_line (path, e.shared, sizeof e.shared))
        snprintf (e.shared, sizeof e.shared, "%d", cpu);
      int seen = 0;
      for (int i = 0; i < k && !seen; i++)
        seen = c [i].level == e.level && c [i].type == e.type
               && !strcmp (c [i].shared, e.shared);
      if (seen)
        continue;
      Mask s, used = { { 0 } };
      e.nshared = mask_parse (e.shared, &s);
      for (int w = 0; w < NWORDS; w++)
        used.w [w] = s.w [w] & cpus->w [w];
      e.nused = mask_count (&used);
      mask_list (&used, e.used, sizeof e.used);
      e.numa = cpu_numa (mask_first (&used));
      if (k == cap)
      {
        Cache * more = alloc ((size_t) 2 * cap * sizeof (Cache));
        memcpy (more, c, (size_t) k * sizeof (Cache));
        free (c);
        c = more, cap *= 2;
      }
      c [k++] = e;
    }
  }
  *n = k;
  return c;
}

static CacheSummary summary_of (const Cache * c, int n)
{
  CacheSummary s = { { 0 }, { 0 } };
  for (int i = 0; i < n; i++)
    if (c [i].level > 0 && c [i].level < LEVELS)
    {
      s.count [c [i].level]++;
      s.kb [c [i].level] += c [i].size_kb > 0 ? c [i].size_kb : 0;
    }
  return s;
}

/* ---------------------------------------------------------------------- */
/* the whole node (mpi_topology_node (), not used by mpi_topology ())     */
/* ---------------------------------------------------------------------- */

/*
.. the GPUs (is_gpu ()). An integrated Intel GPU counts too: a laptop with
.. a discrete GPU shows 2. Every GPU of the node, the job's or not. -1
.. without sysfs
*/
static int count_gpus (void)
{
  DIR * d = opendir ("/sys/bus/pci/devices");
  if (!d)
    return -1;
  struct dirent * e;
  char path [512];
  int n = 0;
  while ((e = readdir (d)))
  {
    if (e->d_name [0] == '.')
      continue;
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/class", e->d_name);
    long cls = read_number (path, 0);
    snprintf (path, sizeof path, "/sys/bus/pci/devices/%s/vendor", e->d_name);
    long vendor = read_number (path, 0);
    n += is_gpu (cls, vendor);
  }
  closedir (d);
  return n;
}

/*
.. the high-speed NICs, their names into buf: RDMA devices (InfiniBand's
.. mlx5_0, Omni-Path's hfi1_0), HPE Slingshot's cxi devices, and the
.. network interfaces of those fabrics (hsnN for Slingshot, ibN for IP
.. over InfiniBand), not Ethernet's
*/
static int list_nics (char * buf, size_t len)
{
  const char * dirs [] = { "/sys/class/infiniband", "/sys/class/cxi",
                           "/sys/class/net" };
  int n = 0;
  buf [0] = '\0';
  for (int i = 0; i < 3; i++)
  {
    DIR * d = opendir (dirs [i]);
    if (!d)
      continue;
    struct dirent * e;
    while ((e = readdir (d)))
    {
      const char * name = e->d_name;
      if (name [0] == '.'
          || (i == 2 && strncmp (name, "hsn", 3) && strncmp (name, "ib", 2)))
        continue;
      if (buf [0])
        append (buf, len, " ");
      append (buf, len, name);
      n++;
    }
    closedir (d);
  }
  return n;
}

void mpi_topology_node (MpiTopologyNode * n, void * topo)
{
#ifdef HAVE_HWLOC
  hwloc_topology_t t = * (hwloc_topology_t *) topo;
  n->sockets = hwloc_get_nbobjs_by_type (t, HWLOC_OBJ_PACKAGE);
  n->numa = hwloc_get_nbobjs_by_type (t, HWLOC_OBJ_NUMANODE);
  n->cores = hwloc_get_nbobjs_by_type (t, HWLOC_OBJ_CORE);
  n->pus = hwloc_get_nbobjs_by_type (t, HWLOC_OBJ_PU);
  n->l3s = hwloc_get_nbobjs_by_type (t, HWLOC_OBJ_L3CACHE);
  n->mem_gb = (double) hwloc_get_root_obj (t)->total_memory / GIB;
#else
  (void) topo;
  /* the online CPUs: the distinct packages, and the distinct (package,
   * core) pairs, core_id being unique within a package only */
  int ncpu = (int) sysconf (_SC_NPROCESSORS_CONF), np = 0, nc = 0;
  int * pkg = alloc ((size_t) ncpu * sizeof (int));
  int * core = alloc ((size_t) ncpu * sizeof (int));
  for (int c = 0; c < ncpu; c++)
  {
    if (cpu_number (c, "online") == 0)    /* cpu0 has no file: online */
      continue;
    int p = cpu_number (c, "topology/physical_package_id");
    int k = cpu_number (c, "topology/core_id");
    if (p < 0 || k < 0 || k > 0xffff)
      continue;
    pkg [np++] = p;
    core [nc++] = p << 16 | k;
  }
  n->sockets = distinct_ints (pkg, np);
  n->cores = distinct_ints (core, nc);
  free (pkg), free (core);
  n->pus = (int) sysconf (_SC_NPROCESSORS_ONLN);
  n->numa = count_paths (SYS_NODE "/node[0-9]*");
  n->l3s = -1;
  n->mem_gb = (double) sysconf (_SC_PHYS_PAGES) * sysconf (_SC_PAGESIZE) / GIB;
#endif
  n->gpus = count_gpus ();
  n->nnics = list_nics (n->nics, sizeof n->nics);
}

/* ---------------------------------------------------------------------- */
/* latency                                                                */
/* ---------------------------------------------------------------------- */

/*
.. Every pair of ranks, once, by the circle method of round-robin
.. tournaments: with m players (n, or n + 1 with a dummy if n is odd),
.. player m - 1 stays and the others turn. In round r, r plays m - 1 and
.. each other i plays (2 r - i) mod (m - 1); the n / 2 pairs of a round
.. are disjoint, so they run at once, and m - 1 rounds cover the n (n -
.. 1) / 2 pairs. A barrier after each round keeps the pairs of different
.. rounds apart. The pairs of a round do share the nodes' memory and
.. links: the latency measured is that of every rank talking at once, as
.. in a halo exchange, not of one pair on a quiet machine.
..
.. A pair warms up (iters / 10 + 10 round trips, not timed), then times
.. iters round trips of 8 bytes; the one-way latency is half a round
.. trip, in microseconds. The lower rank of the pair keeps it in its own
.. row; rank 0 gathers the rows (O(n) memory a rank, n^2 on rank 0) and
.. copies each entry across the diagonal. The diagonal stays 0.
*/
void mpi_topology_latency (MPI_Comm comm, int iters, double * mat)
{
  int me, n;
  MPI_Comm_rank (comm, &me);
  MPI_Comm_size (comm, &n);
  int m = n % 2 ? n + 1 : n, warm = iters / 10 + 10;
  double * row = alloc ((size_t) n * sizeof (double));
  char buf [8] = { 0 };
  for (int r = 0; r < m - 1; r++)
  {
    int partner = me == m - 1 ? r
                : me == r     ? m - 1
                : ((2 * r - me) % (m - 1) + (m - 1)) % (m - 1);
    if (partner < n)                  /* else this rank sits out (dummy) */
    {
      int first = me < partner;
      double t0 = 0;
      for (int k = 0; k < warm + iters; k++)
      {
        if (k == warm)
          t0 = MPI_Wtime ();
        if (first)
        {
          MPI_Send (buf, 8, MPI_CHAR, partner, 0, comm);
          MPI_Recv (buf, 8, MPI_CHAR, partner, 0, comm, MPI_STATUS_IGNORE);
        }
        else
        {
          MPI_Recv (buf, 8, MPI_CHAR, partner, 0, comm, MPI_STATUS_IGNORE);
          MPI_Send (buf, 8, MPI_CHAR, partner, 0, comm);
        }
      }
      if (first)
        row [partner] = (MPI_Wtime () - t0) / iters / 2 * 1e6;
    }
    MPI_Barrier (comm);
  }
  MPI_Gather (row, n, MPI_DOUBLE, mat, n, MPI_DOUBLE, 0, comm);
  free (row);
  if (me == 0)
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        mat [j * n + i] = mat [i * n + j];
}

/* the classes of a pair of ranks, nearest first */
enum { SAME_CORE, SAME_L3, SAME_NUMA, SAME_SOCKET, SAME_NODE, OTHER_NODE,
       CLASSES };
static const char * class_name [CLASSES] = {
  "same core (SMT)", "same L3", "same NUMA domain", "same socket",
  "same node", "another node"
};

/*
.. the class of a pair from their placements. Unknown fields (-1) go to
.. the farther class: a pair whose L3 is not known is "same NUMA domain"
*/
static int pair_class (const MpiTopologyRank * a, const MpiTopologyRank * b)
{
  if (a->node != b->node)
    return OTHER_NODE;
  if (a->socket < 0 || a->socket != b->socket)
    return SAME_NODE;
  if (a->numa < 0 || a->numa != b->numa)
    return SAME_SOCKET;
  if (a->core >= 0 && a->core == b->core)
    return SAME_CORE;
  if (a->l3 >= 0 && a->l3 == b->l3)
    return SAME_L3;
  return SAME_NUMA;
}

static int double_cmp (const void * a, const void * b)
{
  double x = * (const double *) a, y = * (const double *) b;
  return (x > y) - (x < y);
}

/* the latencies by class (pairs, min, median, max), the matrix if small,
.. and every pair in the CSV file */
static void report_latency (const double * mat, const MpiTopologyRank * all,
                            int n, const char * prefix)
{
  if (n < 2)
  {
    printf ("\n  one rank: no latency\n");
    return;
  }
  int np = n * (n - 1) / 2, count [CLASSES] = { 0 };
  double * v [CLASSES];
  for (int c = 0; c < CLASSES; c++)
    v [c] = alloc ((size_t) np * sizeof (double));
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
    {
      int c = pair_class (&all [i], &all [j]);
      v [c][count [c]++] = mat [i * n + j];
    }
  printf ("\n=== One-way latency between ranks, 8-byte messages (us) ===\n");
  printf ("  %-18s %8s %8s %8s %8s\n", "class", "pairs", "min", "median",
          "max");
  for (int c = 0; c < CLASSES; c++)
  {
    if (!count [c])
      continue;
    qsort (v [c], (size_t) count [c], sizeof (double), double_cmp);
    printf ("  %-18s %8d %8.2f %8.2f %8.2f\n", class_name [c], count [c],
            v [c][0], v [c][count [c] / 2], v [c][count [c] - 1]);
  }
  for (int c = 0; c < CLASSES; c++)
    free (v [c]);

  if (n <= 16)
  {
    printf ("\n  matrix (rows and columns: rank):\n      ");
    for (int j = 0; j < n; j++)
      printf ("%7d", j);
    putchar ('\n');
    for (int i = 0; i < n; i++)
    {
      printf ("  %4d", i);
      for (int j = 0; j < n; j++)
        if (i == j)
          printf ("%7s", "-");
        else
          printf ("%7.2f", mat [i * n + j]);
      putchar ('\n');
    }
  }

  char fn [256];
  snprintf (fn, sizeof fn, "%s_latency.csv", prefix);
  FILE * f = fopen (fn, "w");
  if (!f)
    return;
  fprintf (f, "rank_a,rank_b,node_a,node_b,cpu_a,cpu_b,class,us\n");
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      fprintf (f, "%d,%d,%d,%d,%d,%d,\"%s\",%.3f\n", i, j, all [i].node,
               all [j].node, all [i].cpu, all [j].cpu,
               class_name [pair_class (&all [i], &all [j])], mat [i * n + j]);
  fclose (f);
  printf ("  -> %s\n", fn);
}

/* ---------------------------------------------------------------------- */
/* the report                                                             */
/* ---------------------------------------------------------------------- */

static int long_cmp (const void * a, const void * b)
{
  long x = * (const long *) a, y = * (const long *) b;
  return (x > y) - (x < y);
}

/*
.. the pairs of ranks on one core: the ranks keyed by (node, socket, core),
.. sorted, and k ranks on one key make k (k - 1) / 2 pairs. Ranks whose
.. core is not known are left out. Usually a binding mistake
*/
static long shared_cores (const MpiTopologyRank * all, int size)
{
  long * key = alloc ((size_t) size * sizeof (long)), pairs = 0;
  int n = 0;
  for (int i = 0; i < size; i++)
    if (all [i].core >= 0 && all [i].socket >= 0)
      key [n++] = ((long) all [i].node << 40) | ((long) all [i].socket << 20)
                  | all [i].core;
  qsort (key, (size_t) n, sizeof (long), long_cmp);
  for (int i = 0, run = 1; i < n; i++, run++)
    if (i + 1 == n || key [i + 1] != key [i])
    {
      pairs += (long) run * (run - 1) / 2;
      run = 0;
    }
  free (key);
  return pairs;
}

/* the ranks: a table (the first 64 of a large job) and the CSV file */
static void report_ranks (const MpiTopologyRank * all, int size,
                          const char * prefix)
{
  int show = size <= 256 ? size : 64;
  printf ("\n%5s %5s %-20s %5s %5s %6s %5s %5s %4s  %-14s %-6s %-6s %s\n",
          "rank", "node", "host", "lrank", "cpu", "socket", "numa", "core",
          "L3", "affinity", "mems", "gpus", "near");
  for (int i = 0; i < show; i++)
  {
    const MpiTopologyRank * r = &all [i];
    char gpus [MPI_TOPOLOGY_TEXT], near [MPI_TOPOLOGY_TEXT];
    gpu_list (r->gpus, gpus, sizeof gpus);
    gpu_list (r->near, near, sizeof near);
    printf ("%5d %5d %-20.20s %5d %5d %6d %5d %5d %4d  %-14s %-6s %-6s %s\n",
            r->rank, r->node, r->host, r->local_rank, r->cpu, r->socket,
            r->numa, r->core, r->l3, r->affinity, r->mems, gpus, near);
  }
  if (show < size)
    printf ("  ... %d more ranks in the CSV file\n", size - show);
  long pairs = shared_cores (all, size);
  if (pairs)
    printf ("\n  WARNING: %ld pairs of ranks run on one core now (check "
            "--cpu-bind or --bind-to)\n", pairs);

  char fn [256];
  snprintf (fn, sizeof fn, "%s_ranks.csv", prefix);
  FILE * f = fopen (fn, "w");
  if (!f)
    return;
  fprintf (f, "rank,node,host,local_rank,cpu,socket,numa,core,l3,affinity,"
           "mems,gpus,near,gpu_env\n");
  for (int i = 0; i < size; i++)
  {
    const MpiTopologyRank * r = &all [i];
    char gpus [MPI_TOPOLOGY_TEXT], near [MPI_TOPOLOGY_TEXT];
    gpu_list (r->gpus, gpus, sizeof gpus);
    gpu_list (r->near, near, sizeof near);
    fprintf (f, "%d,%d,%s,%d,%d,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\",\"%s\","
             "\"%s\"\n", r->rank, r->node, r->host, r->local_rank, r->cpu,
             r->socket, r->numa, r->core, r->l3, r->affinity, r->mems, gpus,
             near, r->gpu_env);
  }
  fclose (f);
  printf ("  -> %s\n", fn);
}

/* the allocation of every node: a table and the CSV file */
static void report_allocation (const Allocation * a, int nnodes,
                               const char * prefix)
{
  printf ("\n=== Allocation ===\n%5s %-20s %6s %5s %6s  %-16s %-8s %-6s "
          "%s\n", "node", "host", "ranks", "cpus", "cores", "cpu list",
          "numa", "mems", "gpus (of the node's)");
  for (int i = 0; i < nnodes; i++)
    printf ("%5d %-20.20s %6d %5d %6d  %-16s %-8s %-6s %s (%d of %d)\n",
            a [i].node, a [i].host, a [i].nranks, a [i].ncpus, a [i].ncores,
            a [i].cpus, a [i].numa, a [i].mems, a [i].gpus, a [i].ngpus,
            a [i].nodegpus);
  char fn [256];
  snprintf (fn, sizeof fn, "%s_alloc.csv", prefix);
  FILE * f = fopen (fn, "w");
  if (!f)
    return;
  fprintf (f, "node,host,ranks,cpus,cores,cpu_list,numa,mems,gpus,"
           "gpu_list,node_gpus\n");
  for (int i = 0; i < nnodes; i++)
    fprintf (f, "%d,%s,%d,%d,%d,\"%s\",\"%s\",\"%s\",%d,\"%s\",%d\n",
             a [i].node, a [i].host, a [i].nranks, a [i].ncpus, a [i].ncores,
             a [i].cpus, a [i].numa, a [i].mems, a [i].ngpus, a [i].gpus,
             a [i].nodegpus);
  fclose (f);
  printf ("  -> %s\n", fn);
}

/* the memory of NUMA domain d, from "Node d MemTotal: N kB"; -1 if none */
static double numa_memory_gb (int d)
{
  char path [128], line [256];
  double gb = -1;
  snprintf (path, sizeof path, SYS_NODE "/node%d/meminfo", d);
  FILE * f = fopen (path, "r");
  if (!f)
    return gb;
  while (fgets (line, sizeof line, f))
  {
    char * t = strstr (line, "MemTotal:");
    if (t)
    {
      gb = strtod (t + 9, NULL) / (1024. * 1024.);
      break;
    }
  }
  fclose (f);
  return gb;
}

/*
.. Node 0's assigned NUMA domains (those of its allocated CPUs): each
.. one's CPUs and memory, its caches in short (the instances of one level,
.. type, size and sharing counted together), and the distances between
.. the domains. A row of SYS_NODE/nodeK/distance lists every online
.. domain, in the order of SYS_NODE/online, so a domain's column is its
.. place among the online ones. Every instance in the CSV file.
*/
static void report_numa (const Mask * cpus, const Cache * c, int nc,
                         const char * prefix)
{
  Mask numa = { { 0 } }, online;
  for (int cpu = 0; cpu < NBITS; cpu++)
    if (mask_has (cpus, cpu))
      mask_add (&numa, cpu_numa (cpu));
  char line [4096], path [128];
  if (!read_line (SYS_NODE "/online", line, sizeof line)
      || !mask_parse (line, &online))
    mask_add (&online, 0);
  printf ("\n=== Node 0: the assigned NUMA domains ===\n");
  int * done = alloc ((size_t) (nc ? nc : 1) * sizeof (int));
  for (int d = 0; d < NBITS; d++)
  {
    if (!mask_has (&numa, d))
      continue;
    Mask in = { { 0 } };
    for (int cpu = 0; cpu < NBITS; cpu++)
      if (mask_has (cpus, cpu) && cpu_numa (cpu) == d)
        mask_add (&in, cpu);
    char list [MPI_TOPOLOGY_TEXT];
    mask_list (&in, list, sizeof list);
    printf ("\n  numa %d: cpus %s (%d), memory %.1f GB\n", d, list,
            mask_count (&in), numa_memory_gb (d));
    for (int i = 0; i < nc; i++)
    {
      if (done [i] || c [i].numa != d)
        continue;
      int count = 0, used = 0;
      for (int j = i; j < nc; j++)
        if (!done [j] && c [j].numa == d && c [j].level == c [i].level
            && c [j].type == c [i].type && c [j].size_kb == c [i].size_kb
            && c [j].nshared == c [i].nshared)
        {
          done [j] = 1, count++;
          used += c [j].nused;
        }
      printf ("    L%d%c %8d kB  x %-4d shared by %d cpus each, %d "
              "allocated in all\n", c [i].level,
              c [i].type == 'u' ? ' ' : c [i].type, c [i].size_kb, count,
              c [i].nshared, used);
    }
  }
  free (done);

  if (mask_count (&numa) > 1)
  {
    printf ("\n  distances (10: local):\n        ");
    for (int d = 0; d < NBITS; d++)
      if (mask_has (&numa, d))
        printf ("%5d", d);
    putchar ('\n');
    for (int d = 0; d < NBITS; d++)
    {
      if (!mask_has (&numa, d))
        continue;
      static int dist [NBITS];
      int k = 0;
      snprintf (path, sizeof path, SYS_NODE "/node%d/distance", d);
      if (read_line (path, line, sizeof line))
        for (char * s = line, * end; k < NBITS; s = end, k++)
        {
          dist [k] = (int) strtol (s, &end, 10);
          if (end == s)
            break;
        }
      printf ("  %4d  ", d);
      for (int e = 0, col = 0; e < NBITS; e++)
      {
        if (!mask_has (&online, e))
          continue;
        if (mask_has (&numa, e))
          printf ("%5d", col < k ? dist [col] : -1);
        col++;
      }
      putchar ('\n');
    }
  }

  char fn [256];
  snprintf (fn, sizeof fn, "%s_caches.csv", prefix);
  FILE * f = fopen (fn, "w");
  if (!f)
    return;
  fprintf (f, "node,numa,level,type,size_kb,shared,used\n");
  for (int i = 0; i < nc; i++)
    fprintf (f, "%d,%d,%d,%c,%d,\"%s\",\"%s\"\n", c [i].node, c [i].numa,
             c [i].level, c [i].type, c [i].size_kb, c [i].shared,
             c [i].used);
  fclose (f);
  printf ("  -> %s (node 0)\n", fn);
}

/*
.. Node 0's GPUs: each one's model, PCI address, NUMA domain, close CPUs,
.. and the ranks of node 0 that can open it and that run close to it.
.. The GPUs no rank can open are listed too, marked: they are the node's,
.. not the job's. A rank whose GPUs are all far from its CPU is warned
.. of: its data crosses the inter-socket link on every transfer.
*/
static void report_gpus (const MpiTopologyRank * all, int size,
                         const char * prefix)
{
  Gpu g [MPI_TOPOLOGY_GPUS];
  int ng = node_gpus (g, MPI_TOPOLOGY_GPUS), far = 0;
  printf ("\n=== Node 0: GPUs ===\n");
  if (!ng)
  {
    printf ("  none\n");
    return;
  }
  char fn [256];
  snprintf (fn, sizeof fn, "%s_gpus.csv", prefix);
  FILE * f = fopen (fn, "w");
  if (f)
    fprintf (f, "gpu,pci,model,uuid,file,numa,local_cpus,job,opened_by,"
             "near\n");
  printf ("  %3s %-13s %-28s %5s %-12s %-4s %-16s %s\n", "gpu", "pci",
          "model", "numa", "close cpus", "job", "opened by ranks",
          "ranks close");
  for (int i = 0; i < ng; i++)
  {
    Mask by = { { 0 } }, near = { { 0 } };
    for (int r = 0; r < size; r++)
      if (all [r].node == 0)
      {
        if ((all [r].gpus >> i) & 1)
          mask_add (&by, r);
        if ((all [r].near >> i) & 1)
          mask_add (&near, r);
      }
    char lby [MPI_TOPOLOGY_TEXT], lnear [MPI_TOPOLOGY_TEXT];
    mask_list (&by, lby, sizeof lby);
    mask_list (&near, lnear, sizeof lnear);
    int job = mask_count (&by) > 0;
    printf ("  %3d %-13s %-28.28s %5d %-12s %-4s %-16s %s\n", i, g [i].pci,
            g [i].model, g [i].numa, g [i].local, job ? "yes" : "no", lby,
            lnear);
    if (f)
      fprintf (f, "%d,%s,\"%s\",%s,%s,%d,\"%s\",%d,\"%s\",\"%s\"\n", i,
               g [i].pci, g [i].model, g [i].uuid [0] ? g [i].uuid : "-",
               g [i].file [0] ? g [i].file : "-", g [i].numa, g [i].local,
               job, lby, lnear);
  }
  for (int r = 0; r < size; r++)
    far += all [r].node == 0 && all [r].gpus && !(all [r].gpus & all [r].near);
  if (far)
    printf ("  WARNING: %d ranks of node 0 can open only GPUs far from "
            "their CPU (check --gpu-bind=closest)\n", far);
  if (f)
  {
    fclose (f);
    printf ("  -> %s\n", fn);
  }
}

/*
.. 1. the nodes: the ranks sharing memory (MPI_COMM_TYPE_SHARED), keyed by
..    rank so that local rank 0 is the node's lowest rank, its leader; the
..    leaders in a communicator of their own, keyed by rank, so that node
..    ids follow the leaders' ranks (rank 0 leads node 0);
.. 2. each rank's placement and GPUs, gathered on rank 0; the node's
..    allocation, the unions of its ranks' CPUs, memory nodes and GPUs
..    (MPI_BOR on the node), made by the leader with the caches serving
..    those CPUs; the
..    allocations and cache summaries gathered on rank 0, node 0's caches
..    kept there;
.. 3. the reports by rank 0; then the latency between every two ranks.
.. The structs travel as bytes: every rank runs the same build of the
.. same program, and MPI's heterogeneous systems are not handled.
*/
int mpi_topology (MPI_Comm comm, MpiTopologyOptions o)
{
  int rank, size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &size);
  if (!o.prefix)
    o.prefix = "topo";
  if (o.iters < 1)
    o.iters = 1;
  if (o.max_ranks < 1)
    o.max_ranks = 1024;

  /* 1. */
  MPI_Comm nodecomm, leaders;
  MPI_Comm_split_type (comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
                       &nodecomm);
  int lrank, lsize, node = 0, nnodes = 0;
  MPI_Comm_rank (nodecomm, &lrank);
  MPI_Comm_size (nodecomm, &lsize);
  MPI_Comm_split (comm, lrank == 0 ? 0 : MPI_UNDEFINED, rank, &leaders);
  if (lrank == 0)
  {
    MPI_Comm_rank (leaders, &node);
    MPI_Comm_size (leaders, &nnodes);
  }
  MPI_Bcast (&node, 1, MPI_INT, 0, nodecomm);
  MPI_Bcast (&nnodes, 1, MPI_INT, 0, comm);

  /* 2. this rank */
  void * topo = NULL;
#ifdef HAVE_HWLOC
  hwloc_topology_t t;
  hwloc_topology_init (&t);
  hwloc_topology_load (t);
  topo = &t;
#endif
  MpiTopologyRank me = { .rank = rank, .node = node, .local_rank = lrank,
                         .local_size = lsize };
  char name [MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name (name, &len);
  snprintf (me.host, sizeof me.host, "%.*s", MPI_TOPOLOGY_HOST - 1, name);
  mpi_topology_rank (&me, topo);
  MpiTopologyRank * all = NULL;
  if (rank == 0)
    all = alloc ((size_t) size * sizeof (MpiTopologyRank));
  MPI_Gather (&me, (int) sizeof me, MPI_BYTE, all, (int) sizeof me, MPI_BYTE,
              0, comm);

  /* 2. this node: the unions over its ranks */
  Mask mycpus = affinity_mask (), mymems, cpus, mems;
  mask_parse (me.mems, &mymems);
  MPI_Allreduce (mycpus.w, cpus.w, NWORDS, MPI_UNSIGNED_LONG, MPI_BOR,
                 nodecomm);
  MPI_Allreduce (mymems.w, mems.w, NWORDS, MPI_UNSIGNED_LONG, MPI_BOR,
                 nodecomm);
  unsigned long gpus;
  MPI_Allreduce (&me.gpus, &gpus, 1, MPI_UNSIGNED_LONG, MPI_BOR, nodecomm);
  Allocation * allocs = NULL;
  Cache * caches = NULL;
  int ncaches = 0;
  if (lrank == 0)
  {
    Allocation a = { .node = node, .nranks = lsize };
    snprintf (a.host, sizeof a.host, "%s", me.host);
    allocation_of (&cpus, &mems, gpus, &a);
    caches = caches_of (&cpus, node, &ncaches);
    CacheSummary s = summary_of (caches, ncaches), * sums = NULL;
    if (rank == 0)
    {
      allocs = alloc ((size_t) nnodes * sizeof (Allocation));
      sums = alloc ((size_t) nnodes * sizeof (CacheSummary));
    }
    MPI_Gather (&a, (int) sizeof a, MPI_BYTE, allocs, (int) sizeof a,
                MPI_BYTE, 0, leaders);
    MPI_Gather (&s, (int) sizeof s, MPI_BYTE, sums, (int) sizeof s, MPI_BYTE,
                0, leaders);
    if (rank == 0)
    {
      /* 3. */
      printf ("=== MPI job: %d ranks on %d nodes (%s) ===\n", size, nnodes,
#ifdef HAVE_HWLOC
              "hwloc"
#else
              "sysfs"
#endif
              );
      report_ranks (all, size, o.prefix);
      report_allocation (allocs, nnodes, o.prefix);
      report_numa (&cpus, caches, ncaches, o.prefix);
      report_gpus (all, size, o.prefix);
      int differ = 0;
      for (int i = 1; i < nnodes; i++)
        if (memcmp (&sums [i], &sums [0], sizeof (CacheSummary)))
        {
          if (!differ++)
            printf ("\n  caches unlike node 0's (instances, kB a level):");
          printf (" %d", i);
        }
      if (differ)
        printf ("\n");
      free (sums);
    }
  }
  if (rank == 0)
    fflush (stdout);                  /* before the latency's long wait */

  /* 3. the latency between ranks */
  if (o.latency && size <= o.max_ranks)
  {
    double * lat = NULL;
    if (rank == 0)
      lat = alloc ((size_t) size * size * sizeof (double));
    mpi_topology_latency (comm, o.iters, lat);
    if (rank == 0)
      report_latency (lat, all, size, o.prefix);
    free (lat);
  }
  else if (o.latency && rank == 0)
    printf ("\n  no latency: %d ranks, more than %d (--latency-ranks)\n",
            size, o.max_ranks);

  if (leaders != MPI_COMM_NULL)
    MPI_Comm_free (&leaders);
  MPI_Comm_free (&nodecomm);
#ifdef HAVE_HWLOC
  hwloc_topology_destroy (t);
#endif
  free (all), free (allocs), free (caches);
  return 0;
}

int mpi_topology_main (int argc, char ** argv)
{
  MpiTopologyOptions o = { .iters = 1000, .latency = 1, .max_ranks = 1024,
                           .prefix = "topo" };
  MPI_Init (&argc, &argv);
  for (int i = 1; i < argc; i++)
    if (!strcmp (argv [i], "--iters") && i + 1 < argc)
      o.iters = atoi (argv [++i]);
    else if (!strcmp (argv [i], "--no-latency"))
      o.latency = 0;
    else if (!strcmp (argv [i], "--latency-ranks") && i + 1 < argc)
      o.max_ranks = atoi (argv [++i]);
    else if (!strcmp (argv [i], "--out") && i + 1 < argc)
      o.prefix = argv [++i];
  int status = mpi_topology (MPI_COMM_WORLD, o);
  MPI_Finalize ();
  return status;
}

#ifdef MPI_TOPOLOGY_MAIN
int main (int argc, char ** argv)
{
  return mpi_topology_main (argc, argv);
}
#endif

#endif
