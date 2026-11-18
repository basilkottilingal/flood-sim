/*
 * mpi-topology.c — discover the topology of an MPI job on a cluster/supercomputer.
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
 *   mpicc -O2 -o mpi-topology mpi-topology.c    # sysfs only (Linux)
 *   mpicc -O2 -DHAVE_HWLOC -o mpi-topology mpi-topology.c $ (pkg-config --cflags --libs hwloc)  # with hwloc >= 2.0
 *
 * Run (examples)
 *   srun -N 8 --ntasks-per-node=4 --cpu-bind=cores ./mpi-topology
 *   mpirun -np 64 --map-by ppr:8:node --bind-to core ./mpi-topology --iters 2000
 *
 * Options
 *   --iters N    ping-pong iterations per node pair (default 1000)
 *   --no-latency   skip the network latency measurement
 *   --out PREFIX   CSV file prefix (default "topo") -> PREFIX_ranks.csv,
 *          PREFIX_nodes.csv, PREFIX_latency.csv
 */

#ifdef _MPI
#include <mpi.h>
typedef struct
{
  int  rank, node_id, local_rank, local_size;
  int  cpu, socket, numa, core, l3;
  char host[HOSTLEN];
  char affinity[STRLEN];
} rank_info_t;

typedef struct
{
  int  node_id, nranks;
  int  sockets, numa, cores, pus, l3s, gpus, nnics;
  double mem_gb;
  char   host[HOSTLEN];
  char   nics[STRLEN];
} node_info_t;

/* api  */
int  mpi_topology     (int argc, char **argv);

/* lower level api(s) */
void cpuset_to_string (const cpu_set_t *set, char *buf, size_t len);
int  sysfs_cpu_numa   (int cpu);
int  sysfs_count      (const char *pattern);
int  sysfs_count_gpus (void);
int  sysfs_list_nics  (char *buf, size_t len);
void fill_rank_info   (rank_info_t *ri, void *topo_ptr);
void fill_node_info   (node_info_t *ni, void *topo_ptr);
void print_numa_distances (void);
void measure_latency  (MPI_Comm leaders, int iters, double *mat);
void report_latency   (const double *mat, const node_info_t *nodes,
        int n, const char *prefix);
#endif
