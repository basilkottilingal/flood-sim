#ifndef _MPI_TOPOLOGY_H_
#define _MPI_TOPOLOGY_H_

  /*
  .. The part of a machine an MPI job was given: where each rank runs, the
  .. CPUs, NUMA domains and GPUs its ranks may use, the caches serving
  .. those CPUs, and the latency between every two ranks. Experimental: a tool
  .. to look at an allocation before choosing a partition or a binding, not
  .. part of libb-amr's tree (tree.h, tree-mpi.h). Linux only (sysfs,
  .. /proc); hwloc 2 for the placement if built with -DHAVE_HWLOC.
  ..
  .. Only what the job was given. On a cluster, sysfs describes the whole
  .. node, while the batch system confines a job to a cgroup: some CPUs,
  .. some memory nodes. So the allocation is taken from the ranks
  .. themselves: a node's CPUs are the union of its ranks' affinities
  .. (sched_getaffinity ()), its memory nodes the union of their
  .. Mems_allowed_list (/proc/self/status). A CPU no rank may run on is
  .. not reported, nor a cache serving none of the allocated CPUs.
  ..
  .. GPUs. A GPU node has a few (4 to 8 devices), each attached to one
  .. socket or NUMA domain: close to some cores (sysfs: numa_node,
  .. local_cpulist of its PCI device). A job is given some of them by the
  .. cgroup's device rules: it can open those device files only. The
  .. *_VISIBLE_DEVICES variables a batch system sets number devices in the
  .. vendor library's order, which is not always the PCI order, so they
  .. are recorded as they are, not mapped. A node's GPUs are its PCI
  .. display controllers and processing accelerators (classes 0x03, 0x12)
  .. of NVIDIA, AMD or Intel (an integrated Intel GPU counts), numbered 0,
  .. 1, ... in PCI order; a rank's GPUs are those whose device file it can
  .. open (/dev/nvidiaN, the minor from /proc/driver/nvidia; else the
  .. render node /dev/dri/renderDN; for AMD /dev/kfd too, ROCm's compute
  .. driver), the node's allocation their union. A GPU without its driver
  .. has no device file: listed, never the job's. Opening a device file
  .. starts the driver if no one holds it (some 100 ms); nothing else is
  .. done with it.
  ..
  .. What mpi_topology () reports, rank 0 printing and writing CSV files:
  ..   ranks       host, node, local rank, the CPU it runs on now, its
  ..               socket, NUMA domain, core and L3, its CPU affinity and
  ..               memory nodes, the GPUs it can open and those close to
  ..               its CPU, its *_VISIBLE_DEVICES: PREFIX_ranks.csv
  ..   allocation  for each node: its ranks, CPUs, cores, the NUMA domains
  ..               of those CPUs, the memory nodes, the GPUs:
  ..               PREFIX_alloc.csv
  ..   GPUs        node 0's GPUs: model, PCI address, NUMA domain, close
  ..               CPUs, the ranks that can open each and those close to
  ..               it; a warning for a rank whose GPUs are all far from
  ..               its CPU: PREFIX_gpus.csv
  ..   NUMA        node 0's assigned NUMA domains: their CPUs, memory, the
  ..               distances between them (the kernel's, from the ACPI SLIT
  ..               table; 10 is local), and their caches: for each level,
  ..               the size, the CPUs sharing an instance and how many of
  ..               them are allocated. Every instance: PREFIX_caches.csv;
  ..               the nodes whose caches differ from node 0's named
  ..   latency     between every two ranks, by ping-pong of 8 bytes, one
  ..               way, by class: the same core (SMT threads), the same L3,
  ..               the same NUMA domain, the same socket, the same node,
  ..               another node: PREFIX_latency.csv, a pair a line
  ..
  .. A node is the set of ranks sharing memory (MPI_COMM_TYPE_SHARED); its
  .. leader is its lowest rank, and nodes are numbered in the order of
  .. their leaders. A rank not bound to one CPU may move: its cpu, core and
  .. L3, and so the class of its pairs, are where it ran when asked.
  ..
  .. Build (the whole file is #ifdef _MPI, so the library builds without
  .. MPI):
  ..   mpicc -O2 -D_MPI -DMPI_TOPOLOGY_MAIN mpi-topology.c -o mpi-topology
  ..   mpicc -O2 -D_MPI -DMPI_TOPOLOGY_MAIN -DHAVE_HWLOC mpi-topology.c \
  ..     -o mpi-topology $(pkg-config --cflags --libs hwloc)
  .. Run:
  ..   srun -N 8 --ntasks-per-node=4 --cpu-bind=cores ./mpi-topology
  ..   mpirun -np 64 --map-by ppr:8:node --bind-to core ./mpi-topology
  .. Options: --iters N (round trips a pair, 1000), --no-latency,
  .. --latency-ranks N (no latency above N ranks, 1024), --out PREFIX (the
  .. CSV files' prefix, "topo").
  ..
  .. API
  ..   mpi_topology (comm, options)      collective on comm, MPI started by
  ..                                     the caller: everything above, rank
  ..                                     0 of comm reporting. 0
  ..   mpi_topology_main (argc, argv)    the program: MPI_Init, options
  ..                                     from argv, mpi_topology () on
  ..                                     MPI_COMM_WORLD, MPI_Finalize
  .. The parts, for a caller that wants the data rather than the report:
  ..   mpi_topology_rank (r, topo)       this rank's placement and GPUs
  ..                                     into r
  ..   mpi_topology_latency (comm, iters, mat)
  ..                                     collective on comm: the one-way
  ..                                     latency between every two ranks,
  ..                                     n x n (us), on rank 0 (mat unread
  ..                                     elsewhere); O(n) memory a rank
  ..   mpi_topology_cpuset (set, buf, len)
  ..                                     a CPU set (a cpu_set_t *: void *
  ..                                     here, so that the header needs no
  ..                                     _GNU_SOURCE) as "0-3,8,10-11"
  ..   mpi_topology_node (n, topo)       the whole node's resources into n
  ..                                     (sockets, cores, memory, GPUs,
  ..                                     NICs), allocated or not; not used
  ..                                     by mpi_topology ()
  .. topo is a hwloc_topology_t * with -DHAVE_HWLOC, else NULL. Fields a
  .. source cannot tell are -1.
  ..
  .. Not yet: the mapping of *_VISIBLE_DEVICES to devices (the vendor's
  .. library), other vendors' accelerators (Habana Gaudi's /dev/accel,
  .. Huawei Ascend's /dev/davinci), the cgroup's
  .. memory limit, a latency matrix beyond a few thousand ranks (sampled
  .. pairs), bandwidth.
  */

  #ifdef _MPI

  #include <stddef.h>
  #include <mpi.h>

  enum
  {
    MPI_TOPOLOGY_HOST = 64,     /* a host name, its end included       */
    MPI_TOPOLOGY_TEXT = 160,    /* a list of CPUs, nodes or NICs       */
    MPI_TOPOLOGY_GPUS = 64,     /* GPUs a node at most (a bit each)    */
  };

  /*
  .. Where a rank runs. The core is hwloc's logical index (unique on the
  .. node), or sysfs's core_id (unique within a socket only: compare the
  .. pair socket, core). The L3 is hwloc's logical index, or sysfs's cache
  .. id.
  */
  typedef struct
  {
    int rank, node, local_rank, local_size;
    int cpu, socket, numa, core, l3;
    char host [MPI_TOPOLOGY_HOST];
    char affinity [MPI_TOPOLOGY_TEXT];   /* the CPUs it may run on     */
    char mems [MPI_TOPOLOGY_TEXT];       /* the memory nodes it may use */
    unsigned long gpus;                  /* bit g: it can open GPU g   */
    unsigned long near;                  /* bit g: GPU g is close to
                                            its CPU (local_cpulist)    */
    char gpu_env [MPI_TOPOLOGY_TEXT];    /* NAME=value;... of the
                                            *_VISIBLE_DEVICES set      */
  } MpiTopologyRank;

  /* what a whole node has, allocated or not (mpi_topology_node ()) */
  typedef struct
  {
    int node, nranks;
    int sockets, numa, cores, pus, l3s, gpus, nnics;
    double mem_gb;
    char host [MPI_TOPOLOGY_HOST];
    char nics [MPI_TOPOLOGY_TEXT];       /* names, separated by spaces */
  } MpiTopologyNode;

  typedef struct
  {
    int iters;                  /* round trips a pair of ranks         */
    int latency;                /* 0: no latency measured              */
    int max_ranks;              /* no latency above (0: 1024)          */
    const char * prefix;        /* of the CSV files (NULL: "topo")     */
  } MpiTopologyOptions;

  int  mpi_topology         (MPI_Comm comm, MpiTopologyOptions options);
  int  mpi_topology_main    (int argc, char ** argv);
  void mpi_topology_rank    (MpiTopologyRank * r, void * topo);
  void mpi_topology_latency (MPI_Comm comm, int iters, double * mat);
  void mpi_topology_cpuset  (const void * set, char * buf, size_t len);
  void mpi_topology_node    (MpiTopologyNode * n, void * topo);

  #endif

#endif
