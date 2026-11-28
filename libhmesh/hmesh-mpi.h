#ifndef _HMESH_MPI_H_
#define _HMESH_MPI_H_

  /*
  .. The partitioned mesh: each process holds the simplices near its part
  .. of libb-amr's tree. The design is in doc/hmesh.tex, section
  .. "Partitioned hmesh: design"; this is its first stage, the data.
  .. Built with cwrap, -DEULERIAN_LAGRANGIAN (edges shorter than the cells:
  .. what keeps a process's needs within its leaves and those about them),
  .. on libb-amr's tree; with -D_MPI on several processes, without it on
  .. one (every simplex owned).
  ..
  .. Identity. Every point has a global id, the same on every process, for
  .. its life: g = (r0 << 44) | c, r0 the process that made it, c its
  .. counter there (up to 2^20 processes, 2^44 points each). A simplex
  .. across processes is the tuple of its points' global ids.
  ..
  .. Owners. A point's owner is the process whose leaves hold it: libb-amr's
  .. partition (tree_get_partition ()->key_first, K), by the key of the node
  .. of the deepest level at the point (a point on a face shared by two
  .. leaves goes to the upper one, as hmesh-tree.h's cache puts it). A
  .. simplex's owner is the smallest of its points' owners (the min rule).
  ..
  .. Held set. Process r holds the closure of
  ..
  ..   A_r = the simplices r owns, and those whose box, grown by eps (s) =
  ..         half the largest cell size at its points (one step's largest
  ..         move), meets one of r's leaves
  ..
  .. with every simplex its faces and points. The box test walks the Morton
  .. tree from the root, deciding a node by its key range against K: no
  .. lookup, exact.
  ..
  .. Cofaces. Each held simplex below the top keeps n_f, the number of its
  .. cofaces in the whole mesh; a half whose facet has fewer held cofaces
  .. is cut by the partition (hm_part_half_cut ()), as a boundary half is
  .. by the box.
  ..
  ..   hm_part_new (M)            the partition of M, bound to the tree
  ..                              (hm_tree_bind ()); gids from the point
  ..                              ids: a mesh the same on every process
  ..                              (made before the tree was partitioned,
  ..                              then hm_part_distribute ()). NULL if M is
  ..                              not bound or not built for the tree
  ..   hm_part_free (P)
  ..   hm_part_distribute (P)     every process keeps its held set: owners
  ..                              set, n_f counted on the whole mesh first,
  ..                              the rest removed, each level compacted
  ..                              owned first. 0, or -1 if a simplex r owns
  ..                              has a point in no leaf here (edges longer
  ..                              than the design allows)
  ..   hm_part_point_owner (P, x) the owner of a point at x
  ..   hm_part_owner (P, k, s)    a held simplex's owner, as last set
  ..   hm_part_owned (P, k)       after distribute: ids < this are owned
  ..   hm_part_gid (P, v)         a held point's global id
  ..   hm_part_local (P, g)       the local id of global id g, HM_NONE if
  ..                              not held
  ..   hm_part_cofaces (P, k, f)  n_f of a held k-simplex, k < m
  ..   hm_part_half_cut (P, k, h) 1 if the facet of the half h (of a
  ..                              k-simplex) has fewer cofaces held here
  ..                              than in the whole mesh
  ..   hm_part_box_meets (P, lo, hi)  1 if the box [lo, hi] meets one of
  ..                              this process's leaves
  ..
  .. Stage 2, a step. hm_part_advect (P, uf, c, dt, uf_prev, dt_prev):
  .. manifold_advect () (hmesh-vof.h) on the partitioned mesh, collective:
  ..   S1  each process moves the points it owns (all, or none if one
  ..       would break the CFL limit anywhere)
  ..   S2  each sends the new positions, exact, to the peers whose leaves
  ..       touch the point's old leaf
  ..   S3  owners again, by the min rule
  ..   S4  the fields of each simplex changing hands, to its new owner
  ..       (who holds it already)
  ..   S5  the refresh: each owner sends each simplex it owns to the
  ..       peers that must hold it (A_r, closed under faces), as its
  ..       points' global ids in its orientation, with its fields and
  ..       n_f; each receiver adds what it lacks, level by level (faces
  ..       first, by hm_add (), edges not checked: the move stretches
  ..       them), and drops what it neither owns nor was sent
  ..   S6  owned first again, the map; c brought up to date from the
  ..       mesh as it is then (hm_advect_end ()), in the cells of the
  ..       process's own leaves
  .. Returns manifold_advect ()'s values: -1 refused, -2 cannot run, else
  .. the owned edges longer than their cells, over every process; -3 if a
  .. message could not be applied (a bug).
  ..
  .. Stage 3, the tree changed. hm_part_tree_changed (P), collective,
  .. after adapt_wavelet (), tree_mpi_refine (), tree_mpi_coarsen () or
  .. tree_mpi_repartition (); the points stand still:
  ..   K kept (an adaptation: a leaf divided or merged keeps its key
  ..     interval's owner): owners unchanged; the leaves under the mesh
  ..     changed size, and with them the margins: one refresh (S5, S6)
  ..   K changed (a repartition):
  ..     R1  owners from the new K, the old ones from the K kept at the
  ..         last refresh
  ..     R2  each simplex with its closure, from the old owner of its
  ..         carrier (its first point the new owner owns) to its new
  ..         owner; the pairs found by MPI_Alltoall of the sizes (the
  ..         tree's own migration pairs: libb-amr does not give them)
  ..     R3  each maximal simplex with its closure, from its new owner to
  ..         the new peers its grown box reaches: owners then hold their
  ..         stars
  ..     R4  the refresh
  .. 0, or -3 if an entry could not be placed (a bug). c is not touched:
  .. after an adaptation its new leaves hold what libb-amr's prolongation
  .. or restriction gave them, until the next step.
  ..
  .. Every message carries simplices as entries: the level, the global ids
  .. in the sender's orientation, a point's position, n_f below the top,
  .. the fields; the receiver inserts them a level at a time, faces
  .. first (hm_add (), edges not checked), and finds those it holds.
  */

  #include "hmesh.h"
  #include "hmesh-tree.h"

  /* nothing without -DEULERIAN_LAGRANGIAN: the library's plain build
     compiles this file to nothing */
  #ifdef EULERIAN_LAGRANGIAN

  #define HM_GID_SHIFT 44

  typedef struct HmPart HmPart;

  HmPart * hm_part_new (Hmesh * M);
  void     hm_part_free (HmPart * P);
  int      hm_part_distribute (HmPart * P);
  int      hm_part_point_owner (const HmPart * P, const double * x);
  int      hm_part_owner (const HmPart * P, int k, uint32_t s);
  uint32_t hm_part_owned (const HmPart * P, int k);
  uint64_t hm_part_gid (const HmPart * P, uint32_t v);
  uint32_t hm_part_local (const HmPart * P, uint64_t g);
  uint32_t hm_part_cofaces (const HmPart * P, int k, uint32_t f);
  int      hm_part_half_cut (const HmPart * P, int k, uint32_t h);
  int      hm_part_box_meets (const HmPart * P, const double * lo,
                              const double * hi);

  #include "hmesh-vof.h"
  int      hm_part_advect (HmPart * P, face_vector uf, scalar c, double dt,
                           face_vector uf_prev, double dt_prev);
  int      hm_part_tree_changed (HmPart * P);

  #endif
#endif
