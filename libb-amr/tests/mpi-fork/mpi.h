#ifndef _MPI_FORK_H_
#define _MPI_FORK_H_

  /*
  .. A stand-in for MPI, for the tests only: the calls tree-mpi.c and
  .. tests/test_mpi.c make, on processes forked by MPI_Init () and joined by
  .. a socket pair each, where no MPI launcher works (an MPICH built for
  .. PMIx with no PMIx launcher, as Ubuntu's: every rank starts alone).
  ..
  .. MPI_Init () forks MPI_FORK_NP - 1 processes (the environment, 1 if
  .. unset); MPI_Finalize () in rank 0 waits for the others and fails if one
  .. did. Messages between two ranks keep their order (a stream each way),
  .. as MPI's between a pair with one tag do; each carries its tag and
  .. length, and a receive that does not match them stops the run, so a
  .. plan whose two sides disagree is caught. Isend and Irecv progress in
  .. MPI_Waitall () (poll), so no send blocks another's receive.
  ..
  ..   MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size
  ..   MPI_Isend, MPI_Irecv, MPI_Waitall   MPI_BYTE only
  ..   MPI_Allreduce                       MPI_SUM of MPI_INT, MPI_LONG
  */

  typedef int MPI_Comm;
  typedef int MPI_Datatype;
  typedef int MPI_Op;
  typedef struct { int dummy; } MPI_Status;
  typedef int MPI_Request;

  #define MPI_COMM_WORLD      0
  #define MPI_BYTE            1
  #define MPI_INT             2
  #define MPI_LONG            3
  #define MPI_SUM             1
  #define MPI_STATUSES_IGNORE ((MPI_Status *) 0)

  int MPI_Init (int * argc, char *** argv);
  int MPI_Finalize (void);
  int MPI_Comm_rank (MPI_Comm comm, int * rank);
  int MPI_Comm_size (MPI_Comm comm, int * size);
  int MPI_Isend (const void * buf, int count, MPI_Datatype type, int dest,
                 int tag, MPI_Comm comm, MPI_Request * req);
  int MPI_Irecv (void * buf, int count, MPI_Datatype type, int source,
                 int tag, MPI_Comm comm, MPI_Request * req);
  int MPI_Waitall (int n, MPI_Request * req, MPI_Status * status);
  int MPI_Allreduce (const void * in, void * out, int count,
                     MPI_Datatype type, MPI_Op op, MPI_Comm comm);

#endif
