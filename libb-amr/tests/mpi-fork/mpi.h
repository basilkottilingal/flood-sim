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
  ..   MPI_Init, MPI_Init_thread, MPI_Initialized, MPI_Finalize,
  ..   MPI_Finalized, MPI_Comm_rank, MPI_Comm_size
  ..   MPI_Isend, MPI_Irecv, MPI_Waitall   MPI_BYTE only
  ..   MPI_Bcast, MPI_Alltoall,            any type below
  ..   MPI_Allgather, MPI_Allgatherv
  ..   MPI_Allreduce, MPI_Exscan           MPI_SUM, MPI_PROD, MPI_MIN, MPI_MAX
  ..                                       of MPI_INT, MPI_LONG, MPI_LONG_LONG,
  ..                                       MPI_UNSIGNED, MPI_UNSIGNED_LONG,
  ..                                       MPI_UNSIGNED_LONG_LONG, MPI_FLOAT,
  ..                                       MPI_DOUBLE; MPI_Allreduce also a
  ..                                       user's op (MPI_Op_create) of a
  ..                                       contiguous type
  ..   MPI_Type_contiguous, MPI_Type_commit, MPI_Type_free
  ..                                       count items of a type above, as one
  ..   MPI_Op_create, MPI_Op_free          a user's function, commutative
  ..   MPI_File_open, MPI_File_set_size,   on the file system the forked
  ..   MPI_File_get_size,                  processes share: open, ftruncate,
  ..   MPI_File_write_at,                  fstat, pwrite, pread, close;
  ..   MPI_File_read_at, MPI_File_close    MPI_BYTE only
  */

  typedef int MPI_Comm;
  typedef int MPI_Datatype;
  typedef int MPI_Op;
  typedef struct { int dummy; } MPI_Status;
  typedef int MPI_Request;
  typedef int MPI_File;
  typedef int MPI_Info;
  typedef long long MPI_Offset;
  #define MPI_SUCCESS         0
  #define MPI_INFO_NULL       0
  #define MPI_MODE_CREATE     1
  #define MPI_MODE_RDONLY     2
  #define MPI_MODE_WRONLY     4

  #define MPI_COMM_WORLD      0
  #define MPI_BYTE            1
  #define MPI_INT             2
  #define MPI_LONG            3
  #define MPI_DOUBLE          4
  #define MPI_LONG_LONG       5
  #define MPI_UNSIGNED        6
  #define MPI_UNSIGNED_LONG   7
  #define MPI_UNSIGNED_LONG_LONG 8
  #define MPI_FLOAT           9
  #define MPI_SUM             1
  #define MPI_PROD            2
  #define MPI_MIN             3
  #define MPI_MAX             4
  #define MPI_STATUSES_IGNORE ((MPI_Status *) 0)
  #define MPI_THREAD_SINGLE   0
  #define MPI_THREAD_FUNNELED 1

  int MPI_Init (int * argc, char *** argv);
  int MPI_Init_thread (int * argc, char *** argv, int required,
                       int * provided);
  int MPI_Initialized (int * flag);
  int MPI_Finalize (void);
  int MPI_Finalized (int * flag);
  int MPI_Comm_rank (MPI_Comm comm, int * rank);
  int MPI_Comm_size (MPI_Comm comm, int * size);
  int MPI_Isend (const void * buf, int count, MPI_Datatype type, int dest,
                 int tag, MPI_Comm comm, MPI_Request * req);
  int MPI_Irecv (void * buf, int count, MPI_Datatype type, int source,
                 int tag, MPI_Comm comm, MPI_Request * req);
  int MPI_Waitall (int n, MPI_Request * req, MPI_Status * status);
  int MPI_Bcast (void * buf, int count, MPI_Datatype type, int root,
                 MPI_Comm comm);
  int MPI_Alltoall (const void * sbuf, int scount, MPI_Datatype stype,
                    void * rbuf, int rcount, MPI_Datatype rtype,
                    MPI_Comm comm);
  int MPI_Exscan (const void * in, void * out, int count, MPI_Datatype type,
                  MPI_Op op, MPI_Comm comm);
  int MPI_Allreduce (const void * in, void * out, int count,
                     MPI_Datatype type, MPI_Op op, MPI_Comm comm);
  int MPI_Allgather (const void * sbuf, int scount, MPI_Datatype stype,
                     void * rbuf, int rcount, MPI_Datatype rtype,
                     MPI_Comm comm);
  int MPI_Allgatherv (const void * sbuf, int scount, MPI_Datatype stype,
                      void * rbuf, const int * rcounts, const int * displs,
                      MPI_Datatype rtype, MPI_Comm comm);

  typedef void MPI_User_function (void * in, void * inout, int * len,
                                  MPI_Datatype * type);
  int MPI_Type_contiguous (int count, MPI_Datatype old, MPI_Datatype * type);
  int MPI_Type_commit (MPI_Datatype * type);
  int MPI_Type_free (MPI_Datatype * type);
  int MPI_Op_create (MPI_User_function * f, int commute, MPI_Op * op);
  int MPI_Op_free (MPI_Op * op);
  int MPI_File_open (MPI_Comm comm, const char * name, int amode,
                     MPI_Info info, MPI_File * fh);
  int MPI_File_set_size (MPI_File fh, MPI_Offset size);
  int MPI_File_write_at (MPI_File fh, MPI_Offset offset, const void * buf,
                         int count, MPI_Datatype type, MPI_Status * status);
  int MPI_File_get_size (MPI_File fh, MPI_Offset * size);
  int MPI_File_read_at (MPI_File fh, MPI_Offset offset, void * buf,
                        int count, MPI_Datatype type, MPI_Status * status);
  int MPI_File_close (MPI_File * fh);

#endif
