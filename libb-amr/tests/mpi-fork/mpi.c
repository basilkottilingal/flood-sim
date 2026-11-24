/* The stand-in for MPI of the tests (mpi.h): forked ranks, a socket pair
 * between each two, messages framed with their tag and length. */
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <unistd.h>

#include "mpi.h"

enum { NP_MAX = 16, REQ_MAX = 1024 };

static int me = 0, np = 1, fd [NP_MAX], initialized = 0, finalized = 0;
static pid_t child [NP_MAX];

typedef struct
{
  int used, send, peer, tag;
  char * buf;
  size_t len, done;             /* done: of the header and the payload */
  int head [2];                 /* tag, length */
} Req;

static Req reqs [REQ_MAX];

static void fail (const char * what)
{
  fprintf (stderr, "mpi-fork, rank %d: %s\n", me, what);
  abort ();
}

int MPI_Init (int * argc, char *** argv)
{
  (void) argc, (void) argv;
  const char * e = getenv ("MPI_FORK_NP");
  np = e ? atoi (e) : 1;
  if (np < 1 || np > NP_MAX)
    fail ("MPI_FORK_NP out of range");
  int sp [NP_MAX][NP_MAX][2];
  for (int i = 0; i < np; i++)
    for (int j = i + 1; j < np; j++)
      if (socketpair (AF_UNIX, SOCK_STREAM, 0, sp [i][j]))
        fail ("socketpair");
  fflush (stdout), fflush (stderr);
  for (int r = 1; r < np; r++)
  {
    pid_t pid = fork ();
    if (pid < 0)
      fail ("fork");
    if (pid == 0)
    {
      me = r;
      break;
    }
    child [r] = pid;
  }
  for (int i = 0; i < np; i++)
    for (int j = i + 1; j < np; j++)
    {
      int mine = i == me ? 0 : j == me ? 1 : -1;
      for (int s = 0; s < 2; s++)
        if (s != mine)
          close (sp [i][j][s]);
      if (mine >= 0)
      {
        int peer = i == me ? j : i, f = sp [i][j][mine];
        fcntl (f, F_SETFL, fcntl (f, F_GETFL) | O_NONBLOCK);
        fd [peer] = f;
      }
    }
  fd [me] = -1;
  initialized = 1;
  return 0;
}

int MPI_Initialized (int * flag)
{
  *flag = initialized;
  return 0;
}

/* forked processes: any threading the caller does is its own */
int MPI_Init_thread (int * argc, char *** argv, int required, int * provided)
{
  *provided = required;
  return MPI_Init (argc, argv);
}

int MPI_Finalized (int * flag)
{
  *flag = finalized;
  return 0;
}

int MPI_Finalize (void)
{
  finalized = 1;
  for (int r = 0; r < np; r++)
    if (r != me)
      close (fd [r]);
  /* rank 0 waits for the others, and fails if one did: also when called
   * at exit (atexit ()), where only _exit () may end the run */
  int failed = 0;
  if (me == 0)
    for (int r = 1; r < np; r++)
    {
      int status = 0;
      waitpid (child [r], &status, 0);
      if (!WIFEXITED (status) || WEXITSTATUS (status))
      {
        fprintf (stderr, "mpi-fork: rank %d failed\n", r);
        failed = 1;
      }
    }
  if (failed)
  {
    fflush (stdout), fflush (stderr);
    _exit (1);
  }
  return 0;
}

int MPI_Comm_rank (MPI_Comm comm, int * rank)
{
  (void) comm;
  *rank = me;
  return 0;
}

int MPI_Comm_size (MPI_Comm comm, int * size)
{
  (void) comm;
  *size = np;
  return 0;
}

static int post (int send, void * buf, int count, MPI_Datatype type, int peer,
                 int tag, MPI_Request * req)
{
  if (type != MPI_BYTE || count < 0 || peer < 0 || peer >= np || peer == me)
    fail ("a message other than bytes to another rank");
  for (int q = 0; q < REQ_MAX; q++)
    if (!reqs [q].used)
    {
      reqs [q] = (Req) { 1, send, peer, tag, buf, (size_t) count, 0,
                         { tag, count } };
      *req = q;
      return 0;
    }
  fail ("too many requests");
  return 1;
}

int MPI_Isend (const void * buf, int count, MPI_Datatype type, int dest,
               int tag, MPI_Comm comm, MPI_Request * req)
{
  (void) comm;
  return post (1, (void *) buf, count, type, dest, tag, req);
}

int MPI_Irecv (void * buf, int count, MPI_Datatype type, int source,
               int tag, MPI_Comm comm, MPI_Request * req)
{
  (void) comm;
  return post (0, buf, count, type, source, tag, req);
}

/* moves r on as far as its socket lets it: 1 once complete */
static int progress (Req * r)
{
  size_t hl = sizeof (r->head), total = hl + r->len;
  while (r->done < total)
  {
    char * at = r->done < hl ? (char *) r->head + r->done
                             : r->buf + (r->done - hl);
    size_t want = r->done < hl ? hl - r->done : total - r->done;
    ssize_t n = r->send ? write (fd [r->peer], at, want)
                        : read (fd [r->peer], at, want);
    if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR))
      return 0;
    if (n <= 0)
      fail (r->send ? "a send failed" : "a peer closed before sending");
    r->done += (size_t) n;
    if (!r->send && r->done == hl
        && (r->head [0] != r->tag || r->head [1] != (int) r->len))
    {
      char m [160];
      snprintf (m, sizeof m, "from rank %d: tag %d, %d bytes; expected tag "
                "%d, %d bytes", r->peer, r->head [0], r->head [1], r->tag,
                (int) r->len);
      fail (m);
    }
  }
  return 1;
}

int MPI_Waitall (int n, MPI_Request * req, MPI_Status * status)
{
  (void) status;
  for (;;)
  {
    /* the first unfinished request of each peer and way: streams in order */
    struct pollfd pf [2 * NP_MAX];
    int npf = 0, left = 0;
    int taken [NP_MAX][2] = { { 0 } };
    for (int q = 0; q < n; q++)
    {
      Req * r = &reqs [req [q]];
      if (!r->used)
        continue;
      if (progress (r))
      {
        r->used = 0;
        continue;
      }
      left++;
      if (taken [r->peer][r->send])
        continue;
      taken [r->peer][r->send] = 1;
      pf [npf] = (struct pollfd) { fd [r->peer],
                                   r->send ? POLLOUT : POLLIN, 0 };
      npf++;
    }
    if (!left)
      return 0;
    if (poll (pf, (nfds_t) npf, 60000) <= 0)
      fail ("no progress in 60 s");
  }
}

/* the types made by MPI_Type_contiguous (): ids from TYPE_FIRST on, their
 * sizes; the ops by MPI_Op_create (): ids from OP_FIRST on, the functions */
#define TYPE_FIRST 100
#define TYPE_MAX   32
#define OP_FIRST   100
#define OP_MAX     32
static size_t types [TYPE_MAX];
static MPI_User_function * ops [OP_MAX];

/* the size of one item of type, 0 if not one of the stand-in's */
static size_t type_size (MPI_Datatype type)
{
  if (type >= TYPE_FIRST && type < TYPE_FIRST + TYPE_MAX)
    return types [type - TYPE_FIRST];
  switch (type)
  {
  case MPI_BYTE:               return 1;
  case MPI_INT:                return sizeof (int);
  case MPI_LONG:               return sizeof (long);
  case MPI_DOUBLE:             return sizeof (double);
  case MPI_LONG_LONG:          return sizeof (long long);
  case MPI_UNSIGNED:           return sizeof (unsigned);
  case MPI_UNSIGNED_LONG:      return sizeof (unsigned long);
  case MPI_UNSIGNED_LONG_LONG: return sizeof (unsigned long long);
  case MPI_FLOAT:              return sizeof (float);
  default:                     return 0;
  }
}

/* root sends to every other rank in turn */
int MPI_Bcast (void * buf, int count, MPI_Datatype type, int root,
               MPI_Comm comm)
{
  size_t size = type_size (type);
  if (!size || count < 0 || root < 0 || root >= np)
    fail ("Bcast: a type of the stand-in's, from a rank");
  int bytes = (int) ((size_t) count * size), tag = -2;
  MPI_Request q;
  if (me == root)
    for (int r = 0; r < np; r++)
    {
      if (r == root)
        continue;
      MPI_Isend (buf, bytes, MPI_BYTE, r, tag, comm, &q);
      MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
    }
  else
  {
    MPI_Irecv (buf, bytes, MPI_BYTE, root, tag, comm, &q);
    MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
  }
  return 0;
}

/* block r of sbuf to rank r, block r of rbuf from it; one's own copied */
int MPI_Alltoall (const void * sbuf, int scount, MPI_Datatype stype,
                  void * rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm)
{
  size_t size = type_size (stype);
  if (!size || stype != rtype || scount != rcount || scount < 0)
    fail ("Alltoall: the same type and count both ways");
  int bytes = (int) ((size_t) scount * size), tag = -3, nr = 0;
  MPI_Request q [2 * NP_MAX];
  for (int r = 0; r < np; r++)
    if (r == me)
      memcpy ((char *) rbuf + (size_t) r * bytes,
              (const char *) sbuf + (size_t) r * bytes, (size_t) bytes);
    else
    {
      MPI_Irecv ((char *) rbuf + (size_t) r * bytes, bytes, MPI_BYTE, r, tag,
                 comm, &q [nr++]);
      MPI_Isend ((const char *) sbuf + (size_t) r * bytes, bytes, MPI_BYTE,
                 r, tag, comm, &q [nr++]);
    }
  MPI_Waitall (nr, q, MPI_STATUSES_IGNORE);
  return 0;
}

/* out [k] = out [k] op in [k], k < count, for items of type T */
#define COMBINE(T)                                                      \
  for (int k = 0; k < count; k++)                                       \
  {                                                                     \
    T a = ((T *) out) [k], b = ((const T *) in) [k];                    \
    ((T *) out) [k] = op == MPI_SUM ? a + b : op == MPI_PROD ? a * b    \
                    : op == MPI_MIN ? (b < a ? b : a) : (b > a ? b : a);\
  }

static void combine (void * out, const void * in, int count,
                     MPI_Datatype type, MPI_Op op)
{
  switch (type)
  {
  case MPI_INT:                COMBINE (int); break;
  case MPI_LONG:               COMBINE (long); break;
  case MPI_DOUBLE:             COMBINE (double); break;
  case MPI_LONG_LONG:          COMBINE (long long); break;
  case MPI_UNSIGNED:           COMBINE (unsigned); break;
  case MPI_UNSIGNED_LONG:      COMBINE (unsigned long); break;
  case MPI_UNSIGNED_LONG_LONG: COMBINE (unsigned long long); break;
  case MPI_FLOAT:              COMBINE (float); break;
  default:                     fail ("Allreduce: a type of the stand-in's");
  }
}

/* rank r gets the values of ranks 0 .. r - 1 combined, in rank order;
 * rank 0's out is left as it was (undefined, as MPI's) */
int MPI_Exscan (const void * in, void * out, int count, MPI_Datatype type,
                MPI_Op op, MPI_Comm comm)
{
  size_t size = type_size (type);
  if (op < MPI_SUM || op > MPI_MAX || !size || type == MPI_BYTE || count < 0)
    fail ("Exscan: MPI_SUM, MPI_PROD, MPI_MIN or MPI_MAX of a type of the "
          "stand-in's");
  int bytes = (int) ((size_t) count * size), tag = -4;
  char * buf = malloc ((size_t) bytes + 1);
  MPI_Request q;
  for (int r = me + 1; r < np; r++)
  {
    MPI_Isend (in, bytes, MPI_BYTE, r, tag, comm, &q);
    MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
  }
  for (int r = 0; r < me; r++)
  {
    MPI_Irecv (buf, bytes, MPI_BYTE, r, tag, comm, &q);
    MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
    if (r == 0)
      memcpy (out, buf, (size_t) bytes);
    else
      combine (out, buf, count, type, op);
  }
  free (buf);
  return 0;
}

int MPI_Type_contiguous (int count, MPI_Datatype old, MPI_Datatype * type)
{
  size_t size = type_size (old);
  int t = 0;
  while (t < TYPE_MAX && types [t])
    t++;
  if (!size || count < 1 || t == TYPE_MAX)
    fail ("Type_contiguous: of a type of the stand-in's, and room");
  types [t] = (size_t) count * size;
  *type = TYPE_FIRST + t;
  return 0;
}

int MPI_Type_commit (MPI_Datatype * type)
{
  if (!type_size (*type))
    fail ("Type_commit: a type made");
  return 0;
}

int MPI_Type_free (MPI_Datatype * type)
{
  if (*type < TYPE_FIRST || *type >= TYPE_FIRST + TYPE_MAX
      || !types [*type - TYPE_FIRST])
    fail ("Type_free: a type made");
  types [*type - TYPE_FIRST] = 0;
  *type = 0;
  return 0;
}

int MPI_Op_create (MPI_User_function * f, int commute, MPI_Op * op)
{
  int o = 0;
  while (o < OP_MAX && ops [o])
    o++;
  if (!f || !commute || o == OP_MAX)
    fail ("Op_create: a commutative function, and room");
  ops [o] = f;
  *op = OP_FIRST + o;
  return 0;
}

int MPI_Op_free (MPI_Op * op)
{
  if (*op < OP_FIRST || *op >= OP_FIRST + OP_MAX || !ops [*op - OP_FIRST])
    fail ("Op_free: an op made");
  ops [*op - OP_FIRST] = NULL;
  *op = 0;
  return 0;
}

/* out = in op out for a user's op, MPI's order of the arguments */
static void combine_user (void * out, void * in, int count,
                          MPI_Datatype type, MPI_Op op)
{
  ops [op - OP_FIRST] (in, out, &count, &type);
}

/* every rank's to rank 0, combined in rank order, sent back */
int MPI_Allreduce (const void * in, void * out, int count, MPI_Datatype type,
                   MPI_Op op, MPI_Comm comm)
{
  size_t size = type_size (type);
  int user = op >= OP_FIRST && op < OP_FIRST + OP_MAX && ops [op - OP_FIRST];
  if ((!user && (op < MPI_SUM || op > MPI_MAX || type >= TYPE_FIRST))
      || !size || type == MPI_BYTE || count < 0)
    fail ("Allreduce: MPI_SUM, MPI_PROD, MPI_MIN or MPI_MAX of a type of "
          "the stand-in's, or a user's op");
  int bytes = (int) ((size_t) count * size), tag = -1;
  memcpy (out, in, (size_t) bytes);
  char * buf = malloc ((size_t) bytes + 1);
  MPI_Request q;
  if (me == 0)
  {
    for (int r = 1; r < np; r++)
    {
      MPI_Irecv (buf, bytes, MPI_BYTE, r, tag, comm, &q);
      MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
      if (user)
        combine_user (out, buf, count, type, op);
      else
        combine (out, buf, count, type, op);
    }
    for (int r = 1; r < np; r++)
    {
      MPI_Isend (out, bytes, MPI_BYTE, r, tag, comm, &q);
      MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
    }
  }
  else
  {
    MPI_Isend (in, bytes, MPI_BYTE, 0, tag, comm, &q);
    MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
    MPI_Irecv (out, bytes, MPI_BYTE, 0, tag, comm, &q);
    MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
  }
  free (buf);
  return 0;
}

/* block r of rbuf (rcounts [r] items at displs [r]) from rank r, one's
 * own copied: every rank's to every other, in turn */
int MPI_Allgatherv (const void * sbuf, int scount, MPI_Datatype stype,
                    void * rbuf, const int * rcounts, const int * displs,
                    MPI_Datatype rtype, MPI_Comm comm)
{
  size_t size = type_size (stype);
  if (!size || stype != rtype || scount < 0 || scount != rcounts [me])
    fail ("Allgatherv: the same type both ways, this rank's count");
  int tag = -5, nr = 0;
  MPI_Request q [2 * NP_MAX];
  memcpy ((char *) rbuf + (size_t) displs [me] * size, sbuf,
          (size_t) scount * size);
  for (int r = 0; r < np; r++)
    if (r != me)
    {
      MPI_Irecv ((char *) rbuf + (size_t) displs [r] * size,
                 (int) ((size_t) rcounts [r] * size), MPI_BYTE, r, tag, comm,
                 &q [nr++]);
      MPI_Isend (sbuf, (int) ((size_t) scount * size), MPI_BYTE, r, tag,
                 comm, &q [nr++]);
    }
  MPI_Waitall (nr, q, MPI_STATUSES_IGNORE);
  return 0;
}

/* MPI_Allgatherv () with count items from each rank, in rank order */
int MPI_Allgather (const void * sbuf, int scount, MPI_Datatype stype,
                   void * rbuf, int rcount, MPI_Datatype rtype, MPI_Comm comm)
{
  int counts [NP_MAX], displs [NP_MAX];
  if (scount != rcount)
    fail ("Allgather: the same count both ways");
  for (int r = 0; r < np; r++)
    counts [r] = rcount, displs [r] = r * rcount;
  return MPI_Allgatherv (sbuf, scount, stype, rbuf, counts, displs, rtype,
                         comm);
}
