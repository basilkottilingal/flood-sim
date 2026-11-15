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

static int me = 0, np = 1, fd [NP_MAX], initialized = 0;
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

int MPI_Finalize (void)
{
  for (int r = 0; r < np; r++)
    if (r != me)
      close (fd [r]);
  if (me == 0)
    for (int r = 1; r < np; r++)
    {
      int status = 0;
      waitpid (child [r], &status, 0);
      if (!WIFEXITED (status) || WEXITSTATUS (status))
        fprintf (stderr, "mpi-fork: rank %d failed\n", r);
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

/* root sends to every other rank in turn */
int MPI_Bcast (void * buf, int count, MPI_Datatype type, int root,
               MPI_Comm comm)
{
  size_t size = type == MPI_BYTE ? 1 : type == MPI_INT ? sizeof (int)
              : type == MPI_LONG ? sizeof (long)
              : type == MPI_DOUBLE ? sizeof (double) : 0;
  if (!size || count < 0 || root < 0 || root >= np)
    fail ("Bcast: MPI_BYTE, MPI_INT, MPI_LONG or MPI_DOUBLE from a rank");
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

int MPI_Allreduce (const void * in, void * out, int count, MPI_Datatype type,
                   MPI_Op op, MPI_Comm comm)
{
  (void) comm;
  if (op != MPI_SUM || (type != MPI_INT && type != MPI_LONG))
    fail ("Allreduce: MPI_SUM of MPI_INT or MPI_LONG only");
  size_t size = type == MPI_INT ? sizeof (int) : sizeof (long);
  int bytes = (int) (count * size), tag = -1;
  memcpy (out, in, (size_t) bytes);
  char * buf = malloc ((size_t) bytes);
  MPI_Request q;
  if (me == 0)
  {
    for (int r = 1; r < np; r++)
    {
      MPI_Irecv (buf, bytes, MPI_BYTE, r, tag, comm, &q);
      MPI_Waitall (1, &q, MPI_STATUSES_IGNORE);
      for (int k = 0; k < count; k++)
        if (type == MPI_INT)
          ((int *) out) [k] += ((int *) buf) [k];
        else
          ((long *) out) [k] += ((long *) buf) [k];
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
