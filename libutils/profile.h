#ifndef _PROFILE_H_
#define _PROFILE_H_

  /*
  .. Where the time goes: regions of code timed on every thread and every
  .. MPI process, reported as a table, a flame graph or a timeline.
  ..
  .. Marking regions (active with -DPROFILE_ENABLE; without it, nothing):
  ..   PROFILE ();              the rest of this function, named after it
  ..   PROFILE_SCOPE ("name");  the rest of this block, named name
  .. Both are declarations: put them first in a block. The region ends
  .. where the block does, by any way out (return, break, goto), as the
  .. compiler calls profile_leave () there (GCC's cleanup attribute; GCC
  .. and clang). Regions nest: each is timed on its call path, the
  .. regions open around it on this thread. A name may be used at several
  .. places: they are one region. A place keeps the name it was first
  .. given (its id is found once): give each place one name, a constant
  .. or a string that does not change.
  ..
  .. What is measured, a call path a thread: calls, inclusive time (from
  .. entering to leaving), self time (inclusive less the regions inside
  .. it), the longest call. The clock is CLOCK_MONOTONIC; an enter and a
  .. leave cost some 20 to 60 ns, measured at start and reported: time
  .. functions and loops over blocks, not cells. A thread's paths start at
  .. the outermost region it entered itself: in an OpenMP loop, a worker's
  .. regions are roots of their own, not under the master's.
  ..
  .. A process's value of a path: calls summed over its threads, times the
  .. largest of its threads' (the wall time of the path on that process).
  .. Over the MPI processes: mean, largest, smallest, and the imbalance,
  .. largest over mean (1: balanced).
  ..
  .. Calls (collective on MPI_COMM_WORLD with -D_MPI, MPI started; every
  .. process the same calls in the same order):
  ..   profile_snapshot ()   a snapshot of everything so far, for a report
  ..                         of the run as it goes: begun here, finished
  ..                         by the next calls, never waiting on the
  ..                         network (non-blocking collectives on a
  ..                         communicator of its own). Only if the
  ..                         previous one is still not done does it wait
  ..                         for it. Rank 0 appends each snapshot's table
  ..                         to PREFIX_snapshots.txt
  ..   profile_query (path, &calls, &incl, &self, &imbalance)
  ..                         rank 0 (or a process alone), after a snapshot
  ..                         or the end: the means over the processes of
  ..                         a path "a;b;c" (calls, seconds), and its
  ..                         imbalance (largest inclusive time over the
  ..                         mean); 0 if none. Any pointer may be NULL
  ..   profile_finish ()     the end: the last snapshot, waited for, and
  ..                         the outputs. Call it before MPI_Finalize; at
  ..                         exit it is called if not (with MPI still
  ..                         running, every process must exit that way).
  ..                         After MPI_Finalize, each process writes its
  ..                         own table, PREFIX_rankN.txt
  ..
  .. Outputs, by the environment (read once):
  ..   PROFILE      a list of: table, flame, trace; or off (no timing).
  ..                Default: table
  ..   PROFILE_OUT  the files' prefix. Default: profile
  ..   table   the call tree on stdout (rank 0) and in PREFIX.txt: calls a
  ..           process, self and inclusive time (mean, max), imbalance, %
  ..           of the wall time
  ..   flame   PREFIX.folded: "a;b;c us" a path, the mean self time, for
  ..           speedscope.app or flamegraph.pl
  ..   trace   PREFIX.trace.json: every call as a bar, a row a process and
  ..           thread (Chrome's trace events: ui.perfetto.dev or
  ..           chrome://tracing); at most PROFILE_TRACE_EVENTS calls a
  ..           thread (default 65536), the first ones, the rest counted.
  ..           Gathered on rank 0 at the end: for runs of a few processes
  ..
  .. Example:
  ..   void step (void) { PROFILE (); ... }
  ..   for (int n = 0; n < steps; n++) {
  ..     step ();
  ..     if (n % 100 == 0) profile_snapshot ();
  ..   }
  ..   profile_finish ();
  ..   cc -DPROFILE_ENABLE ... profile.c;  PROFILE=table,flame ./run
  */

  /* a region left at the end of its block: by the cleanup attribute */
  typedef struct { int on; } ProfileScope;

  int  profile_enter    (int * id, const char * name);
  void profile_leave    (ProfileScope * scope);
  void profile_snapshot (void);
  void profile_finish   (void);
  int  profile_query    (const char * path, double * calls, double * incl,
                         double * self, double * imbalance);

  #define PROFILE_CAT_(a, b)  a##b
  #define PROFILE_CAT(a, b)   PROFILE_CAT_(a, b)

  #ifdef PROFILE_ENABLE
    /* the place's region id, found once (-1: not yet) */
    #define PROFILE_SCOPE(name)                                        \
      static int PROFILE_CAT (profile_id_, __LINE__) = -1;             \
      ProfileScope PROFILE_CAT (profile_scope_, __LINE__)              \
        __attribute__ ((cleanup (profile_leave), unused))              \
        = { profile_enter (&PROFILE_CAT (profile_id_, __LINE__), (name)) }
  #else
    #define PROFILE_SCOPE(name)  do { } while (0)
  #endif
  #define PROFILE()  PROFILE_SCOPE (__func__)

#endif
