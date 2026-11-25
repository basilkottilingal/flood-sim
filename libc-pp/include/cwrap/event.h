/*
.. Events: code that runs at given steps or times, in an order the program
.. states. cwrap turns
..
..   @event name (schedule, position) { body }   a definition
..   @event name (position) { body }             an overload, no schedule
..   @event name;                                a declaration
..
.. into the descriptors below, which every file puts in the linker section
.. cwrap_events. events_run () gathers them from every file of the program,
.. sorts them once, and runs them. See libc-pp/README.md, "Events".
..
.. A name is a slot with one schedule, written on its definition only. Its
.. bodies are the definition's and the overloads'. Three names are reserved,
.. with fixed schedules: defaults, init (once each, before the first step,
.. in that order) and end (once, after the last step, also after a failure).
..
.. A body returns 0, or a failure code that stops the run. It reads the step
.. i and the time t as its parameters.
..
.. The runtime, in event.c, linked into every program by cwrap:
..
..   events_run ()            defaults, init, the steps, end; the status
..   events_dtnext (dt)       the step's dt, cut to land on the next time
..   events_stop ()           ends the run after this step, status 0
..   event_cancel ()          the running body runs no more
..   event_status ()          0, or the failure code (for end's bodies)
..   events_print (fp)        the order, as resolved
..   events_set_reduce (max)  combines a status over the processes (MPI)
..   events_i, events_t, events_dt   the step, the time, the last step's dt
*/
#ifndef CWRAP_EVENT_H
#define CWRAP_EVENT_H

#include <stdio.h>

  /*
  .. A slot's schedule, as numbers: filled by the slot's plan () when the
  .. steps start, after init, so that it may read what init has set.
  */
  enum {
    CWRAP_EVENT_EVERY = 0,      /* no i or t schedule: every step (when ()) */
    CWRAP_EVENT_STEP  = 1,      /* i = i0, i += di */
    CWRAP_EVENT_TIME  = 2,      /* t = t0, t += dt */
  };
  typedef struct {
    int kind;
    long i0, di;                /* di 0: once, at i0 */
    double t0, dt;              /* dt 0: once, at t0 */
    int has_imax, has_tmax, tmax_strict;
    long imax;                  /* i <= imax (i < n is i <= n - 1) */
    double tmax;                /* t <= tmax, or t < tmax if strict */
  } cwrap_event_plan;

  /* where a slot or a body is: first, nothing said, last */
  enum { CWRAP_EVENT_FIRST = 0, CWRAP_EVENT_MIDDLE = 1, CWRAP_EVENT_LAST = 2 };
  /* the reserved slots, in the order they run */
  enum {
    CWRAP_EVENT_STEPS = 0, CWRAP_EVENT_DEFAULTS = 1, CWRAP_EVENT_INIT = 2,
    CWRAP_EVENT_END = 3,
  };

  typedef struct cwrap_event_slot {
    const char *name, *schedule;    /* the schedule as written */
    const char *file;               /* where it is defined */
    int line;
    const char *unit;               /* the file compiled (translation unit) */
    int seq;                        /* its order among that unit's events */
    int phase;                      /* CWRAP_EVENT_STEPS, or reserved */
    int rank;                       /* first, middle, last */
    void (*plan) (cwrap_event_plan *p);
    int (*when) (int i, double t);  /* an extra condition, or NULL */
    const char *const *next;        /* names it runs after, NULL-ended */
    const char *const *prev;        /* names it runs before */
  } cwrap_event_slot;

  typedef struct cwrap_event_body {
    const cwrap_event_slot *slot;
    int (*run) (int i, double t);
    const char *file;
    int line;
    const char *unit;
    int seq;
    int rank;                       /* first, middle, last, in its slot */
  } cwrap_event_body;

  /* the reserved slots, defined in event.c */
  extern const cwrap_event_slot __cwrap_event_slot_defaults,
    __cwrap_event_slot_init, __cwrap_event_slot_end;

  int events_run (void);
  double events_dtnext (double dt);
  void events_stop (void);
  void event_cancel (void);
  int event_status (void);
  void events_print (FILE *fp);
  void events_set_reduce (int (*max) (int status));

  extern int events_i;              /* the step, from 0 */
  extern double events_t;           /* the time, from 0 */
  extern double events_dt;          /* the last dt given by events_dtnext () */

#endif
