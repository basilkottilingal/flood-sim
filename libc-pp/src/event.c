/*
.. The runtime of @event: see include/cwrap/event.h.
..
.. Every body cwrap compiles puts a pointer to its descriptor in the linker
.. section cwrap_events; the linker gathers them from every file between
.. __start_cwrap_events and __stop_cwrap_events. events_run () then
..
..   1. groups the bodies by slot, and checks them: a body linked from two
..      files, a name in next. or prev. that is no event's or is reserved;
..   2. sorts the slots (Kahn): next.X and prev.X are edges, then the rank
..      (first, nothing, last), then the order in their file. Two slots the
..      edges leave unordered, at the same rank, from two files: an error,
..      as is a cycle;
..   3. sorts each slot's bodies: by rank, then those of the slot's own file
..      in their order there (the solver's first), then those of one other
..      file in theirs. Bodies of two other files at one rank: an error;
..   4. runs defaults, init, the steps while a slot may still run, then end,
..      which runs after a failure too.
..
.. No libm: cwrap does not add -lm to a program.
*/
#include <stdlib.h>
#include <string.h>

#include "../include/cwrap/event.h"

const cwrap_event_slot __cwrap_event_slot_defaults = {
  .name = "defaults", .schedule = "", .file = "", .unit = NULL,
  .phase = CWRAP_EVENT_DEFAULTS, .rank = CWRAP_EVENT_MIDDLE,
};
const cwrap_event_slot __cwrap_event_slot_init = {
  .name = "init", .schedule = "", .file = "", .unit = NULL,
  .phase = CWRAP_EVENT_INIT, .rank = CWRAP_EVENT_MIDDLE,
};
const cwrap_event_slot __cwrap_event_slot_end = {
  .name = "end", .schedule = "", .file = "", .unit = NULL,
  .phase = CWRAP_EVENT_END, .rank = CWRAP_EVENT_MIDDLE,
};

/* weak: a program with no event has no such section */
extern const cwrap_event_body *const __start_cwrap_events []
  __attribute__ ((weak));
extern const cwrap_event_body *const __stop_cwrap_events []
  __attribute__ ((weak));

int events_i;
double events_t, events_dt;

struct slot {
  const cwrap_event_slot *s;
  const cwrap_event_body **bodies;  /* in the order they run */
  int nbodies, *cancelled, ncancelled;
  cwrap_event_plan plan;
  long next_i, k;                   /* the next step; the next time's index */
  int done;
};

static struct {
  struct slot *slots;               /* as gathered */
  struct slot **order;              /* the steps' slots, sorted */
  int nslots, nsteps;
  struct slot *phase [4];            /* by CWRAP_EVENT_DEFAULTS ... END */
  int status, stop, cancel, dt_given;
  double t_new;
  int (*reduce) (int);
} ev;

static int same (const char *a, const char *b)
{
  return a == b || (a && b && strcmp (a, b) == 0);
}

static const char *unit_name (const char *u)
{
  return u ? u : "(none)";
}

/* two items of one rank: "both first", "neither first nor last" */
static const char *rank_name (int rank)
{
  return rank == CWRAP_EVENT_FIRST ? "both first"
    : rank == CWRAP_EVENT_LAST ? "both last" : "neither first nor last";
}

static int reserved (const char *name)
{
  return !strcmp (name, "defaults") || !strcmp (name, "init")
    || !strcmp (name, "end");
}

static void where (const cwrap_event_slot *s)
{
  if (s->line)
    fprintf (stderr, "  '%s' is defined at %s:%d\n", s->name, s->file,
             s->line);
}

static struct slot *find_slot (const char *name)
{
  for (int k = 0; k < ev.nslots; k++)
    if (strcmp (ev.slots [k].s->name, name) == 0)
      return &ev.slots [k];
  return NULL;
}

static void release (void)
{
  for (int k = 0; k < ev.nslots; k++) {
    free (ev.slots [k].bodies);
    free (ev.slots [k].cancelled);
  }
  free (ev.slots);
  free (ev.order);
  ev.slots = NULL;
  ev.order = NULL;
  ev.nslots = ev.nsteps = 0;
  memset (ev.phase, 0, sizeof (ev.phase));
}

/* 1. the slots of the bodies linked in, and the three phases */
static int gather (void)
{
  const cwrap_event_body *const *b0 = __start_cwrap_events,
    *const *b1 = __stop_cwrap_events;
  const cwrap_event_slot *phases [4] = {
    NULL, &__cwrap_event_slot_defaults, &__cwrap_event_slot_init,
    &__cwrap_event_slot_end,
  };
  int n = b0 && b1 ? (int) (b1 - b0) : 0, nslots = 0, err = 0;
  struct slot *slots = calloc (n + 3, sizeof (*slots));

  for (int p = CWRAP_EVENT_DEFAULTS; p <= CWRAP_EVENT_END; p++)
    slots [nslots++].s = phases [p];
  for (int j = 0; j < n; j++) {
    const cwrap_event_body *b = b0 [j];
    int k;

    for (int i = 0; i < j; i++)
      if (b0 [i]->line == b->line && same (b0 [i]->file, b->file)) {
        fprintf (stderr, "cwrap: event error: the body of '%s' at %s:%d is "
                 "linked twice, from %s and %s: compile event definitions "
                 "in one file of the program only%s%s%s\n", b->slot->name,
                 b->file, b->line, unit_name (b0 [i]->unit),
                 unit_name (b->unit),
                 b->slot->phase ? "" : ", and declare them (@event ",
                 b->slot->phase ? "" : b->slot->name,
                 b->slot->phase ? "" : ";) in the others");
        err = 1;
      }
    for (k = 0; k < nslots && slots [k].s != b->slot; k++)
      if (strcmp (slots [k].s->name, b->slot->name) == 0) {
        fprintf (stderr, "cwrap: event error: two definitions of '%s'\n",
                 b->slot->name);
        where (slots [k].s);
        where (b->slot);
        err = 1;
      }
    if (k == nslots)
      slots [nslots++].s = b->slot;
    slots [k].nbodies++;
  }
  for (int k = 0; k < nslots; k++) {
    slots [k].bodies = calloc (slots [k].nbodies + 1,
                               sizeof (*slots [k].bodies));
    slots [k].cancelled = calloc (slots [k].nbodies + 1, sizeof (int));
    slots [k].nbodies = 0;
  }
  for (int j = 0; j < n; j++)
    for (int k = 0; k < nslots; k++)
      if (slots [k].s == b0 [j]->slot)
        slots [k].bodies [slots [k].nbodies++] = b0 [j];
  ev.slots = slots;
  ev.nslots = nslots;
  for (int p = CWRAP_EVENT_DEFAULTS; p <= CWRAP_EVENT_END; p++)
    ev.phase [p] = &slots [p - 1];
  return err;
}

/* the order of two bodies of slot s: < 0 if a runs first. Two bodies of
 * two files, neither the slot's, at one rank: by their files' names, which
 * check_bodies () reports */
static int body_order (const cwrap_event_slot *s, const cwrap_event_body *a,
                       const cwrap_event_body *b)
{
  int oa = !same (a->unit, s->unit), ob = !same (b->unit, s->unit);

  if (a->rank != b->rank)
    return a->rank - b->rank;
  if (oa != ob)
    return oa - ob;                 /* the slot's own file first */
  if (!same (a->unit, b->unit))
    return strcmp (unit_name (a->unit), unit_name (b->unit));
  return a->seq - b->seq;
}

/* 3. a slot's bodies, sorted by insertion (there are few), then checked */
static int sort_bodies (struct slot *sl)
{
  const cwrap_event_body **bs = sl->bodies;
  int err = 0;

  for (int j = 1; j < sl->nbodies; j++)
    for (int i = j; i > 0 && body_order (sl->s, bs [i - 1], bs [i]) > 0; i--) {
      const cwrap_event_body *t = bs [i];

      bs [i] = bs [i - 1];
      bs [i - 1] = t;
    }
  for (int j = 0; j < sl->nbodies; j++)
    for (int i = 0; i < j; i++)
      if (bs [i]->rank == bs [j]->rank && !same (bs [i]->unit, bs [j]->unit)
          && !same (bs [i]->unit, sl->s->unit)
          && !same (bs [j]->unit, sl->s->unit)) {
        fprintf (stderr, "cwrap: event error: the order of two bodies of "
                 "'%s' is not defined: %s:%d (in %s) and %s:%d (in %s), "
                 "%s, from two files neither of which defines '%s'. "
                 "Make one of them first or last\n", sl->s->name,
                 bs [i]->file, bs [i]->line, unit_name (bs [i]->unit),
                 bs [j]->file, bs [j]->line, unit_name (bs [j]->unit),
                 rank_name (bs [i]->rank), sl->s->name);
        err = 1;
      }
  return err;
}

/* the steps' slot named in next. or prev. of s, or NULL after an error */
static struct slot *edge (const cwrap_event_slot *s, const char *name,
                          const char *word)
{
  struct slot *t = find_slot (name);

  if (reserved (name)) {
    fprintf (stderr, "cwrap: event error: '%s' is %s.%s, but %s is a phase "
             "of its own, which no step's event is before or after\n",
             s->name, word, name, name);
    where (s);
    return NULL;
  }
  if (!t) {
    fprintf (stderr, "cwrap: event error: '%s' is %s.%s, but no event is "
             "named '%s'; the events are:", s->name, word, name, name);
    for (int k = 0; k < ev.nslots; k++)
      if (ev.slots [k].s->phase == CWRAP_EVENT_STEPS)
        fprintf (stderr, " %s", ev.slots [k].s->name);
    fprintf (stderr, "\n");
    where (s);
    return NULL;
  }
  if (t->s == s) {
    fprintf (stderr, "cwrap: event error: '%s' is %s.%s, itself\n",
             s->name, word, name);
    where (s);
    return NULL;
  }
  return t;
}

/* 2. the steps' slots, by Kahn's algorithm on the edges, choosing among the
 * slots free to go the lowest rank, then the first in its file */
static int sort_slots (void)
{
  int n = 0, err = 0, *ind, (*pred) [2] = NULL, npred = 0, placed = 0;
  struct slot **steps, **order;
  char *done;

  for (int k = 0; k < ev.nslots; k++)
    n += ev.slots [k].s->phase == CWRAP_EVENT_STEPS;
  steps = calloc (n + 1, sizeof (*steps));
  order = calloc (n + 1, sizeof (*order));
  ind = calloc (n + 1, sizeof (*ind));
  done = calloc (n + 1, 1);
  n = 0;
  for (int k = 0; k < ev.nslots; k++)
    if (ev.slots [k].s->phase == CWRAP_EVENT_STEPS)
      steps [n++] = &ev.slots [k];
  /* the edges a -> b, a before b, as indices into steps */
  for (int pass = 0; pass < 2; pass++) {
    npred = 0;
    for (int b = 0; b < n; b++)
      for (int w = 0; w < 2; w++) {
        const char *const *names = w ? steps [b]->s->prev : steps [b]->s->next;

        for (int j = 0; names && names [j]; j++) {
          struct slot *t = edge (steps [b]->s, names [j], w ? "prev" : "next");
          int a = 0;

          if (!t) {
            err |= pass == 0;
            continue;
          }
          while (steps [a] != t)
            a++;
          if (pass) {
            pred [npred][0] = w ? b : a;
            pred [npred][1] = w ? a : b;
          }
          npred++;
        }
      }
    if (pass == 0) {
      if (err)
        break;
      pred = calloc (npred + 1, sizeof (*pred));
    }
  }
  if (err) {
    free (steps);
    free (order);
    free (ind);
    free (done);
    return err;
  }
  for (int e = 0; e < npred; e++)
    ind [pred [e][1]]++;
  while (placed < n) {
    int best = -1;

    for (int a = 0; a < n; a++) {
      const cwrap_event_slot *s, *t;

      if (done [a] || ind [a])
        continue;
      if (best < 0) {
        best = a;
        continue;
      }
      s = steps [a]->s;
      t = steps [best]->s;
      if (s->rank < t->rank
          || (s->rank == t->rank && same (s->unit, t->unit)
              && s->seq < t->seq))
        best = a;
    }
    if (best < 0) {
      /* a cycle: from any slot left, back along its edges until a slot
       * comes again */
      int a = 0, *seen = calloc (n + 1, sizeof (int)), step = 1;

      while (done [a])
        a++;
      while (!seen [a]) {
        seen [a] = step++;
        for (int e = 0; e < npred; e++)
          if (pred [e][1] == a && !done [pred [e][0]]) {
            a = pred [e][0];
            break;
          }
      }
      fprintf (stderr, "cwrap: event error: the order of events is a "
               "cycle:");
      for (int k = 0, b = a; k <= n; k++) {
        fprintf (stderr, " %s%s", steps [b]->s->name, k && b == a ? "" : " ->");
        if (k && b == a)
          break;
        for (int e = 0; e < npred; e++)
          if (pred [e][0] == b && !done [pred [e][1]] && seen [pred [e][1]]) {
            b = pred [e][1];
            break;
          }
      }
      fprintf (stderr, "\n");
      for (int b = 0; b < n; b++)
        if (seen [b])
          where (steps [b]->s);
      free (seen);
      err = 1;
      break;
    }
    /* another slot free to go, of the same rank, from another file */
    for (int a = 0; a < n; a++) {
      const cwrap_event_slot *s = steps [a]->s, *t = steps [best]->s;

      if (a != best && !done [a] && !ind [a] && s->rank == t->rank
          && !same (s->unit, t->unit)) {
        fprintf (stderr, "cwrap: event error: the order of '%s' and '%s' "
                 "is not defined: %s, nothing orders them, and they "
                 "come from two files (%s, %s). Give one of them next. or "
                 "prev. the other, or first or last\n", t->name, s->name,
                 rank_name (s->rank), unit_name (t->unit),
                 unit_name (s->unit));
        where (t);
        where (s);
        err = 1;
      }
    }
    done [best] = 1;
    order [placed++] = steps [best];
    for (int e = 0; e < npred; e++)
      if (pred [e][0] == best)
        ind [pred [e][1]]--;
  }
  free (steps);
  free (ind);
  free (done);
  free (pred);
  ev.order = order;
  ev.nsteps = placed;
  return err;
}

static int resolve (void)
{
  int err;

  release ();
  if (gather ())
    return 1;                       /* what follows would only repeat it */
  err = sort_slots ();
  for (int k = 0; k < ev.nslots; k++)
    err |= sort_bodies (&ev.slots [k]);
  return err;
}

void events_print (FILE *fp)
{
  const char *names [4] = { "", "defaults", "init", "end" };

  if (resolve ()) {
    fprintf (fp, "events: not resolved, see the errors\n");
    return;
  }
  for (int p = CWRAP_EVENT_DEFAULTS; p <= CWRAP_EVENT_END; p++) {
    struct slot *sl = ev.phase [p];

    if (p == CWRAP_EVENT_END) {
      fprintf (fp, "steps\n");
      for (int k = 0; k < ev.nsteps; k++) {
        struct slot *s = ev.order [k];

        fprintf (fp, "  %s (%s) %s:%d\n", s->s->name, s->s->schedule,
                 s->s->file, s->s->line);
        for (int j = 0; j < s->nbodies; j++)
          fprintf (fp, "    %s:%d%s\n", s->bodies [j]->file,
                   s->bodies [j]->line,
                   s->bodies [j]->rank == CWRAP_EVENT_FIRST ? " (first)"
                   : s->bodies [j]->rank == CWRAP_EVENT_LAST ? " (last)" : "");
      }
    }
    fprintf (fp, "%s\n", names [p]);
    for (int j = 0; j < sl->nbodies; j++)
      fprintf (fp, "    %s:%d%s\n", sl->bodies [j]->file, sl->bodies [j]->line,
               sl->bodies [j]->rank == CWRAP_EVENT_FIRST ? " (first)"
               : sl->bodies [j]->rank == CWRAP_EVENT_LAST ? " (last)" : "");
  }
}

/* Times: compared with a tolerance relative to their size. */
static double tol (double x)
{
  return 1e-10 * (x < 0 ? (-x > 1 ? -x : 1) : (x > 1 ? x : 1));
}

static double slot_time (const struct slot *sl, long k)
{
  return sl->plan.t0 + k * sl->plan.dt;
}

/* whether a slot's bounds are past, at i, t */
static int past (const struct slot *sl, int i, double t)
{
  const cwrap_event_plan *p = &sl->plan;

  if (p->has_imax && i > p->imax)
    return 1;
  if (p->has_tmax && (p->tmax_strict ? t >= p->tmax - tol (p->tmax)
                      : t > p->tmax + tol (p->tmax)))
    return 1;
  if (p->kind == CWRAP_EVENT_STEP && p->has_imax && sl->next_i > p->imax)
    return 1;
  if (p->kind == CWRAP_EVENT_TIME && p->has_tmax) {
    double tk = slot_time (sl, sl->k);

    if (p->tmax_strict ? tk >= p->tmax - tol (p->tmax)
        : tk > p->tmax + tol (p->tmax))
      return 1;
  }
  return 0;
}

static int alive (struct slot *sl)
{
  if (!sl->done && (sl->ncancelled == sl->nbodies
                    || past (sl, events_i, events_t)))
    sl->done = 1;
  return !sl->done;
}

/* whether a slot ends of itself: a bound, or a schedule that runs once */
static int finite_slot (const struct slot *sl)
{
  const cwrap_event_plan *p = &sl->plan;

  return p->has_imax || p->has_tmax
    || (p->kind == CWRAP_EVENT_STEP && p->di <= 0)
    || (p->kind == CWRAP_EVENT_TIME && p->dt <= 0);
}

/* The next time a slot or a bound waits for, after t: the time events'
 * next times, and the bounds t <= T still ahead. */
static double tnext (double t)
{
  double tn = __builtin_inf ();

  for (int k = 0; k < ev.nsteps; k++) {
    struct slot *sl = ev.order [k];
    const cwrap_event_plan *p = &sl->plan;

    if (sl->done)
      continue;
    if (p->kind == CWRAP_EVENT_TIME) {
      long j = sl->k;
      double tj = slot_time (sl, j);

      while (tj <= t + tol (tj) && p->dt > 0)
        tj = slot_time (sl, ++j);
      if (tj > t + tol (tj) && tj < tn
          && !(p->has_tmax && tj > p->tmax + tol (p->tmax)))
        tn = tj;
    }
    if (p->has_tmax && p->tmax > t + tol (p->tmax) && p->tmax < tn)
      tn = p->tmax;
  }
  return tn;
}

double events_dtnext (double dt)
{
  double t = events_t, tn = tnext (t);

  if (dt > 0 && tn < __builtin_inf ()) {
    double rem = tn - t;
    long n = (long) (rem / dt);

    if (n * dt < rem * (1 - 1e-9))
      n++;                          /* the fewest steps no longer than dt */
    if (n <= 1) {
      dt = rem;
      ev.t_new = tn;                /* lands on tn exactly */
    } else {
      dt = rem / n;
      ev.t_new = t + dt;
    }
  } else
    ev.t_new = t + dt;
  ev.dt_given = 1;
  events_dt = dt;
  return dt;
}

void events_stop (void)
{
  ev.stop = 1;
}

void event_cancel (void)
{
  ev.cancel = 1;
}

int event_status (void)
{
  return ev.status;
}

void events_set_reduce (int (*max) (int))
{
  ev.reduce = max;
}

/* runs the bodies of a slot not cancelled, until one fails (all of them if
 * keep_going); the status, the first failure kept as the run's */
static int run_slot (struct slot *sl, int keep_going)
{

  for (int j = 0; j < sl->nbodies; j++) {
    const cwrap_event_body *b = sl->bodies [j];
    int rc;

    if (sl->cancelled [j])
      continue;
    ev.cancel = 0;
    rc = b->run (events_i, events_t);
    if (rc)
      fprintf (stderr, "event '%s' (%s:%d) failed with %d at i = %d, "
               "t = %g\n", sl->s->name, b->file, b->line, rc, events_i,
               events_t);
    if (ev.reduce)
      rc = ev.reduce (rc);          /* every process learns of a failure */
    if (ev.cancel) {
      sl->cancelled [j] = 1;
      sl->ncancelled++;
    }
    if (rc && !ev.status)
      ev.status = rc;
    if (ev.status && !keep_going)
      break;
  }
  return ev.status;
}

/* whether a step's slot is due at i, t; if so, its schedule moves on */
static int due (struct slot *sl, int i, double t)
{
  const cwrap_event_plan *p = &sl->plan;
  int now = 0;

  if (!alive (sl))
    return 0;
  if (p->kind == CWRAP_EVENT_STEP) {
    if (i >= sl->next_i) {
      now = 1;
      if (p->di > 0)
        while (sl->next_i <= i)
          sl->next_i += p->di;
      else
        sl->done = 1;
    }
  } else if (p->kind == CWRAP_EVENT_TIME) {
    double tk = slot_time (sl, sl->k);

    if (t >= tk - tol (tk)) {
      now = 1;
      if (p->dt > 0)
        while (slot_time (sl, sl->k) <= t + tol (slot_time (sl, sl->k)))
          sl->k++;
      else
        sl->done = 1;
    }
  } else
    now = 1;
  return now && (!sl->s->when || sl->s->when (i, t));
}

int events_run (void)
{
  ev.status = ev.stop = 0;
  events_i = 0;
  events_t = events_dt = 0;
  if (resolve ()) {
    release ();
    return -1;
  }
  ev.status = run_slot (ev.phase [CWRAP_EVENT_DEFAULTS], 0);
  if (!ev.status)
    ev.status = run_slot (ev.phase [CWRAP_EVENT_INIT], 0);
  /* the schedules, now that init has run */
  for (int k = 0; k < ev.nsteps; k++) {
    struct slot *sl = ev.order [k];

    memset (&sl->plan, 0, sizeof (sl->plan));
    if (sl->s->plan)
      sl->s->plan (&sl->plan);
    sl->next_i = sl->plan.i0;
    sl->k = 0;
    sl->done = 0;
  }
  while (!ev.status && !ev.stop) {
    int any = 0, finite = 0, finite_alive = 0;

    /* the run goes on while a finite slot (bounded, or once) may run; with
     * none in the program, while any may, until events_stop () */
    for (int k = 0; k < ev.nsteps; k++) {
      int f = finite_slot (ev.order [k]), a = alive (ev.order [k]);

      any |= a;
      finite |= f;
      finite_alive |= f && a;
    }
    if (finite ? !finite_alive : !any)
      break;
    ev.dt_given = 0;
    for (int k = 0; k < ev.nsteps && !ev.status; k++)
      if (due (ev.order [k], events_i, events_t))
        ev.status = run_slot (ev.order [k], 0);
    if (ev.status)
      break;
    if (ev.dt_given)
      events_t = ev.t_new;
    else {
      double tn = tnext (events_t);

      if (tn < __builtin_inf ())
        events_t = tn;
    }
    events_i++;
  }
  run_slot (ev.phase [CWRAP_EVENT_END], 1);
  release ();
  return ev.status;
}
