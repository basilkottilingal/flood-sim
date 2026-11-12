#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scalar.h"

#define nn __CWRAP_SCALAR_WORDS

/* bit s % 64 of word s / 64 is set while the id s is free: see scalar.h */
unsigned long long __cwrap_scalar_free_bits [nn];
#define bits __cwrap_scalar_free_bits

/* the lowest free id, now taken, or -1 */
static inline int popbit ()
{
  for (unsigned i = 0; i < nn; i++)
  {
    unsigned long long v = bits [i];
    if (v != 0)
    {
      bits [i] = v & (v - 1);
      return (int) (i * 64 + __cwrap_scalar_ctz (v));
    }
  }
  return -1;
}

static inline void pushbit (int b)
{
  bits [b/64] |= 1ull << (b & 63);
}

/*
.. Names: a pointer each, never a copy, "" for none (scalar.h); and the order
.. in which the ids were taken, so that scalar_find () gives the last made of
.. a name.
*/
static const char * names [NSCALAR_MAX];
static unsigned long long born [NSCALAR_MAX], births;

scalar scalar_new ()
{
  int s = popbit ();
  if (s < 0)
  {
    fprintf (stderr, "scalar_new: all %d scalars are in use (NSCALAR_MAX)\n",
             NSCALAR_MAX);
    abort ();
  }
  names [s] = "";
  born [s] = ++births;
  if (scalar_new_hook)
    scalar_new_hook ((scalar) s);
  return (scalar) s;
}

scalar __cwrap_scalar_new_named (const char * name)
{
  scalar s = scalar_new ();
  names [s] = name ? name : "";
  return s;
}

const char * scalar_name (scalar s)
{
  return scalar_in_use (s) && names [s] ? names [s] : "";
}

scalar scalar_find (const char * name)
{
  int found = -1;
  if (!name || !*name)
    return -1;
  for (int s = 0; s < NSCALAR_MAX; s++)
    if (scalar_in_use (s) && names [s] && strcmp (names [s], name) == 0
        && (found < 0 || born [s] > born [found]))
      found = s;
  return (scalar) found;
}

int scalar_in_use (scalar s)
{
  return s >= 0 && s < NSCALAR_MAX
         && !((bits [s/64] >> (s & 63)) & 1);
}

void (* scalar_free_hook) (scalar) = NULL;
void (* scalar_new_hook) (scalar) = NULL;

void scalar_free (scalar s)
{
  /*
  .. A global scalar (at file scope, or static in a function) is numbered
  .. once for the whole run, and never freed: another scalar's operators may
  .. depend on it (block-ghosts.h). cwrap refuses such a free where it sees
  .. one; this catches the others, as through a pointer.
  */
  if (s >= 0 && s < n_global_scalars)
  {
    fprintf (stderr, "scalar_free: scalar %s (%d) is global, numbered once "
             "for the whole run: it is never freed\n",
             names [s] && *names [s] ? names [s] : "#unnamed", (int) s);
    abort ();
  }
  /*
  .. an id not in use is left alone: freed twice, or never set, as when a goto
  .. jumps past the declaration of a local scalar.
  */
  if (scalar_in_use (s))
  {
    if (scalar_free_hook)
      scalar_free_hook (s);
    names [s] = "";
    pushbit ((int) s);
  }
}

int scalar_count ()
{
  int n = 0;
  for (int w = 0; w < nn; w++)
    for (unsigned long long v = __cwrap_scalar_used (w); v; v &= v - 1)
      n++;
  return n;
}

/*
.. cwrap puts a slot for each scalar at file scope (and each static scalar in
.. a function) in the section "scalars". The linker joins the slots of every
.. file between __start_scalars and __stop_scalars, which are weak for a
.. program with none. Priority 101 runs this before other constructors, so
.. every id is free before any code can ask for one.
*/
extern const struct __cwrap_scalar_slot __start_scalars [] __attribute__ ((weak));
extern const struct __cwrap_scalar_slot __stop_scalars [] __attribute__ ((weak));
int n_global_scalars = 0;

__attribute__ ((constructor (101))) static void number_global_scalars ()
{
  for (int s = 0; s < NSCALAR_MAX; s++)       /* all ids free */
    pushbit (s);
  for (const struct __cwrap_scalar_slot * p = __start_scalars; p < __stop_scalars; p++)
    *p->s = __cwrap_scalar_new_named (p->name);
  n_global_scalars = (int) (__stop_scalars - __start_scalars);
}
