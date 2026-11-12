#include <stdio.h>

scalar h, zb;
vector u;
static const scalar * const global_h = &h;   /* h, where a local hides it */

static void show (scalar s)
{
  const char * n = scalar_name (s);
  printf ("%s ", *n ? n : "-");
}

/* a tmp of its own while inner's is alive: the same name, the last found */
static void second (void)
{
  scalar tmp;
  printf ("%d ", scalar_find ("tmp") == tmp);
}

static void inner (void)
{
  scalar h;                         /* hides the global h */
  show (h);
  {
    scalar h;                       /* hides inner.h */
    show (h);
    printf ("%d %d ", scalar_find ("inner.h.2") == h, scalar_find ("h") == *global_h);
  }
  static scalar kept;               /* static: its function's, always */
  show (kept);
  scalar tmp;                       /* hides no scalar: its own name */
  show (tmp);
  second ();
  vector u;                         /* hides the global u */
  show (u.x);
  show (u.y);
}

int main (void)
{
  show (h);
  show (zb);
  show (u.x);
  show (u.y);
  printf ("| ");
  inner ();
  printf ("| ");
  scalar m = scalar_new ();         /* not cwrap's: named "" */
  show (m);
  printf ("%d ", scalar_find (""));
  scalar_free (m);
  printf ("%s\n", *scalar_name (m) ? "?" : "-");
  return 0;
}
