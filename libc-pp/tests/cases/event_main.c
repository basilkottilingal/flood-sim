/* The program: the solver's header, an overload of adv defined elsewhere,
 * an event after it. */
#include "event_solver.h"

char log_ [4096];
int n_;

@event adv;
@event adv { LOG ("adv.main "); }
@event adv (first) { LOG ("adv.first "); }
@event gauge (i += 2, next.adv) { LOG ("gauge%d ", i); }
@event init (last) { LOG ("init.main "); }

int main (void)
{
  int rc = events_run ();

  printf ("%s\nrc %d\n", log_, rc);
  return 0;
}
