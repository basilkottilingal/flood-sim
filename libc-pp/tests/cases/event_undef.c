/* An overload of an event declared, whose definition is linked from no
 * file: the linker's undefined reference. */
#include <cwrap/event.h>

@event adv;
@event adv {}

int main (void)
{
  return events_run ();
}
