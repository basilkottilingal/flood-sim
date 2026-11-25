/* gcc's errors in an event's body point at its own lines. */
#include <cwrap/event.h>

@event flux (i++,
             i <= 3,
             first)
{
  int a = 1;

  return a + undeclared;
}
