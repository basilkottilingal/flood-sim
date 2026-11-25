#include <cwrap/event.h>
@event vof (i++) {}
@event vof (i += 2) {}
@event gauge () {}
@event init (i = 0) {}
@event end (next.vof) {}
@event vof (next.adapt) {}
@event both (i++, t += 1) {}
@event tpp (t++) {}
@event minus (i -= 1) {}
@event fl (i++, first, last) {}
@event twice (i = 1, i = 2) {}
@event init;
@event ret (i++) { return; }
@event ok (i++, i % 2 == 0) { return 0; }
