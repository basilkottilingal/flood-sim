@loop none (int n) { for (;;) break; }
@loop twice (int n) { {...} {...} }
void outside (void) { {...} }
@loop returns (int n) { if (n) return; {...} }
@loop unnamed (int) { {...} }
@loop varargs (int n, ...) { {...} }
@loop hidden (int n) : z { { int z; } {...} }
@loop ok (int n) { {...} }
@loop ok (long n) { {...} }
@loop self (int n) { self (n) { } {...} }
void calls (void) { ok (1, 2) { } ok () { } }
