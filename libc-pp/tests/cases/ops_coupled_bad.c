typedef struct { int s; } Pr;
enum { GHOST_PROLONGATION, GHOST_RESTRICTION, GHOST_BOUNDARY };
scalar zb;
@prolongation needs (const Pr * pr) (scalar a) { (void) pr; (void) a; }
@prolongation wrong (const Pr * pr) (int a) { (void) pr; (void) a; }
void f (const Pr * pr, int id, scalar local)
{
  @prolongation (id) (pr);
  ghost_ops (zb, prolongation (needs, local));
  ghost_ops (zb, smoothing (needs));
}
