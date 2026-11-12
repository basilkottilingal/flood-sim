typedef struct { int s; } Pr;
enum { GHOST_PROLONGATION, GHOST_RESTRICTION, GHOST_BOUNDARY };
scalar zb;
@prolongation needs (const Pr * pr) (scalar a) { (void) pr; (void) a; }
@prolongation wrong (const Pr * pr) (int a) { (void) pr; (void) a; }
@boundary b_one (int s) { (void) s; }
void f (const Pr * pr, int id, scalar local)
{
  @prolongation (id) (pr);
  ghost_ops (zb, prolongation (needs, local));
  ghost_ops (zb, smoothing (needs));
  ghost_ops (zb, boundary (b_one));
  ghost_ops (zb, boundary.back (b_one));
  ghost_ops (zb, boundary.middle (b_one));
  ghost_ops (zb, boundary.left (b_one), boundary.all (b_one));
}
