/* Constructs that are easy to get wrong in a C grammar. The test passes when
 * cwrap reproduces this file byte for byte and gcc accepts the result, so each
 * line is valid GNU C. */
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>

/* Typedef names and scopes */
typedef int T;
typedef T U, *P, A[3], (*Fn) (T);
typedef int T;                          /* a repeated typedef is allowed */
typedef struct node node;
struct node { node *next; T T; int U : 3; int : 2; };
static T shadow_block (void) { int T = 1; T = T * 2; return T; }
static int shadow_param (int T) { return T * 2; }
static int shadow_nested_fn (int (*T) (void)) { return T (); }
static int (*fn_returning (int T, int x)) (int) { (void)x; (void)T; T = 0; return NULL; }
static void uses_type_after_scope (void) { { int T; (void)T; } T x = 0; (void)x; }
static int type_in_parens (int (T)) { return 0; }
static int after_param_list (int a, T b) { T c = a + b; return c; }
static int for_scope (void)
{
  int s = 0;
  for (int T = 0; T < 3; T++)
    s += T;
  T r = s;
  return r;
}
/* foreach is a keyword only until it is declared, and in is always a name */
struct foreach { int foreach; };
static int foreach_member (struct foreach *in) { return in->foreach; }
static int foreach_param (int foreach) { return foreach * 2; }
static int foreach (int in) { return in; }
static int foreach_call (int (*f) (int))
{
  int (*foreach_ptr) (int) = foreach;
  return foreach (1) + foreach (foreach (2)) + (foreach) (3) + (*foreach) (4)
    + foreach_ptr (5) + f (6) + (foreach == f);
}
static int foreach_var (int foreach)
{
  int a[2] = { 0 }, *p = &foreach;
  foreach = 1;
  foreach++;
  a[foreach - 1] = foreach ? foreach : -foreach;
  goto foreach;
foreach:
  return *p + (int) sizeof foreach + (foreach) + a[0];
}
static int foreach_enum (void) { enum { foreach = 2 }; return foreach * foreach; }
static int in_is_a_type (void) { typedef int in; in s = 0; return s; }
static T *pointer_then_name (T *U) { return U; }
static void pointer_decl (void) { T *x = NULL; T *y = x; (void)y; }
enum { E0, E1 = 3, E2, };
static void enum_shadow (void) { enum { T = 2 }; int x = T * 3; (void)x; }
static int label_named_type (void) { goto T; T: return 1; }
static int member_named_type (struct node *n) { return n->T + (*n).T; }
static int casts (long v) { return (T)v + (int)(T)(long)v + sizeof (T) + sizeof v; }

/* Declarators */
int *(*complicated[2]) (char (*) [3], ...);
static int vla (int n, int a[static n], int b[const], int c[*]);
static int vla (int n, int a[static n], int b[const], int c[n]) { return a[0] + b[0] + c[0]; }
char *restrict_ptr (const char *restrict s);

/* Expressions */
static int expressions (int a, int b, int *p, struct node *n)
{
  a = b ? a : b;
  a = b ?: a;
  a += a << 1 >> 1 & 3 | 4 ^ 5;
  a = (a, b);
  a = -a + +b - ~a + !b;
  (*p)++, ++*p, p[0]--, --p[0];
  a = (struct node) { .next = n, .T = 1 }.T;
  a = ((int [2]) { 1, 2 })[1];
  a = sizeof (int) + sizeof (int [2]) + sizeof a + _Alignof (int);
  return _Generic (a, int: 1, long: 2, default: 3);
}

/* Initializers */
struct pair { int a, b; };
static struct pair pairs[4] = { [0] = { 1, 2 }, [1 ... 2] = { .a = 3 }, { b: 4 }, };
static int empty[1] = { };
static const char *strings = "a" "b" L"" "c";

/* Statements */
static int statements (int a)
{
  __label__ out;
  if (a)
    if (a > 1)
      a = 1;
    else
      a = 2;
  while (a < 0)
    a++;
  do a--; while (a > 10);
  switch (a) {
  case 0:
    break;
  case 1 ... 3:
    a = 1;
    __attribute__ ((fallthrough));
  default:
    ;
  }
  switch (a) {
  case 5:
    a = 6;
    int after_label = 1;
    (void)after_label;
  default:
  }
  void *target = &&out;
  goto *target;
out:
  return a;
}

/* GNU extensions */
__extension__ typedef long long big;
__extension__ static inline int ext_fn (void) { return 1; }
static int attr_fn (int a __attribute__ ((unused))) __attribute__ ((noinline, const));
static int attr_fn (int a) { return ({ int b = a; b; }); }
struct __attribute__ ((packed)) packed { char c; __extension__ int i; };
typedef struct { int x; } __attribute__ ((aligned (8))) aligned_t;
static __typeof__ (int) typeof_type = 1;
static typeof (typeof_type) typeof_expr = 2;
extern int renamed (void) __asm__ ("real_name");
static inline int asm_stmt (int a)
{
  __asm__ volatile ("nop" : [out] "=r" (a) : "0" (a) : "cc");
  __asm__ ("nop");
  return a;
}
static int va (int n, ...)
{
  va_list ap;
  va_start (ap, n);
  n += va_arg (ap, int);
  va_end (ap);
  return n + (int) offsetof (struct pair, b)
    + __builtin_types_compatible_p (int, T);
}
_Static_assert (sizeof (int) >= 2, "int is too small");
_Alignas (16) static char aligned_buf[16];
_Noreturn void stop (void);
_Atomic int atomic_q;
_Atomic (int) atomic_s;
_Thread_local int tls;
__int128 wide;
double _Complex z;
;

int main (void)
{
  return 0 * (shadow_block () + shadow_param (1) + for_scope ());
}
