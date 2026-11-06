/* Nested @loop: each case against the same computation in plain loops */
#include <stdio.h>
#include "loop_nest.h"

static int checks, fails;
#define CHECK(c)                                                        \
  do {                                                                  \
    checks++;                                                           \
    if (!(c)) {                                                         \
      fails++;                                                          \
      printf ("FAIL line %d: %s\n", __LINE__, #c);                      \
    }                                                                   \
  } while (0)

static int side;
static int counted_arg (int v) { side++; return v; }

int main (void)
{
  /* every name the loops declare, as the caller's own: none may change */
  int i = -1, k = -2, count = -3, r = -4, c = -5, s = -6, before = -7;
  int idx = -8, a = -9, n = -10, sq = -11;

  /* 1. a loop in its own block: the inner k hides the outer */
  {
    int got = 0, want = 0;
    repeat (3) {
      int outer = k;
      repeat (4)
        got += outer * 10 + k;
    }
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 4; y++)
        want += x * 10 + y;
    CHECK (got == want);
  }

  /* 2. five deep, as single statements */
  {
    int got = 0;
    repeat (3) repeat (3) repeat (3) repeat (3) repeat (3) got++;
    CHECK (got == 243);
  }

  /* 3. an argument that reads the outer loop's i: evaluated before the
   * inner i exists */
  {
    int got = 0, want = 0;
    range (0, 5)
      range (0, i)
        got += i;
    for (int x = 0; x < 5; x++)
      for (int y = 0; y < x; y++)
        want += y;
    CHECK (got == want);
  }

  /* 4. each argument evaluated once per use */
  {
    int got = 0;
    side = 0;
    range (counted_arg (0), counted_arg (3))
      range (counted_arg (0), counted_arg (2))
        got++;
    CHECK (got == 6 && side == 2 + 3 * 2);
  }

  /* 5. a loop that nests another loop in itself in its body */
  {
    int got = 0, want = 0;
    grid (3, 4)
      got += r * 100 + c;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 4; y++)
        want += x * 100 + y;
    CHECK (got == want);
  }

  /* 6. three levels of definitions */
  {
    long got = 0, want = 0;
    cube (3)
      got += s * 10000 + r * 100 + c;
    for (int z = 0; z < 3; z++)
      for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++)
          want += z * 10000 + x * 100 + y;
    CHECK (got == want);
  }

  /* 7. a name passed on from two levels down */
  {
    int got = 0, want = 0;
    cube_r (3)
      got += r;
    for (int z = 0; z < 3; z++)
      for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++)
          want += x;
    CHECK (got == want);
  }

  /* 8. inner loops before, around and after {...} */
  {
    int got = 0;
    pre_post (4)
      got += k + 1;
    CHECK (got == 1);
    got = 0;
    pre_post (3)
      got += 10;
    CHECK (got == 10);
  }

  /* 9. a break in an inner loop of the body stops that loop only */
  {
    int got = 0, want = 0;
    pairs (4)
      got += a * 10 + b;
    for (int x = 0; x < 4; x++)
      for (int y = 0; y <= x; y++)
        want += x * 10 + y;
    CHECK (got == want);
  }

  /* 10. break and continue in nested blocks of uses */
  {
    int got = 0, want = 0;
    range (0, 5) {
      if (i == 3)
        continue;
      range (0, 5) {
        if (i == 2)
          break;
        got++;
      }
      if (i == 4)
        break;
      got += 100;
    }
    for (int x = 0; x < 5; x++) {
      if (x == 3)
        continue;
      for (int y = 0; y < 5; y++) {
        if (y == 2)
          break;
        want++;
      }
      if (x == 4)
        break;
      want += 100;
    }
    CHECK (got == want);
  }

  /* 11. a parameter the block sees, nested in itself */
  {
    Bag outer = { { 1, 2, 3 }, 3 }, inner = { { 10, 20 }, 2 };
    const Bag * bag = &outer;
    int got = 0, again = 0;
    each (bag) {
      int ox = x;
      each (&inner)
        got += ox * x;
    }
    CHECK (got == (1 + 2 + 3) * (10 + 20));
    /* the inner argument is the outer loop's bag, by the same name */
    each (bag) {
      (void) x;               /* the outer x, which the inner one hides */
      each (bag) {
        again += bag->n;
        (void) x;
      }
    }
    CHECK (again == 3 * 3 * 3);
    CHECK (bag == &outer);
  }

  /* 12. the caller's names are seen in nested blocks, and never changed */
  {
    int got = 0;
    cube (2) {
      got += i + count + idx;
      (void) r, (void) c, (void) s;
    }
    CHECK (got == 8 * (-1 - 3 - 8));
    grid (2, 2) {
      got += k + before;
      (void) r, (void) c;
    }
    CHECK (got == 8 * (-1 - 3 - 8) + 4 * (-2 - 7));
  }

  /* 13. different loops nested as single statements */
  {
    Bag b6 = { { 1, 2, 3, 4, 5, 6 }, 6 };
    int got = 0;
    repeat (2) range (0, 3) each (&b6) got += x;
    CHECK (got == 2 * 3 * 21);
  }

  /* 14. a static local in a body is one per use */
  {
    int c1 = 0, c2 = 0;
    repeat (3)
      counted (1)
        c1 = calls;
    counted (1)
      c2 = calls;
    CHECK (c1 == 3 && c2 == 1);
  }

  /* 15. caches and foreach, nested */
  {
    int got = 0;
    squares (3) {
      int o = v;
      squares (2)
        got += o * 10 + v;
    }
    CHECK (got == 2 * (0 + 10 + 40) + 3 * (0 + 1));
  }

  /* 16. scalars at every level, freed however each level ends */
  {
    int base = scalar_count (), deepest = 0, rounds = 0;
    with_scalar (2) {
      with_scalar (3) {
        scalar t;
        deepest = scalar_count () - base;
        rounds++;
        if (k == 1)
          break;
        (void) t;
        (void) s;
      }
    }
    CHECK (deepest == 3 && rounds == 2 * 2 && scalar_count () == base);
  }

  /* 17. a use inside a statement expression, as the argument of a use */
  {
    int got = 0;
    range (0, ({ int t = 0; repeat (3) t += 2; t; }))
      got++;
    CHECK (got == 6);
  }

  /* 18. many uses in one function, each with hidden names of its own */
  {
    int got = 0;
    range (0, 2) got++;
    range (0, 2) range (0, 2) got++;
    range (0, 2) range (0, 2) range (0, 2) got++;
    CHECK (got == 2 + 4 + 8);
  }

  CHECK (i == -1 && k == -2 && count == -3 && r == -4 && c == -5 && s == -6);
  CHECK (before == -7 && idx == -8 && a == -9 && n == -10 && sq == -11);
  printf ("%d checks, %d failed\n", checks, fails);
  return fails != 0;
}
