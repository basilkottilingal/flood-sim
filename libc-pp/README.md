# c wrapper

*Documentation as of 2026-11-10.*

This subdirectory is written with the help of claude-code

`cwrap` is a front end that works like `gcc`. It accepts gcc's command line,
and for each `.c` file it runs:

    gcc -E   ->   cwrap's translator   ->   gcc -x cpp-output

Then it links the way gcc would. The translator parses the preprocessed C into
a lossless syntax tree: every token keeps its whitespace and linemarkers.
cwrap rewrites the tree and writes it back out. Line numbers are preserved, so
gcc's errors, warnings and debug info still point at your source.

On top of C11 with GNU extensions, cwrap adds:

- **`scalar`, `vector` and `tensor`**, with `-dimension=1`, `2` or `3`:
  integer ids that cwrap hands out, at file scope and in functions, freed
  automatically as their scope ends
- **`cache (T)`**: growable arrays that rebuild themselves on demand
- **`foreach (v in V)`**: a loop over a cache, or over any struct with fields
  `n` and `a`
- **`@loop`**: block macros, written as ordinary code with typed arguments,
  for loops that `#define` cannot express cleanly
- **`foreach ()` and `s[]`, `s[i,j,k]`**: Basilisk-style stencils, every cell
  of every leaf block of libb-amr, with `-dimension`

---

## Contents

1. [Building](#building)
2. [Using cwrap](#using-cwrap)
3. [Diagnostics](#diagnostics)
4. [Scalars](#scalars)
5. [Caches](#caches)
6. [foreach](#foreach)
7. [@loop: block macros](#loop-block-macros)
8. [Stencils: foreach () and s[]](#stencils-foreach--and-s)
9. [Ghost operators: @prolongation, @restriction, @boundary](#ghost-operators-prolongation-restriction-boundary)
10. [How the translator works](#how-the-translator-works)
11. [Source layout](#source-layout)
12. [Tests](#tests)
13. [Limitations](#limitations)
14. [Future: a single build, and portable C99 output](#future-a-single-build-and-portable-c99-output)
15. [Future: any object as the collection in foreach](#future-any-object-as-the-collection-in-foreach)

---

## Building

You need:
- gcc, or a compatible compiler, set with `CWRAP_CC`
- GNU bison and flex
- GNU ld on ELF, which means Linux

ld is needed because global scalars are numbered through a linker section; see
[Scalars](#scalars).

    make            # build/cwrap and build/scalar.o
    make test       # end-to-end tests
    make clean

The Makefile builds the full paths of `src/`, `build/` and `include/` into the
binary, so `build/cwrap` finds its runtime wherever you run it from. If you
move the tree, rebuild with `make clean && make`.

To build a project with cwrap, use it as the compiler:

    make CC=/path/to/libc-pp/build/cwrap

---

## Using cwrap

cwrap takes the same options as gcc:

    cwrap -Wall -O2 -c foo.c -o foo.o
    cwrap foo.c bar.c -o prog -lm
    cwrap foo.o bar.o -o prog          # linking objects alone also works

### -translate

`-translate` stops after the translator and writes C instead of compiling:

    cwrap -translate foo.c -o foo.out.c   # write to a file
    cwrap -translate foo.c                # write to standard output
    cwrap -translate a.c b.c > all.c      # several files, joined in order

- **Flags are still used for preprocessing:** `-D`, `-I` and `-include` apply
  as usual. Flags for compiling and linking have no effect.
- **Several files are translated one by one and joined.** They are not merged
  into one translation unit; see [Future](#future-a-single-build-and-portable-c99-output).
- **The output is preprocessed C** for this platform, with the system headers
  expanded.

### Passed straight to gcc

These run gcc unchanged:
- `-E`, `-M`, `-MM` and `-x`
- reading from `-`
- a command with no `.c` file, such as linking

cwrap still adds its include directory, and its runtime when linking.

### -dimension=N

`-dimension=1`, `-dimension=2` or `-dimension=3` gives the euclidean dimension,
and turns on `scalar`, and for 2 and 3 also `vector` and `tensor`; see
[Scalars](#scalars).
- **Without it,** these are ordinary names that your program may define itself.
- **It is cwrap's alone:** gcc never sees it.
- **Any other value is an error.**
- **Every file of a program must use the same dimension.**

### What cwrap adds to every command

- **`-DCWRAP_DIMENSION=N -include src/scalar.h`** in each preprocessing run,
  with `-dimension=N` only.
- **`-I include/`**, placed after your own `-I` directories, so
  `#include <cwrap/cache.h>` resolves.
- **`build/scalar.o`** in each link.

### Environment variables

| Variable | Effect |
|---|---|
| `CWRAP_CC` | the compiler to run (default `gcc`) |
| `CWRAP_VERBOSE` | print each command cwrap runs |
| `CWRAP_KEEP` | keep the directory of intermediate files, and print its path |
| `CWRAP_DUMP_AST` | print each parsed tree to stderr |

In the kept directory, `N-name.i` is gcc's preprocessed output and
`N-name.t.i` is the translation of it.

---

## Diagnostics

cwrap's messages look like gcc's:

    foo.c:8: error: '=' copies cache 'w', and two caches must not share an array
    foo.c:7: note: 'w' is declared here

### Colors

| Part | Color |
|---|---|
| `file:line:` | bold |
| `error:` | bold red |
| `warning:` | bold magenta |
| `note:` | bold cyan |
| text in `'quotes'` | bold |

- **By default, color is used only on a terminal**, and not when `TERM=dumb`.
- **`-fdiagnostics-color`, `=always`, `=never` and `-fno-diagnostics-color`**
  work as in gcc; the last one given counts. gcc gets them too.
- **`GCC_COLORS`** sets the colors with the keys `error`, `warning`, `note`,
  `locus` and `quote`. If it is set but empty, color is off.

### What happens to a file

| Situation | Result |
|---|---|
| The file parses | translated, then compiled by gcc |
| A syntax error that cwrap's grammar does not cover | `warning: not translated, leaving this file to gcc`. cwrap lists up to 20 syntax errors, then passes the original preprocessed file to gcc, which accepts it or explains the problem in its own words. |
| A cwrap error, such as copying a cache or a `foreach` while `foreach` is declared | reported, and the file never reaches gcc. With `-translate`, nothing is written and the exit status is 1. |
| A cwrap warning, such as changing a scalar or passing a cache by value | reported, and compilation goes on |

---

## Scalars

With `-dimension=N`, cwrap force-includes `src/scalar.h` into every file, with
`CWRAP_DIMENSION` defined as N. It declares these types:

| `-dimension=` | Types |
|---|---|
| 1 | `typedef int scalar;` |
| 2 | that, and `typedef struct { scalar x, y; } vector;` and `typedef struct { vector x, y; } tensor;` |
| 3 | that, and `typedef struct { scalar x, y, z; } vector;` and `typedef struct { vector x, y, z; } tensor;` |

A scalar is an id, from 0 to `NSCALAR_MAX - 1` (64), taken from one pool
shared by the whole program. A vector or tensor is a set of scalars, and
everything below about scalars applies to each of them, in the order
`x, y (, z)`, and for a tensor `x.x, x.y, ..., y.x, ...`. `scalar.o` is linked
into every program, whatever the dimension.

**Where the values are.** In every dimension, 1 included, two more types say
where a scalar's values sit in a cell:

- `face_vector F;`: a struct of scalars `F.x` (`F.y`, `F.z`), each at the
  cell's lower face along its axis (at i - 1/2 for `F.x`): fluxes, face
  velocities;
- `vertex_scalar v;`: a scalar at the cell's lower corner.

They are numbered, named and freed like the others (`face_vector_free ()`;
a global one never). `scalar_location (s)` gives where a scalar is:
`CWRAP_CENTRE`, `CWRAP_FACE_X`, `_Y`, `_Z` or `CWRAP_VERTEX`. The grid
(libb-amr) reads it to treat such scalars as faces or vertices. A slot
carries the location, `{ &F.x, "F.x", 1 }`, and a local face vector gets
`__cwrap_face_vector_new_named ()`, a vertex scalar
`__cwrap_scalar_new_located ("v", 4)`.

### File-scope scalars

    scalar h, u, v;      // ids 0, 1, 2 ... in the order the linker sees them

- **Each one gets a slot,** its address and its name, on its own line, so no
  line numbers move:

      scalar h; static const struct __cwrap_scalar_slot __cwrap_scalar_slot_0
        __attribute__ ((section ("scalars"), used)) = { &h, "h" };

  A vector or tensor gets one slot for each of its scalars. In 2 dimensions,
  `tensor t;` gets four, all on its line:

      static const struct __cwrap_scalar_slot __cwrap_scalar_slot_0 ... = { &t.x.x, "t.x.x" };
      static const struct __cwrap_scalar_slot __cwrap_scalar_slot_1 ... = { &t.x.y, "t.x.y" };
      static const struct __cwrap_scalar_slot __cwrap_scalar_slot_2 ... = { &t.y.x, "t.y.x" };
      static const struct __cwrap_scalar_slot __cwrap_scalar_slot_3 ... = { &t.y.y, "t.y.y" };

- **The linker joins the slots of every file** between `__start_scalars` and
  `__stop_scalars`.
- **A constructor in `scalar.c`, at priority 101, numbers them** before `main`
  and before any other constructor, and names them. So scalars across
  separately compiled files get ids 0 to n-1 with no gaps.
- **`n_global_scalars`** holds n.
- **Not numbered:**
  - `extern` declarations
  - `typedef`s
  - pointers, arrays and functions
  - names with an initializer of their own, which also get a warning

### Scalars in functions

    void f (void)
    {
      scalar x;                 // a new id, named "x", freed as the scope ends
      static scalar s;          // numbered at startup, like a global
      scalar k = scalar_new (); // yours: freed only when you call scalar_free
      vector u = vector_new (); // the same: vector_free (u) when done
    }

- **`scalar x;` becomes**

      scalar x __attribute__ ((cleanup (__cwrap_scalar_free))) = __cwrap_scalar_new_named ("x");

  so `scalar_free (x)` runs however the scope ends: the end of the block,
  `return`, `break`, `continue` or `goto`. Scalars are freed in reverse order.
- **A vector or tensor gets an id for each scalar,** and a cleanup that frees
  them, last first:

      tensor w __attribute__ ((cleanup (__cwrap_tensor_free)))
        = __cwrap_tensor_new_named ((const char *const []) { "w.x.x", "w.x.y", "w.y.x", "w.y.y" });

- **`static scalar s;`** becomes `static scalar s;` followed by one slot for each
  of its scalars, as declarations of their own on the same line. They are
  numbered with the globals, and named after their function: `f.s`.
- **Left alone:** pointers, arrays, names with an initializer, and `extern`,
  `register` and thread-local scalars.

### Changing a scalar

Assigning to a scalar is a **warning**, because cwrap gives scalars their
values. This covers `=`, the compound assignments, and `++` and `--`. It also
covers a whole vector or tensor, or part of one, as in `v = w`, `t.y = v` and
`t.x.y = 1`. The warning names the operator and what it changes, and a note
shows the declaration.

### Names

The scalars cwrap numbers or makes, and only those, are named after their
declaration, so that messages, output and visualisation can say `h` and not
`scalar 7`. A program never names a scalar itself:

    scalar h;                  /* "h" */
    vector u;                  /* "u.x", "u.y" */
    void f (void)
    {
      scalar tmp;              /* "tmp": it hides no scalar */
      scalar h;                /* "f.h": it hides the global h */
      {
        scalar h;              /* "f.h.2": it hides f.h */
      }
      static scalar kept;      /* "f.kept": static, always its function's */
      scalar k = scalar_new ();  /* "": not cwrap's */
    }

- **No allocation, no copy:** a name is a pointer to a string literal that
  cwrap writes, in the slot of a global or in `__cwrap_scalar_new_named
  ("tmp")`; `scalar.c` keeps the pointer.
- **A local that hides another scalar of its name** (a global, or a local of
  an outer block) is named after its function, `f.h`, and a deeper one is
  numbered, `f.h.2`, `f.h.3`, so that the scalars alive in one scope never
  share a name. cwrap sees what each declaration hides when it reads it. In a
  `@loop` definition, the loop's name stands for the function's.
- **Two functions may each have a `tmp` alive at once** (one calls the
  other): names are labels, not keys, and `scalar_find ()` gives the last
  made.
- **API (`src/scalar.h`):**

  | | |
  |---|---|
  | `scalar_name (s)` | its name: `""` if it has none, or is not in use; never `NULL` |
  | `scalar_find (name)` | the scalar in use with that name, the last made if several; -1 if none, or for `""` |

  (`__cwrap_scalar_new_named`, `__cwrap_vector_new_named` and
  `__cwrap_tensor_new_named` are cwrap's, for the locals it makes.)
- **Any other scalar is named `""`:** one from `scalar_new ()`, `vector_new
  ()` or `tensor_new ()`, and any id once freed.
- **libb-amr's messages use them:** `ghost_ops: scalar zbl (1): its operators
  would depend on itself, through scalar hl (0)`.

### A global scalar is never freed

A scalar at file scope, `static` in a function, or `extern`, and each scalar
of such a vector or tensor, is numbered once, before `main ()`, and keeps its
id for the whole run: another scalar's operators may depend on it, as
coupled prolongations will. So:

- **`scalar_free ()`, `vector_free ()` or `tensor_free ()` of one, or of a
  member of one, is a compile error,** with a note at its declaration:

      scalar zb;
      ...
      scalar_free (zb);   /* error: scalar_free () of 'zb', a global scalar: it
                             is numbered once for the whole run and never freed */

- **A free cwrap cannot see,** as through a pointer, stops the run:
  `scalar_free ()` of an id below `n_global_scalars` prints the id and
  aborts.

A local scalar, and one from `scalar_new ()`, is freed as before.

### Runtime API (`src/scalar.h`)

| Function | |
|---|---|
| `scalar scalar_new (void)` | the lowest free id. Aborts with a message if all `NSCALAR_MAX` are in use. |
| `void scalar_free (scalar)` | frees an id. An id not in use is left alone, so freeing twice is harmless. |
| `int scalar_in_use (scalar)` | whether an id is taken |
| `int scalar_count (void)` | the number of ids in use |
| `int n_global_scalars` | the number of file-scope and static scalars |
| `vector vector_new (void)` | in 2 and 3 dimensions: a vector with a new id for each scalar, `x` first |
| `void vector_free (vector)` | frees its scalars, the last first |
| `tensor tensor_new (void)` | a tensor of new vectors, `x` first |
| `void tensor_free (tensor)` | frees its vectors, the last first |

The pool is a bit array in which a set bit means a free id. The lowest free id
is found with `ctz`: `__builtin_ctz` where it exists, otherwise a de Bruijn
table. The constructor marks every id free before it numbers the globals.

---

## Caches

    #include <cwrap/cache.h>

Including the header turns on `cache (T)` and `foreach` for the rest of the
file. cwrap notices the header through the linemarker gcc writes when the
header starts. Before that line, and in files that never include it, `cache`
and `foreach` are ordinary names.

### The type

`cache (T) v;` is an array of T: `v.a[0]` to `v.a[v.n - 1]`. cwrap turns
`cache (T)` into `struct cache_T`, and defines it once per file, on the line of
the first top-level declaration that uses it:

    struct cache_T {
      __typeof__ (T) *a;
      int n, cap, is_dirty;
      void (*update) (struct cache_T *self);
    };

- **T may be any type:** `cache (int)`, `cache (point)`, `cache (char *)`,
  `cache (cache (int))`.
- **The tag comes from T's words,** with `*` written as `p` and any other
  punctuation as its hex code: `cache (char *)` is `struct cache_char_p`. If
  two different types would get the same tag, that is an error.
- **Different spellings make different types:** `cache (unsigned)` and
  `cache (unsigned int)` are not the same.
- **Use `cache (T)` itself, not a typedef of it.** `typedef cache (int) ints;`
  is an error.

### Starting values

Every field starts at zero: an empty cache with no `update`.

| Where | Who zeroes it |
|---|---|
| File scope, or `static` | C does |
| A local `cache (T) v;` | cwrap adds `= { 0 }`, and a cleanup that frees `v.a` as the scope ends |
| A local with a brace initializer, such as `cache (int) v = { .update = f };` | yours, and cwrap adds the cleanup |
| A member of a struct, or memory from `malloc` | you: `cache_init (v)`, or use `calloc` |

### update and is_dirty

Setting both is up to you. Before a `foreach` over a cache, cwrap runs

    if (v.update && v.is_dirty) { v.is_dirty = 0; v.update (&v); }

once, before the first element is read. An `update` function that sets
`is_dirty` again gets called on the next loop.

    static void rebuild (cache (int) * c)
    {
      cache_clear (*c);
      for (int i = 0; i < 20; i++)
        cache_push (*c, i * i);
    }

    cache (int) squares = { .update = rebuild };
    ...
    squares.is_dirty = 1;
    foreach (s in squares)      // rebuild runs here, once
      total += s;

### API (macros in `include/cwrap/cache.h`)

Each one takes the cache itself, not its address, and evaluates it once.

| Macro | |
|---|---|
| `cache_init (v)` | zero every field |
| `cache_clear (v)` | `n = 0`, keeping the memory |
| `cache_reserve (v, m)` | room for at least `m` elements |
| `cache_push (v, x)` | append `x`. `x` is read after the array grows, so it may be an element of `v`. It may contain commas, as in `cache_push (pts, (point) { 1, 2 })`. |
| `cache_free (v)` | free the array, and set `n` and `cap` to 0 |

The capacity starts at `CACHE_MIN` (8) and doubles. Running out of memory, or
needing more than `INT_MAX` elements, aborts with a message.

### Copies

Two caches must not share one array: a push to one could `realloc` it out from
under the other, and two cleanups would free it twice. So:

| Code | |
|---|---|
| `w = v;` (either side a cache) | error |
| `cache (int) w = v;` | error. Use a brace initializer or `cache_push`. |
| `return v;` | error |
| `f (v)`, a cache passed by value | warning: pass `&v` |

Only variables known by name are checked. `*p = *q;`, or copying a struct that
contains a cache, is not caught.

---

## foreach

    foreach (v in V)
      statement

This becomes, on as many lines as it took:

    for (int __i_v = 0; __i_v < V.n; ++__i_v) { __auto_type v = V.a[__i_v]; statement }

If V is a cache, the index starts at
`V.update && V.is_dirty ? (V.is_dirty = 0, V.update (&V), 0) : 0`, which
updates a dirty cache once before the loop.

- **`foreach` is a keyword only after `#include <cwrap/cache.h>`.**
- **`v` is a copy of each element,** declared in the loop's own scope. It hides
  any outer `v`, including a typedef named `v`, until the loop ends.
- **`in` is not a keyword.** It is recognised by its spelling, so `in` and
  `each` stay usable as names, and even as type names.
- **`foreach (` is always the loop or a call.** Elsewhere, `foreach` can still
  name a variable, a member, a parameter or a function.
- **Declaring `foreach` and also writing a loop is an error.** If a file has a
  `foreach` loop while something named `foreach` is visible, cwrap stops with
  an error and a note pointing at the declaration.
- **V must be a plain name.** `foreach (x in *p)` and `foreach (x in s->list)`
  do not parse yet. See
  [the plan for any object](#future-any-object-as-the-collection-in-foreach).

---

## @loop: block macros

A block macro wraps code around a block that the caller supplies, like a
traversal that visits every node. With `#define` it needs a begin and an end
macro, a backslash on every line, and untyped arguments evaluated as often as
they appear. `@loop` is written as ordinary code instead:

    @loop name (T1 p1, T2 p2) : a, b
    {
      ...
      {...}           /* the caller's block goes here */
      ...
    }

and used like `for`, with any statement as its block:

    name (e1, e2) {
      ... a, b ...
    }

The use becomes, on the lines it took:

    { T1 p1' = e1; T2 p2' = e2; <the body, with {...} replaced by the block> }

For example:

    @loop range (int lo, int hi, int step) : i
    {
      int count = 0;
      for (int i = lo; i < hi; i += step) {
        count++;
        {...}
      }
    }

    range (0, 10, 3)
      total += i;     /* i is the loop's; count is still the caller's */

### The names after `:`

The list after `:` is the set of names the caller's block may see. Every
other name the loop declares, parameters included, is private to the loop.

It is not "the iterators", though the iterators are usually on it. It is the
loop's interface to its block, much as the parameters are its interface to the
caller:

| | The caller gives the loop | The loop gives the block |
|---|---|---|
| declared as | parameters: `(Index root)` | names after `:`: `: IJKL, node, sblocks` |
| at the use | arguments: `morton (TREE_ROOT_INDEX)` | names the block uses: `level`, `node` |
| everything else | hidden | hidden |

In libb-amr's traversal, for instance:

    @loop morton (Index root) : IJKL, node, sblocks
    {
      struct mstack { Index p; int visited; } stack [...];
      int depth = 0, visited = 0;
      Index cur = root;
      ...
          INDEX_UNPACK (cur);     /* int i, j, k, level */
          BAMR_NODE_ENV ();       /* BAMRNode * node; void *** sblocks */
          {...}
      ...
    }

| Name | Listed | What it is |
|---|---|---|
| `i`, `j`, `k`, `level` | yes, through `IJKL` | the index of the node visited: the iterators |
| `node`, `sblocks` | yes | data of the node visited, not iterators |
| `root` | no | the parameter |
| `stack`, `depth`, `cur`, `visited` | no | the traversal's own state |

A block may then use any of these names for itself:

    int depth = 0;                    /* the caller's own depth */
    morton (TREE_ROOT_INDEX) {
      depth++;                        /* counts nodes: the loop's depth is
                                         another variable */
      if (node->flags & BAMR_NODE_LEAF)
        ...
    }

#### Why the names are hidden unless listed

A block macro made with `#define` pastes its body around the block as text,
so any name the body declares is seen by the block. A caller's `depth` above
would silently become the traversal's stack depth, and the block would break
the traversal by changing it. With `@loop`, cwrap renames the loop's own names
at each use. The `range` example above becomes:

    { int __cwrap_range_lo_1 = 0 ; int __cwrap_range_hi_1 = 10 ; int __cwrap_range_step_1 = 3 ;
    {
      int __cwrap_range_count_1 = 0;
      for (int i = __cwrap_range_lo_1; i < __cwrap_range_hi_1; i += __cwrap_range_step_1) {
        __cwrap_range_count_1++;
        total += i;
      }
    }
    }

`i` is listed and keeps its name. `lo`, `hi`, `step` and `count` get names of
their own, `__cwrap_<loop>_<name>_<use>`, which no code can spell by accident.
(The linemarkers cwrap also puts in are left out here.)

#### What can be listed

- **Any name the loop declares that is in scope where `{...}` is:** a loop
  variable, a helper, a typedef. cwrap checks each one: a name the loop does
  not declare, or one whose scope has ended before `{...}`, is an error.
- **A parameter.** `@loop each (const Bag * b) : x, b` lets the block use `b`.
  The argument may then be the caller's own `b`, as in `each (b) ...`: every
  argument is copied to a hidden local first, so `const Bag * b = b;`, which
  would read the new `b`, never happens.
- **Names that a loop used in the body gives it.** In

      @loop walk2 (int d) : level, sum
      {
        int sum = 0;
        walk ((P) { d, 0 }) {     /* walk gives its block level and q */
          sum += level;
          {...}
        }
      }

  the block of `walk2` sees `level`, passed on from `walk`, and `sum`, its own.
  It does not see `q`, which `walk` gives `walk2` but `walk2` does not list:
  inside a definition, a name an inner loop gives is the outer loop's own, and
  hidden like the rest.
- **Nothing.** Without `:`, the block sees none of the loop's names, as for a
  loop that only repeats its block:

      @loop twice ()
      {
        for (int n = 0; n < 2; n++)
          {...}
      }

- **Names that depend on the dimension.** The list is read after the
  preprocessor, so a macro can stand for part of it. `IJKL` gives `morton`
  `i, level`, `i, j, level` or `i, j, k, level` as the dimension has them.

#### The block's view

- **A listed name is the same variable as in the body:** its type and its
  value where `{...}` is, and a change the block makes is seen by the body.
  Changing an iterator in the block changes the loop, as it would in
  hand-written code.
- **A name that is not listed** means whatever it means at the use: a variable
  of the caller, a global, or nothing, in which case gcc reports it as
  undeclared at the line of the block.
- **A listed name that the block does not use** may draw gcc's warning for an
  unused variable, as in hand-written code. Cast it to `(void)` in the body if
  that is expected, as `INDEX_UNPACK` does.

### Rules

- **Arguments are typed.** Each parameter becomes a local of its own type,
  initialised from its argument, so the argument is converted and evaluated
  once, and gcc reports a wrong type at the line of the use. The number of
  arguments must match.
- **The body's names are hidden.** Every name the loop declares, parameters
  included, is renamed for each use, where it is declared and where it is
  used. So the block cannot see them, and they cannot capture the block's own
  names. Renaming follows C's scopes: each name is matched to its declaration,
  not by spelling, and members, designators, labels and tags are left alone.
- **The names after `:` are what the block sees,** and nothing else of the
  loop: see [above](#the-names-after-).
- **`{...}` appears exactly once,** anywhere a statement can go in the body,
  however deeply nested.
- **`break` and `continue` in the block** act on the loop around `{...}` in
  the body, as they would in hand-written code.
- **`return` in the block** returns from the function that uses the loop.
  `return` in the body is an error, as it would leave that function unseen.
- **The body may use other loops,** and pass on names they give it:
  `walk2 (...) : level, sum` can list the `level` that a `walk` in its body
  gives it. A name an inner loop gives but the outer one does not list is the
  outer loop's own, and hidden like the rest.
- **The body is parsed after the preprocessor,** so macros in it are expanded
  and `#if dimension` works inside it.
- **Pragmas in the body go with each use.** `_Pragma ("omp simd")` or
  `#pragma GCC ivdep` before a loop in the body is written before that loop
  in every expansion, and not where the loop is defined.
- **Definitions emit nothing,** so they belong in headers. A definition can
  appear more than once if it is the same each time, as from a header
  included twice. A different one under the same name is an error.
- **Line numbers:** cwrap puts linemarkers around the expansion, so errors and
  the debugger show the body at its lines in the header, the block at its lines
  in your file, and everything after the use keeps its line.
- **`@` needs no header.** C has no other use for it.
- **A loop's name is reserved once defined:** it cannot also name a variable
  or function.

### Errors

cwrap stops with an error, and the file never reaches gcc, for:
- no `{...}`, or more than one
- `{...}` outside a definition
- `return` in a body
- a parameter with no name, or `...`
- a name after `:` that the loop does not declare where `{...}` is
- the wrong number of arguments
- a loop used in its own definition
- a second, different definition with the same name

---

## Stencils: foreach () and s[]

With `-dimension=N` and libb-amr's `<block-loops.h>`, as in Basilisk:

    foreach ()
      lap[] = s[1] + s[-1] + s[0,1] + s[0,-1] - 4*s[];

- **`foreach ()` visits every interior cell of every leaf block.** It is two
  loops, the `@loop`s `foreach_leaf ()` over the leaves, in Morton order, and
  `foreach_cell (tree_block_shape (), 0)` over the cells of each.
- **Stencil loops are declared with `@foreach`,** in `<block-loops.h>`:

      @foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 0);

  the loop over the blocks and the loop over the cells of each, both
  `@loop`s, with their arguments. The block loop must give `sblocks`, the
  cell loop `c` and `shape` (and `level`, `i`, `j`, `k`, `ci`, `cj`, `ck`
  for `x`, `y`, `z` and `Delta`). A loop may take parameters, which its two
  loops' arguments use, and each use binds once:
  `@foreach foreach_in (int l) = leaf_at (l) foreach_cell (...);` then
  `foreach_in (2) s[] = 1;`. Everything below holds for each. The cell loop
  is to visit interior cells only.
- **No stencil loop fills ghosts.** What a ghost is depends on what is
  beyond the block: a leaf of the same level, a coarser one, the domain's
  boundary, another rank later. `ghosts_update ()` (block-ghosts.h) does it,
  by the scalar's operations, and `foreach ()` calls it as needed (below).
  A derived scalar is computed on the interior with `foreach ()`, and its
  ghosts filled like any other's. Loops over the ghosts, by kind of
  boundary, are left for later, with MPI.
- **`s[]` is the cell of scalar `s`; `s[i]`, `s[i,j]`, `s[i,j,k]` the cell at
  that offset.** Offsets left out are 0, so `s[1]` is `s[1,0,0]` in any
  dimension. More offsets than the dimension is an error.
- **No access leaves its block, checked at compile time.** Each offset must
  be an integer constant expression (`1`, `-2`, `1 - 2`, an enum constant)
  within `BAMR_GHOSTS`, the ghost layers, which `<block-loops.h>` exports as
  an enum constant, since `GHOSTS` is gone after `gcc -E`. In
  `foreach_neighbor (r)` an offset plus `r` must stay within them too. gcc
  checks it, at the access's line and column, with a static assertion in the
  index (`(int) (0 * sizeof (struct { _Static_assert (...); ... }))`), which
  costs nothing:

      foreach () s[] = t[3];     /* error: static assertion failed: "t[3]: an
                                    offset is not a constant within
                                    BAMR_GHOSTS, the ghost layers" */
      foreach () s[] = t[n];     /* error: expression in static assertion is
                                    not constant */

  A data-dependent stencil is written with constant offsets:
  `u[] > 0 ? s[-1] : s[1]`, which also vectorises (a select, not a gather).
- **Only a scalar is expanded.** `s` must be an expression of type `scalar`
  that cwrap can follow: a scalar name or parameter, a vector's or tensor's
  member (`u.x[]`, `t.x.y[]`), an element of an array of scalars or of a
  pointer to them (`list[k][]`, `(*p)[]`). `x[i]` of anything else is the
  ordinary subscript; `x[]` of anything else is an error.
- **Each scalar the body reads or writes gets its block once per leaf,** before
  the cells, and no other scalar does. A scalar the body only reads gets a
  `const` block:

      foreach_leaf () {
        Real *const       __cwrap_block_1_0 = block_data (sblocks, lap);
        const Real *const __cwrap_block_1_1 = block_data (sblocks, s);
        foreach_cell (tree_block_shape (), 0)
          __cwrap_block_1_0[c] = __cwrap_block_1_1[c + (1)] + ...
                               + __cwrap_block_1_1[c + (0) + (1) * shape.sy] ...;
      }

  `s[] = `, `s[] += `, `s[]++` and `&s[]` count as writes. Two spellings of
  one scalar (`s` and `u.x` with the same id) get a block each.
- **A scalar the body may change is found at every cell instead,** as
  `block_data (sblocks, list[k])[...]`: one that uses a name declared in the
  body, such as a loop index, or a name the body assigns.
- **`x`, `y`, `z` and `Delta` are the cell's,** as in Basilisk: the
  coordinates of its centre and its size, on the domain that
  `tree_shift_resize ()` sets. `y` and `z` beyond the dimension are 0.

      foreach ()
        h[] = exp (-(x*x + y*y)/(Delta*Delta));

  - **They hide any outer `x`, `y`, `z` or `Delta`** inside the body, while a
    name the body declares itself (`double x = ...;`) is its own. `u.x` is a
    member, not the coordinate.
  - **Only those the body uses are computed:** `Delta` and the centre of the
    leaf's first cell once per leaf, then one multiply-add per cell and
    coordinate:

        const double __cwrap_Delta_1 = tree_get_size () / ((double) tree_get_block_m ()
                                         * (double) (1L << level));
        const double __cwrap_x_1 = tree_get_origin (0) + ((double) i * tree_get_block_m ()
                                     + 0.5) * __cwrap_Delta_1;
        ...
          const double x = __cwrap_x_1 + (ci - shape.ox) * __cwrap_Delta_1;

  - **They cannot name a scalar** in `foreach ()`: `Delta[]` is an error,
    and so is any use of an outer scalar, vector or tensor called `x`, `y`,
    `z` or `Delta` (`f (x)`, `u = y`), which the coordinate would hide.
    Rename it. An outer variable of another type is hidden quietly, as in
    Basilisk.
  - **They need `foreach_leaf` to give `level`, `i`, `j`, `k`, and
    `foreach_cell` to give `ci`, `cj`, `ck`,** as in `<block-loops.h>`.
- **The body sees none of the loops' names.** `c`, `shape`, `sblocks`, `node`,
  `level`, `i`, `j`, `k`, ... are renamed in the expansion, so the body's
  own, or those around the loop, keep their meaning.
- **The loop is serial over the leaves,** and `foreach_cell` vectorises each
  row as `BLOCK_BACKEND` says (see `block-loops.h`): a sum into a variable
  is fine; a cell must not be read along i in the row that writes it.
- **No `break`:** it leaves only the row of cells.
- **Ghosts are kept up to date lazily,** as in Basilisk. After the loop,
  each scalar it writes is marked stale (`ghosts_mark_dirty ()`, a bit).
  Before a loop that reads a scalar at an offset (`p[1]`, `p[0,k]`: any
  offset not written as `0`), those of them marked are updated, all in one
  `ghosts_update_dirty ()`; a scalar read only at its cell (`p[]`) needs
  none. Any update clears the bits, `ghosts_update ()` by hand or
  `adapt_wavelet ()` too.

      foreach () p[] = x*x;            /* p marked */
      foreach () q[] = p[1] - p[-1];   /* p updated first; q marked */

  - **A scalar that changes in the body,** `list[k]` with `k` the body's,
    stands for every element of `list` if it is an array of a known size,
    and for every scalar otherwise (`ghosts_mark_all_dirty ()`,
    `ghosts_update_dirty (-1, 0)`).
  - **Cells written by hand,** through `block_data ()`, are not seen: call
    `ghosts_update ()` as before, or `ghosts_mark_dirty ()`.
  - **Writing a cell at an offset** (`p[1] = ...`) **is an error:** a cell
    writes only its own.
  - **Reading a neighbour of a scalar the same loop writes is a warning:**
    `foreach () { p[] = x; q[] = p[1,0] - p[-1,0]; }` reads the ghosts of
    `p` stale, and its neighbours old or new by the order of the cells.
    Write `p` in a `foreach ()` of its own first.
  - **Loops whose ghosts the caller keeps:** a stencil loop whose block
    loop gives `ghosts_by_caller`, as libb-amr's `foreach_level (l)` and
    `foreach_level_or_leaf (l)` (one level of the tree, for multigrid),
    keeps every rule above, but cwrap writes no update or mark around it.
    Its caller updates the ghosts of that level (`ghosts_update_level (s,
    l)`) and marks what it wrote.
  - **Loops over one colour of a checkerboard:** a stencil loop whose cell
    loop gives `checkerboard`, as libb-amr's
    `foreach_level_or_leaf_colour (l, colour)` (red-black Gauss-Seidel),
    may read a scalar it writes at an odd sum of offsets: `p[1]`,
    `p[0,-1]`, `p[2,1]`. Those cells have the other colour, which the loop
    does not change, so there is no warning. An even sum (`p[1,1]`,
    `p[2]`) is the cell's own colour, possibly changed already: a compile
    error, from a static assertion gcc checks as it checks the reach, so
    offsets written as constant expressions are checked too.
    `foreach_neighbor ()` reading such a scalar is an error: it reaches
    both colours. With `BAMR_BOTH` the loop visits every cell in a fixed
    order, and the same reads are Gauss-Seidel in place.
- **`foreach_dimension ()` writes a statement once for every axis,** as in
  Basilisk. The statement is written for x. In the next copy x becomes y,
  y becomes z, z becomes x (in 2D x and y swap), and so on. This applies
  to the names `x`, `y`, `z` (coordinates, members `.x`, `child.x`), to
  names ending in `_x`, `_y`, `_z`, and to the offsets of `s[...]` and
  `coarse (s, ...)`, which are permuted the same way:

      foreach ()
        foreach_dimension ()
          g.x[] = (p[1] - p[-1]) / (2 * Delta);  /* g.y[] from p[0,1], p[0,-1] */

  The copies follow the statement on its last line, so no line moves. In
  a stencil loop their accesses are the loop's, with every rule. It works
  outside stencil loops too (`u_x`, `v.x`), and in 1D it is the statement
  as written. libb-amr's `prolong_linear` and `prolong_minmod` are written
  with it.
- **Loops by axis:** where `N_x`, `N_y` (, `N_z`) are `@foreach` loops and
  `N` is not, `N () statement` runs `N_x ()` with the statement, then
  `N_y ()` with it rotated as `foreach_dimension ()` rotates, and so on.
  Each is a loop of its own, with every rule. libb-amr's `foreach_face ()`
  is one. A cell loop that gives `face_axis` has `x`, `y`, `z` at the
  face: half a cell down that axis.
- **A comma in `s[i,j]` splits a macro's argument:** the preprocessor
  does not group brackets. Read the cell into a local first, or put it in
  parentheses: `fabs ((v[1,1]) - w)`.
- **Faces and vertices read at their cell are kept up to date too.** A
  face vector's own values change when the fluxes are made to agree at
  level jumps. So before a loop, the scalars it reads at their cell are
  passed to `ghosts_update_dirty_located ()`, where that is declared. It
  updates those that are stale and not at the centre (`scalar_location
  ()`), and is cheap for the others.
- **`foreach_neighbor (r)` moves the point,** as in Basilisk: its statement
  runs for each cell within `r` of the point along every axis, the point's
  own included, and `s[]`, `s[i,j]`, `x`, `y`, `z` and `Delta` are then
  those of that cell.

      foreach () {
        double sum = 0;
        int n = 0;
        foreach_neighbor (1) {
          sum += s[];
          n++;
        }
        mean[] = sum / n;
      }

  - **It reads the ghosts,** which hold the neighbours: no tree is walked.
    Its reads are reads at an offset, so the lazy ghost update applies.
  - **It only reads:** writing a cell in it (`s[] = ...`, `s[]++`, `&s[]`)
    is an error, as it would write a neighbour. Its accesses are also
    written as values, `((Real) block[...])`, not lvalues, so gcc refuses a
    write there too.
  - **Only in a stencil loop, and not in another `foreach_neighbor ()`:**
    either is an error.
  - **`r` is an integer constant from 0 to `BAMR_GHOSTS`,** checked at
    compile time like the offsets, and each offset inside it plus `r` too.
  - The point movers to come, as `foreach_child ()` in loops over parents,
    follow the same rules.
- **Stencil loops do not nest, which is checked at compile time:**
  - one inside another, of the same `@foreach` or not, is an error;
  - so is a call, inside one, of a function of the same file that runs one,
    itself or through other functions of the file. A function of another
    file, or one called through a pointer, cannot be seen.
  - `s[]` outside a stencil loop is an error.
- **`foreach ()` is the loop only while no `foreach` is declared** (a
  function, a variable). Otherwise `foreach ()` is the call it was.
- **Errors:** without the `@foreach` of a loop, as without
  `<block-loops.h>`, its use is an error that says to include it; a loop
  that writes, or reads at an offset,
  without `ghosts_mark_dirty ()`, `ghosts_mark_all_dirty ()` and
  `ghosts_update_dirty ()` declared, one that says to include
  `<block-ghosts.h>`.

---

## Ghost operators: @prolongation, @restriction, @boundary

The operators that fill ghosts and refined blocks are chosen per scalar by
an integer id, and called through a `switch` whose cases are inlined: no
function pointer, so they inline, vectorise, and can be offloaded.

    @prolongation prolong_minmod (const Prolongation * pr) { ... }
          /* static inline void __cwrap_prolongation_prolong_minmod (...) { ... }
             const int prolong_minmod = 2; */
    @prolongation prolong_minmod;
          /* extern const int prolong_minmod; */
    @prolongation (id) (&pr);
          /* switch (id) { case 1: __cwrap_prolongation_prolong_linear ((&pr)); break;
                           case 2: __cwrap_prolongation_prolong_minmod ((&pr)); break;
                           ... default: __builtin_trap (); } */

and the same with `@restriction` and `@boundary`.

- **A definition** is at file scope, like a function: its parameters, its
  body. Its name is the operator's id, an `int`: ids count from 1 per kind,
  in the order of the definitions in the file; 0 is left for "the default".
  A name defined twice is an error.
- **A declaration** lets any other file use the id by name.
- **A dispatch** is a statement: the id, then the arguments, evaluated
  once. It reaches every operator of its kind defined before it in the
  file, and none defined after; with none, it is an error. An id that is
  none of them stops the run (`__builtin_trap`).
- **So one file of a program dispatches,** and must see every definition.
  In libb-amr that is `block-ghosts.c`, which includes the operators given
  (`block-ops.h`, the defaults first, so that their ids are the same in
  every program) and then the core (`block-ghosts-impl.h`). A program with
  operators of its own builds its own such file instead of
  `block-ghosts.c`; see libb-amr's `block-ghosts.h` and `tests/test_ghosts.c`.
- **A scalar is given operators anywhere, at any time,** by their ids:
  `scalar_set_ops (h, (ScalarOps) { .prolongation = prolong_minmod })`. That
  marks its ghosts stale, to be made again by the new ones.

### Coupled operators

An operator that also reads other scalars names them in a second parameter
list: they are bound per scalar, and passed to it at each call.

    @prolongation prolong_lake (const Prolongation * pr) (scalar zb)
    {
      scalar s = pr->s;
      foreach_prolongation (pr)
        s[] = coarse (s) + coarse (zb) - zb[];   /* the surface, less the bed */
    }

    scalar zb, h;                 /* global: see below */
    ...
    ghost_ops (h, prolongation (prolong_lake, zb), boundary.all (bc_odd),
                  boundary.left (bc_inflow));

- **`ghost_ops (s, kind (op, scalars...), ...)`** binds operators to `s`,
  anywhere, at any time, each with the scalars it depends on: kinds
  `prolongation`, `restriction`, and a boundary condition per face of the
  domain, `boundary.left`, `.right`, `.bottom`, `.top`, `.back`, `.front`,
  or `boundary.all` for every face. It becomes a call of libb-amr's
  `ghost_ops_set ()` per item, in order, so a face after `boundary.all`
  replaces it there.
- **A boundary names its face:** a plain `boundary (op)` is an error, as is
  a face the dimension lacks (`.back` in 2D) or no face at all
  (`.middle`); `boundary.all` after a face, in one call, is a warning: it
  binds that face again.
- **The scalars depended on must be global** (at file scope, `static` or
  `extern`, or a member of a global vector or tensor), checked by cwrap:
  a global scalar is never freed (see [Scalars](#a-global-scalar-is-never-freed)),
  so a binding cannot outlive what it reads. A local one is an error.
- **The dispatch passes them:** `@prolongation (id, uses) (args)` calls one
  that depends on scalars with `uses [0]`, `uses [1]`, ... after `args`;
  without `uses`, it is an error if one of its kind depends on any.
  `@prolongation uses (id)` is how many scalars that operator takes.
- **What libb-amr then checks, as `ghost_ops ()` runs:** that the operator
  takes that many scalars, that they are global, and that no scalar's
  operators would depend on it itself through others. Either stops the run
  with a message, binding nothing. `scalar_set_ops ()` refuses an operator
  that depends on scalars.
- **And keeps:** wherever several scalars' ghosts are made (a ghost update,
  a refinement, a coarsening), those depended on are made first, whatever
  the ids, and a stale one is updated with the scalar that needs it;
  marking one stale (`ghosts_mark_dirty`) marks those that depend on it.

### Operator bodies: stencils over two levels

A prolongation or a restriction is written with the stencil syntax, in a
loop over two levels that libb-amr declares with `@foreach`:

    @prolongation prolong_linear (const Prolongation * pr)
    {
      scalar s = pr->s;
      foreach_prolongation (pr) {           /* the fine cells pr asks for */
        RealSum v = coarse (s);
        v += child.x * (coarse (s, 1) - coarse (s, -1)) / 8.;
        v += child.y * (coarse (s, 0, 1) - coarse (s, 0, -1)) / 8.;
        s[] = (Real) v;
      }
    }

    @restriction restriction_mean (Index p, scalar s)
    {
      foreach_restriction (p) {             /* the parent's cells, child by child */
        RealSum sum = 0;
        foreach_child ()
          sum += s[];
        s[] = (Real) (sum / TREE_NCHILDREN);
      }
    }

- **What a loop over two levels is, cwrap tells by the names its two loops
  give**, as `sblocks`, `c` and `shape`:
  - `coarse_sblocks` (block loop), `cc` and `sgn` (cell loop): the point
    has a coarse cell. `coarse (s)` is it, `coarse (s, i, j, k)` its
    neighbours at that level; `child.x`, `.y`, `.z` the half of it the
    point is in, -1 or +1 (0 beyond the dimension).
  - `child_sblocks` (block loop) and `fc` (cell loop): the point has
    children. `foreach_child ()` moves the point to each, in order, i
    fastest, as `foreach_neighbor ()` does.
- **`s[]` is the point's cell, and the only cell of its level:** an offset
  there (`s[1]`) is an error, those cells being the ones filled, and so is
  `foreach_neighbor ()`.
- **`coarse ()` and `foreach_child ()` only read,** and are written as
  values; a write is an error. `coarse ()` offsets are constants within
  `BAMR_COARSE_REACH` (1: only the coarse block's first ghost layer is valid
  at a level jump), checked at compile time as the others.
- **`coarse ()` outside such a loop, as in `foreach ()`, and `foreach_child
  ()` in a loop whose points have no children, are errors.**
- **No ghost is tracked:** these loops run inside the ghost update itself.
- **Each scalar's blocks are found once per block,** fine and coarse, or
  parent and child, as in `foreach ()`, and the restriction runs child by
  child, so that the parent block and the child's are fixed in the inner
  loop.

libb-amr's given operators are written so (`block-ops.h`): the same results
as before, bit for bit, and as fast or faster, as nothing is called through
a pointer. The boundary conditions stay plain C, with `boundary_cells ()`.

---

## How the translator works

1. **Lexer** (`src/lexer.l`, reentrant flex)
   - Each token keeps the whitespace, linemarkers and `#pragma`s before it as
     trivia, along with its file and line.
   - It asks the scope table whether an identifier is a typedef name or a
     loop, and records on its token which declaration it named, which is
     what renaming in a `@loop` follows.
   - A linemarker naming `.../cwrap/cache.h` turns `cache` and `foreach` into
     keywords.
2. **Parser** (`src/parser.y`, pure bison, `%expect 0`)
   - It covers C11 with the GNU extensions found in system headers.
   - Empty marker rules open and close scopes at the right token.
   - On syntax errors, error recovery lets it go on and find more, up to 20.
3. **Scope table** (`src/scope.c`)
   - It holds a hash table of ordinary names for each scope, which is what lets
     `T * x;` be parsed correctly.
   - The parser checks names against it for the scalar and cache rules.
4. **Translation** (`src/translate.c`)
   - Some rewrites happen while the parser reduces a rule:
     - `foreach`, before the loop's scope closes
     - `foreach ()`: each `s[...]` of a scalar as it reduces, then the two
       `@loop` expansions once the body is read
     - each use of a `@loop`, before the scope with the names it gives the
       block closes
     - local scalars and caches, as each declaration in a block is reduced
     - `cache (T)` itself
   - The rest happens in a pass over the finished tree: the global scalar slots,
     and the `struct cache_T` definitions.
   - Code that cwrap adds is written as C text and parsed with the same parser
     (`parse_snippet`), so it becomes real subtrees.
5. **Output:** the tree is written back token by token, with each token's
   trivia. A file that cwrap did not change comes out byte for byte the same as
   it went in.

---

## Source layout

| File | |
|---|---|
| `src/driver.c` | the command line, running gcc for each stage, temporary files, linking |
| `src/lexer.l`, `src/lexer.h` | tokens and their trivia |
| `src/parser.y` | the grammar, error recovery, and calls into the translation |
| `src/ast.c`, `src/ast.h` | the tree, its editing functions (`ast_replace`, `ast_insert_before`), output and dump |
| `src/scope.c`, `src/scope.h` | names by scope |
| `src/translate.c`, `src/translate.h` | the scalar, cache, foreach and `@loop` rewrites and checks |
| `src/diag.c`, `src/diag.h` | gcc-style colored diagnostics |
| `src/allocator.c`, `src/allocator.h` | the arena that holds each file's tree |
| `src/scalar.c`, `src/scalar.h` | the scalar runtime |
| `include/cwrap/cache.h` | the cache API |
| `tests/run.sh`, `tests/cases/` | end-to-end tests and their files |
| `tests/scope_test.c` | unit tests of the scope table |
| `tests/grammar.c` | a file using the whole grammar, which must come out byte for byte the same |

---

## Tests

    make test

`tests/run.sh` copies `tests/cases/` into a scratch directory and runs cwrap
on it. The tests cover:
- compiling and linking as gcc does
- byte-for-byte output for files cwrap does not change
- line numbers in errors
- error recovery
- every `foreach` case
- numbering scalars across files
- local scalars and their cleanup
- vectors and tensors in 2 and 3 dimensions, and scalars off without
  `-dimension`
- caches: running them, their translation, the errors and warnings, and
  having them off without the header
- colored diagnostics
- `@loop`: running uses of it, every error, line numbers around an expansion,
  and its translation
- `foreach ()` on a small grid of its own (`stencil.h`): offsets, a vector,
  the body's own names, `x`, `y` and `Delta`, the ghosts updated only when
  stale, the dependency warning, `foreach_neighbor ()` and its checks,
  `@foreach` with parameters, nesting through calls, loops over two levels
  (`coarse ()`, `child.x`, `foreach_child ()`) and their errors, ghost
  operators (`ops_run.c`, `ops_bad.c`) and coupled ones (`ops_coupled.c`,
  `ops_coupled_bad.c`), a global scalar never freed, a block per scalar
  used, `const` when only read, line numbers, every error, and a function
  named `foreach`. libb-amr's
  `tests/test_stencil.c` checks it on the real tree against `cell_at`.
- nested `@loop`: 24 self-checking cases (`loop_nest.c`), each against the
  same computation in plain loops: a loop in its own block, five deep,
  arguments that read an outer loop's names, each argument evaluated once,
  definitions three levels deep, names passed on from two levels down, break
  and continue at every level, parameters the block sees, statics, caches and
  scalars; with their errors and line numbers

It prints `ok` or `FAIL` for each check, then `all passed`.

---

## Limitations

- **`-E` with `-dimension`** goes straight to gcc without `scalar.h`, so
  `scalar`, `vector` and `tensor` are unknown there.
- **Linux only:** global scalars rely on GNU ld's `__start_`/`__stop_` symbols
  for an ELF section.
- **Not portable source:** the output of `-translate` is preprocessed, so it
  contains this system's headers and fully expanded macros.
- **GCC extensions in the generated code:** `__attribute__ ((cleanup))`,
  `constructor`, `section`, `__auto_type`, `__typeof__`, and statement
  expressions in `cache.h`.
- **`longjmp` skips the frees** of local scalars and caches, as it skips any
  `cleanup`.
- **The collection in a `foreach` must be a name.** See
  [the plan](#future-any-object-as-the-collection-in-foreach).
- **Cache copies are checked only for variables known by name.**
- **The element type of a cache must be visible at file scope** where its
  struct is defined. A type defined inside a function does not work.
- **The tag table does not hold struct members,** so a `foreach` over a struct
  that is not a cache uses `__auto_type` rather than the real element type.
- **Scalars in a struct** (`p.h[]`, with `struct { scalar h; } p`) are not
  followed: the tag table has no members.
- **In a `@loop` body:**
  - **Labels are not renamed,** so a loop with a label can be used once per
    function. GNU's `__label__` declarations make them local.
  - **A global the body uses can be hidden by a local of the caller** with
    the same name, as the body is pasted where the loop is used. cwrap does
    not check for it.
  - **A `static` local is one per use,** as the body is copied each time.
  - **Columns can shift** after a renamed name on the same line, as
    `__cwrap_range_count_1` is longer than `count`. Lines are always right.
  - **A parameter whose type defines a struct,** as in
    `@loop f (struct { int a; } s)`, defines it at each use, and twice if the
    block sees the parameter, as two types that do not match. Name the
    struct, or use a typedef.

---

## Future: a single build, and portable C99 output

The goal is a mode, perhaps `-export`, that turns a whole program into **one
self-contained C99 file**. That file could then be built on machines without
cwrap or GCC extensions, such as supercomputers. It is planned, not started.

### Why one build

Several of today's mechanisms exist only because files are compiled
separately:
- the linker section and constructor that number global scalars
- `scalar.o` as a separate runtime

If cwrap sees every file at once, both can be replaced by plain C.

### Plan

1. **Number scalars as constants.**
   - With every file in view, cwrap writes `scalar h = 0, u = 1, ...;` itself.
   - The pool's first n bits are a static initializer, or a call placed first
     in `main`.
   - No constructor, `section`, `weak` or `used` is needed.
2. **Fake standard headers before preprocessing.**
   - `#include <stdio.h>` is pointed at a stand-in that declares only the type
     names the parser needs (`typedef int FILE;` and so on), plus a
     `#pragma cwrap include <stdio.h>` marker.
   - After translation, cwrap replaces each stand-in's region with the real
     `#include`. Standard macros (`NULL`, `EOF`, `assert`, `INT_MAX`) are left
     unexpanded, for the target machine's headers to expand.
   - Macros that take a type, such as `va_arg` and `offsetof`, are mapped to
     builtins the grammar knows, and written back afterwards.
   - `#if` tests on header macros are reported, since the stand-ins decide
     them.
   - cwrap's own `scalar.h` and `cache.h` stay inlined: they are the runtime.
3. **Join the files.**
   - All `.c` files go through one preprocessing pass, so include guards remove
     repeated headers. That also avoids C99's error on a repeated typedef.
   - File-scope `static` names are renamed per file.
   - The runtime from `scalar.c` is emitted once, as `static` code.
4. **Find types through the scope stack,** to drop `__auto_type` and
   `__typeof__`.
   - `foreach` over a cache uses T, which cwrap already knows.
   - `foreach` over any other struct follows V's declaration to its struct, then
     to the member `a`. That needs a tag table with member lists.
   - The `cache.h` macros become typed functions that cwrap generates, such as
     `__cwrap_cache_int_push (&v, x)`. That also removes the statement
     expressions. It needs a small type checker for names, `*`, `&`, `.`, `->`
     and `[]`.
   - A type that cannot be found is an error, not a guess.
5. **Free locals at every exit, instead of using `cleanup`.**
   - The end of a block: when its scope is popped.
   - `return e;`: store `e` in a temporary first, typed from the function's
     declarator, then free, then return.
   - `break` and `continue`: free each scope up to the nearest loop, or
     `switch` for `break`. Scopes will record which kind they are.
   - `goto`: free the scopes the jump leaves. Forward labels are fixed up when
     the function ends.
   - Errors:
     - jumping into a scope past a declaration that needs freeing, by `goto` or
       a `switch` `case`
     - `longjmp`
6. **Report what cannot be lowered,** such as GNU C in the user's own code
   (`__attribute__`, case ranges, `?:` with the middle operand left out),
   instead of passing it through.

The order of work: constant ids, then the fake headers, then types, then frees
at every exit. Frees at every exit is the largest piece, mostly because of
`goto` and jumps past declarations.

---

## Future: any object as the collection in foreach

Today V in `foreach (v in V)` must be a plain name. The plan is to accept any
expression that names an object, such as:
- `foreach (x in mesh->cells)`
- `foreach (x in grid.levels[l])`
- `foreach (x in *pc)`
- `foreach (x in caches[k])`

It is planned, not started. It waits for the same type work as step 4 of the
C99 plan above: types and declarations followed through the whole scope stack.

### Grammar

V becomes a `unary_expression`, the same nonterminal as the left side of `=`:

    FOREACH '(' IDENTIFIER IN unary_expression ')' statement

That covers names, `a[i]`, `a.x`, `a->b` and `*a`, and rejects `a + b` or
`x = y` while parsing. `IN` and `)` close it off, so it should keep `%expect 0`.

### Lowering: evaluate V once

The current lowering repeats V for `update`, `is_dirty`, `n` on every pass, and
`a[i]`. That is wrong for `caches[k++]` or `next_block ()->c`. V will be bound
once to a pointer:

    for (__typeof__ (V) *__c_v = &(V); __c_v; __c_v = 0)
      for (int __i_v = (__c_v->update && __c_v->is_dirty
                        ? (__c_v->is_dirty = 0, __c_v->update (__c_v), 0) : 0);
           __i_v < __c_v->n; ++__i_v) {
        __auto_type v = __c_v->a[__i_v];
        statement
      }

- `break` leaves the inner loop, and the outer one stops because the pointer
  is now null. `continue` is unchanged.
- `__c_v->a` is read on every pass, so pushing to V inside the loop stays safe.
- V must be an lvalue. Returning a cache is already an error, so that costs
  nothing.
- It fixes `foreach (c in c.kids)`. Today that expands to
  `__auto_type c = c.kids.a[...]`, where `c` in the initializer is already the
  new `c`. Binding the pointer before `v` is declared avoids it.

### Knowing that V is a cache

Today `is_cache` comes from looking up V's name in the scope stack. For
`p->mesh.cells[2]`, cwrap needs the type of every step. In stages:

1. **Assume a cache** when V is not a plain name. Every `struct cache_T` has
   the same layout, so the lowering works without knowing T. If V is not a
   cache, gcc reports a missing member `n`: unclear, but caught at compile
   time. A plain name keeps today's lookup.
2. **Record cache members.** cwrap already sees every `cache (T)`, including
   inside struct definitions, so it can note "struct mesh has a cache member
   `cells`". That checks `m.cells` and `p->cells` whenever the type of `m` or
   `p` is known from a declaration in scope. Any step it cannot resolve falls
   back to stage 1.
3. **Follow every type through the scope stack:** a tag table with members,
   pointer types of locals and parameters, and a small type checker for names,
   `*`, `&`, `.`, `->` and `[]`. This is the same piece the C99 export needs,
   so the two should be built together.

Still open: whether a V that is not a name means "cache only", or also covers
other structs with `n` and `a`, and arrays reached through a path, such as
`foreach (x in s.arr)`. Arrays need their length from the type, which again
means stage 3.
