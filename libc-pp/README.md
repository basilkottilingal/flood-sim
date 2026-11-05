# c wrapper

*Documentation as of 2026-11-05.*

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

---

## Contents

1. [Building](#building)
2. [Using cwrap](#using-cwrap)
3. [Diagnostics](#diagnostics)
4. [Scalars](#scalars)
5. [Caches](#caches)
6. [foreach](#foreach)
7. [How the translator works](#how-the-translator-works)
8. [Source layout](#source-layout)
9. [Tests](#tests)
10. [Limitations](#limitations)
11. [Future: a single build, and portable C99 output](#future-a-single-build-and-portable-c99-output)

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

### File-scope scalars

    scalar h, u, v;      // ids 0, 1, 2 ... in the order the linker sees them

- **Each one gets a slot,** on its own line, so no line numbers move:

      scalar h; static scalar *const __cwrap_scalar_slot_0
        __attribute__ ((section ("scalars"), used)) = &h;

  A vector or tensor gets one slot for each of its scalars. In 2 dimensions,
  `tensor t;` gets four, all on its line:

      static scalar *const __cwrap_scalar_slot_0 __attribute__ ((section ("scalars"), used)) = &t.x.x;
      static scalar *const __cwrap_scalar_slot_1 __attribute__ ((section ("scalars"), used)) = &t.x.y;
      static scalar *const __cwrap_scalar_slot_2 __attribute__ ((section ("scalars"), used)) = &t.y.x;
      static scalar *const __cwrap_scalar_slot_3 __attribute__ ((section ("scalars"), used)) = &t.y.y;

- **The linker joins the slots of every file** between `__start_scalars` and
  `__stop_scalars`.
- **A constructor in `scalar.c`, at priority 101, numbers them** before `main`
  and before any other constructor. So scalars across separately compiled
  files get ids 0 to n-1 with no gaps.
- **`n_global_scalars`** holds n.
- **Not numbered:**
  - `extern` declarations
  - `typedef`s
  - pointers, arrays and functions
  - names with an initializer of their own, which also get a warning

### Scalars in functions

    void f (void)
    {
      scalar x;                 // x = scalar_new (), freed as the scope ends
      static scalar s;          // numbered at startup, like a global
      scalar k = scalar_new (); // yours: freed only when you call scalar_free
    }

- **`scalar x;` becomes**

      scalar x __attribute__ ((cleanup (__cwrap_scalar_free))) = scalar_new ();

  so `scalar_free (x)` runs however the scope ends: the end of the block,
  `return`, `break`, `continue` or `goto`. Scalars are freed in reverse order.
- **A vector or tensor gets an id for each scalar,** and a cleanup that frees
  them, last first:

      tensor w __attribute__ ((cleanup (__cwrap_tensor_free)))
        = { { scalar_new (), scalar_new () }, { scalar_new (), scalar_new () } };

- **`static scalar s;`** becomes `static scalar s;` followed by one slot for each
  of its scalars, as declarations of their own on the same line. They are
  numbered with the globals.
- **Left alone:** pointers, arrays, names with an initializer, and `extern`,
  `register` and thread-local scalars.

### Changing a scalar

Assigning to a scalar is a **warning**, because cwrap gives scalars their
values. This covers `=`, the compound assignments, and `++` and `--`. It also
covers a whole vector or tensor, or part of one, as in `v = w`, `t.y = v` and
`t.x.y = 1`. The warning names the operator and what it changes, and a note
shows the declaration.

### Runtime API (`src/scalar.h`)

| Function | |
|---|---|
| `scalar scalar_new (void)` | the lowest free id. Aborts with a message if all `NSCALAR_MAX` are in use. |
| `void scalar_free (scalar)` | frees an id. An id not in use is left alone, so freeing twice is harmless. |
| `int scalar_in_use (scalar)` | whether an id is taken |
| `int scalar_count (void)` | the number of ids in use |
| `int n_global_scalars` | the number of file-scope and static scalars |

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
| `cache_push (v, x)` | append `x`. `x` is read after the array grows, so it may be an element of `v`. |
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
  do not parse yet.

---

## How the translator works

1. **Lexer** (`src/lexer.l`, reentrant flex)
   - Each token keeps the whitespace, linemarkers and `#pragma`s before it as
     trivia, along with its file and line.
   - It asks the scope table whether an identifier is a typedef name.
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
| `src/translate.c`, `src/translate.h` | the scalar, cache and foreach rewrites and checks |
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
- **The collection in a `foreach` must be a name.**
- **Cache copies are checked only for variables known by name.**
- **The element type of a cache must be visible at file scope** where its
  struct is defined. A type defined inside a function does not work.
- **The tag table does not hold struct members,** so a `foreach` over a struct
  that is not a cache uses `__auto_type` rather than the real element type.

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
