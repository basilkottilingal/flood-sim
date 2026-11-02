#include <stdint.h>
#include <string.h>

#include "allocator.h"
#include "scope.h"

/* Bucket counts are powers of two. An empty scope has no table, since most
 * blocks declare nothing. The limit keeps a table within what allocate () will
 * grant: past it the chains just get longer. */
enum { MIN_BUCKETS = 8, MAX_BUCKETS = 1 << 15, NSIZES = 13 };

struct scope {
  struct scope *up;
  struct var **bucket;
  unsigned nbuckets, n;
};

struct scope_table {
  struct scope *top;
  struct scope *fn_params;
  int param_depth;
  int typedef_seen;
  enum scope_kind declaring;
  /* Tables no scope uses, by size: free_tables[i] has MIN_BUCKETS << i buckets.
   * The arena cannot take them back, so they wait here to be used again. The
   * first slot of a free table links to the next one. */
  struct var **free_tables[NSIZES];
};

/* Types gcc knows without any declaration, which reach us unexplained from
 * system headers. */
static const char *const builtin_types[] = {
  "__builtin_va_list", "__int128_t", "__uint128_t",
  "__float80", "__float128", "__ibm128",
  "_Float16", "_Float32", "_Float64", "_Float128",
  "_Float32x", "_Float64x", "_Float128x",
  "_Decimal32", "_Decimal64", "_Decimal128",
};

static unsigned hash_name (const char *name)
{
  uint32_t h = 2166136261u;       /* FNV-1a */

  for (const unsigned char *p = (const unsigned char *)name; *p; p++)
    h = (h ^ *p) * 16777619u;
  return h;
}

static struct var *find_in (const struct scope *sc, const char *name,
                            unsigned hash)
{
  if (!sc->nbuckets)
    return NULL;
  for (struct var *v = sc->bucket[hash & (sc->nbuckets - 1)]; v; v = v->chain)
    if (v->hash == hash && strcmp (v->name, name) == 0)
      return v;
  return NULL;
}

static unsigned size_index (unsigned nbuckets)
{
  unsigned i = 0;

  while ((unsigned)MIN_BUCKETS << i < nbuckets)
    i++;
  return i;
}

/* An empty table of n buckets: a free one if there is any. */
static struct var **take_table (struct scope_table *st, unsigned n)
{
  unsigned i = size_index (n);
  struct var **table = st->free_tables[i];

  if (!table)
    return allocate (n * sizeof (*table));
  st->free_tables[i] = (struct var **)table[0];
  memset (table, 0, n * sizeof (*table));
  return table;
}

static void release_table (struct scope_table *st, struct var **table,
                           unsigned n)
{
  unsigned i;

  if (!table)
    return;
  i = size_index (n);
  table[0] = (struct var *)st->free_tables[i];
  st->free_tables[i] = table;
}

static void grow (struct scope_table *st, struct scope *sc)
{
  unsigned n = sc->nbuckets ? sc->nbuckets * 2 : MIN_BUCKETS;
  struct var **bucket;

  if (n > MAX_BUCKETS)
    return;
  bucket = take_table (st, n);
  for (unsigned i = 0; i < sc->nbuckets; i++) {
    struct var *next;

    for (struct var *v = sc->bucket[i]; v; v = next) {
      unsigned b = v->hash & (n - 1);

      next = v->chain;
      v->chain = bucket[b];
      bucket[b] = v;
    }
  }
  release_table (st, sc->bucket, sc->nbuckets);
  sc->bucket = bucket;
  sc->nbuckets = n;
}

/* The scope's entries stay in the arena, so a struct var handed out earlier
 * stays valid. */
static void drop_scope_table (struct scope_table *st, struct scope *sc)
{
  release_table (st, sc->bucket, sc->nbuckets);
  sc->bucket = NULL;
  sc->nbuckets = 0;
  sc->n = 0;
}

struct scope_table *scope_table_new (void)
{
  struct scope_table *st = allocate (sizeof (*st));

  scope_push (st);
  for (size_t i = 0; i < sizeof (builtin_types) / sizeof (*builtin_types); i++)
    scope_add (st, builtin_types[i], SCOPE_TYPEDEF, NULL);
  return st;
}

struct var *scope_add (struct scope_table *st, const char *name,
                       enum scope_kind kind, const struct ast *decl)
{
  struct scope *sc = st->top;
  struct var *v;
  unsigned hash, b;

  if (!name)
    return NULL;
  hash = hash_name (name);
  v = find_in (sc, name, hash);
  if (!v) {
    if (sc->n >= sc->nbuckets)
      grow (st, sc);
    v = allocate (sizeof (*v));
    v->name = allocate_str (name, (int)strlen (name));
    v->hash = hash;
    b = hash & (sc->nbuckets - 1);
    v->chain = sc->bucket[b];
    sc->bucket[b] = v;
    sc->n++;
  }
  v->kind = kind;
  v->decl = decl;
  return v;
}

const struct var *scope_find (const struct scope_table *st, const char *name)
{
  unsigned hash;

  if (!name)
    return NULL;
  hash = hash_name (name);
  for (const struct scope *sc = st->top; sc; sc = sc->up) {
    const struct var *v = find_in (sc, name, hash);

    if (v)
      return v;
  }
  return NULL;
}

int scope_is_type (const struct scope_table *st, const char *name)
{
  const struct var *v = scope_find (st, name);

  return v && v->kind == SCOPE_TYPEDEF;
}

void scope_push (struct scope_table *st)
{
  struct scope *sc = allocate (sizeof (*sc));

  sc->up = st->top;
  st->top = sc;
}

void scope_pop (struct scope_table *st)
{
  struct scope *sc = st->top;

  if (!sc->up)
    return;
  st->top = sc->up;
  drop_scope_table (st, sc);
}

void scope_push_params (struct scope_table *st)
{
  scope_push (st);
  st->param_depth++;
}

void scope_pop_params (struct scope_table *st)
{
  /* The body of a function reopens this scope, so it keeps its table. */
  if (st->param_depth > 0 && --st->param_depth == 0 && !st->fn_params) {
    st->fn_params = st->top;
    st->top = st->top->up;
    return;
  }
  scope_pop (st);
}

void scope_resume_params (struct scope_table *st)
{
  struct scope *sc = st->fn_params;

  st->fn_params = NULL;
  if (!sc) {
    scope_push (st);
    return;
  }
  sc->up = st->top;
  st->top = sc;
}

void scope_declarator_named (struct scope_table *st)
{
  /* The kept parameters belonged to a prototype: no body is coming for them. */
  if (st->param_depth == 0 && st->fn_params) {
    drop_scope_table (st, st->fn_params);
    st->fn_params = NULL;
  }
}

void scope_typedef_seen (struct scope_table *st)
{
  st->typedef_seen = 1;
}

void scope_specifiers_done (struct scope_table *st)
{
  st->declaring = st->typedef_seen ? SCOPE_TYPEDEF : SCOPE_VARIABLE;
  st->typedef_seen = 0;
}

void scope_specifiers_reset (struct scope_table *st)
{
  st->typedef_seen = 0;
}

enum scope_kind scope_declaring (const struct scope_table *st)
{
  return st->declaring;
}
