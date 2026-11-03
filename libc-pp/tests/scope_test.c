/* Unit checks of scope.c: lookups, shadowing, table reuse, and that nothing it
 * is handed can stop it. Exits non-zero on the first failure.
 *
 * Includes scope.c itself, to see the free lists without widening its API. */
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "allocator.h"
#include "scope.c"

static void test_lookup_and_shadowing (void)
{
  struct scope_table *st = scope_table_new ();

  assert (!scope_find (st, "x"));
  assert (scope_is_type (st, "__builtin_va_list"));

  scope_add (st, "T", SCOPE_TYPEDEF, NULL);
  assert (scope_is_type (st, "T"));

  scope_push (st);
  assert (scope_is_type (st, "T"));
  scope_add (st, "T", SCOPE_VARIABLE, NULL);
  assert (!scope_is_type (st, "T"));
  assert (scope_find (st, "T")->kind == SCOPE_VARIABLE);
  scope_pop (st);

  assert (scope_is_type (st, "T"));
}

/* Each of these is an error in C, for gcc to report and not for us. */
static void test_bad_input_is_accepted (void)
{
  struct scope_table *st = scope_table_new ();

  scope_add (st, "a", SCOPE_VARIABLE, NULL);
  assert (scope_add (st, "a", SCOPE_TYPEDEF, NULL));
  assert (scope_is_type (st, "a"));
  assert (!scope_add (st, NULL, SCOPE_VARIABLE, NULL));
  assert (!scope_find (st, NULL));

  for (int i = 0; i < 5; i++) {
    scope_pop (st);
    scope_pop_params (st);
  }
  assert (scope_is_type (st, "a"));
  scope_add (st, "b", SCOPE_VARIABLE, NULL);
  assert (scope_find (st, "b"));
}

static void test_function_parameters (void)
{
  struct scope_table *st = scope_table_new ();

  /* int (*f (int a)) (int b) { ... } */
  scope_declarator_named (st);
  scope_push_params (st);
  scope_add (st, "a", SCOPE_VARIABLE, NULL);
  scope_pop_params (st);
  assert (!scope_find (st, "a"));
  scope_push_params (st);
  scope_add (st, "b", SCOPE_VARIABLE, NULL);
  scope_pop_params (st);
  assert (!scope_find (st, "b"));

  scope_resume_params (st);
  assert (scope_find (st, "a"));
  assert (!scope_find (st, "b"));
  scope_pop (st);
  assert (!scope_find (st, "a"));
}

static void test_closed_scope_table_is_reused (void)
{
  struct scope_table *st = scope_table_new ();
  struct var **table;

  scope_push (st);
  scope_add (st, "old", SCOPE_VARIABLE, NULL);
  table = st->top->bucket;
  scope_pop (st);
  assert (st->free_tables[0] == table);

  scope_push (st);
  scope_add (st, "new", SCOPE_VARIABLE, NULL);
  assert (st->top->bucket == table);
  assert (st->free_tables[0] == NULL);
  assert (scope_find (st, "new"));
  assert (!scope_find (st, "old"));     /* the table came back empty */
  scope_pop (st);
}

/* An entry outlives its scope's table, for whoever kept a pointer to it. */
static void test_entry_survives_its_table (void)
{
  struct scope_table *st = scope_table_new ();
  const struct var *v;

  scope_push (st);
  v = scope_add (st, "kept", SCOPE_TYPEDEF, NULL);
  scope_pop (st);
  scope_push (st);
  scope_add (st, "other", SCOPE_VARIABLE, NULL);
  assert (strcmp (v->name, "kept") == 0 && v->kind == SCOPE_TYPEDEF);
  scope_pop (st);
}

static void test_outgrown_table_is_reused (void)
{
  struct scope_table *st = scope_table_new ();
  char name[32];

  scope_push (st);
  for (int i = 0; i <= MIN_BUCKETS; i++) {
    snprintf (name, sizeof (name), "v%d", i);
    scope_add (st, name, SCOPE_VARIABLE, NULL);
  }
  assert (st->top->nbuckets == MIN_BUCKETS * 2);
  assert (st->free_tables[0] != NULL);
  for (int i = 0; i <= MIN_BUCKETS; i++) {
    snprintf (name, sizeof (name), "v%d", i);
    assert (scope_find (st, name));
  }
  scope_pop (st);
  assert (st->free_tables[1] != NULL);
}

/* A table waiting for the function body must not go to another scope. */
static void test_parameters_keep_their_table (void)
{
  struct scope_table *st = scope_table_new ();
  struct var **table;

  scope_declarator_named (st);
  scope_push_params (st);
  scope_add (st, "a", SCOPE_VARIABLE, NULL);
  table = st->top->bucket;
  scope_pop_params (st);
  assert (st->free_tables[0] != table);

  scope_push (st);
  scope_add (st, "other", SCOPE_VARIABLE, NULL);
  assert (st->top->bucket != table);
  scope_pop (st);

  scope_resume_params (st);
  assert (scope_find (st, "a"));
  scope_pop (st);
  assert (st->free_tables[0] == table);
  assert (!scope_find (st, "a"));
}

/* After a prototype, the next declarator gives its parameters' table back. */
static void test_prototype_gives_its_table_back (void)
{
  struct scope_table *st = scope_table_new ();
  struct var **table;

  scope_declarator_named (st);
  scope_push_params (st);
  scope_add (st, "a", SCOPE_VARIABLE, NULL);
  table = st->top->bucket;
  scope_pop_params (st);
  assert (st->free_tables[0] != table);
  scope_declarator_named (st);
  assert (st->free_tables[0] == table);
}

/* Past the largest table, which only lengthens the chains. */
static void test_many_names (void)
{
  struct scope_table *st = scope_table_new ();
  enum { N = 200000 };
  char name[32];

  for (int i = 0; i < N; i++) {
    snprintf (name, sizeof (name), "name%d", i);
    scope_add (st, name, i % 2 ? SCOPE_TYPEDEF : SCOPE_VARIABLE, NULL);
  }
  for (int i = 0; i < N; i++) {
    snprintf (name, sizeof (name), "name%d", i);
    assert (scope_find (st, name));
    assert (scope_is_type (st, name) == (i % 2));
  }
}

int main (void)
{
  test_lookup_and_shadowing ();
  test_bad_input_is_accepted ();
  test_function_parameters ();
  test_closed_scope_table_is_reused ();
  test_entry_survives_its_table ();
  test_outgrown_table_is_reused ();
  test_parameters_keep_their_table ();
  test_prototype_gives_its_table_back ();
  test_many_names ();
  free_arena ();
  return 0;
}
