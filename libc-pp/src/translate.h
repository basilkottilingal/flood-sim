#ifndef TRANSLATE_H
#define TRANSLATE_H

enum {
  TRANSLATE_OK,
  TRANSLATE_FAILED,       /* I/O error, or an error reported: no gcc after this */
  TRANSLATE_SKIPPED,      /* did not parse, so copied as it is for gcc to judge */
};

/* Translates preprocessed C from in_path to out_path. Returns a TRANSLATE_ value. */
int translate_file (const char *in_path, const char *out_path);

struct ast;
struct scope_table;

/* Warns if target, the operand of an assignment, of ++ or of --, names a
 * variable declared as a scalar, whose value cwrap gives it. op is the operator.
 * The parser calls this, as it knows which declaration of a name is visible. */
void check_scalar_assignment (const struct scope_table *scopes,
                              const struct ast *target, const struct ast *op);

/* The euclidean dimension, 1 to 3, from -dimension=; 0 turns scalars off. */
void translate_set_dimension (int dimension);

/* Rewrites the scalars that decl, a declaration in a function, declares:
 *
 *   scalar x;          becomes  scalar x __attribute__ ((cleanup (...))) = scalar_new ();
 *   static scalar s;   becomes  static scalar s; static scalar *const slot ... = &s;
 *
 * so that x gets an id from scalar_new (), and scalar_free (x) runs however
 * its scope ends, while s is numbered at startup like a global scalar. A
 * vector or tensor gets the same for each of its scalars. Pointers, arrays,
 * names with an initializer, and extern, register and thread-local scalars are
 * left alone. The parser calls this as it reduces a declaration in a block,
 * and uses what it returns: decl, or decl with its slots after it. */
struct ast *local_scalars (struct ast *decl);

/* The loop `foreach (v in V) body`, given its tokens and body, rewritten as
 *
 *   for (int __i_v = 0; __i_v < V.n; ++__i_v) { __auto_type v = V.a[__i_v]; body }
 *
 * If V is a cache, the index starts at
 *
 *   V.update && V.is_dirty ? (V.is_dirty = 0, V.update (&V), 0) : 0
 *
 * The body is moved into the result, and the other arguments are dropped. Lines
 * are kept: the result spans as many as the loop did. */
struct ast *foreach_to_for (struct ast *kw, struct ast *lparen, struct ast *var,
                            struct ast *in, struct ast *coll, struct ast *rparen,
                            struct ast *body, int is_cache);

/* Caches, from <cwrap/cache.h>. The functions that return an int return the
 * number of errors they reported: the parser adds them up, and an error keeps
 * the file from gcc. */

/* cache (T), given its tokens, as the struct specifier `struct cache_T`. The
 * translation defines the struct before the first external declaration that
 * uses it. */
struct ast *cache_specifier (struct ast *kw, struct ast *lparen, struct ast *type,
                             struct ast *rparen, int *rejected);

/* Gives each plain cache that decl, a declaration in a function, declares
 * `= { 0 }`, unless it has a brace initializer, and a cleanup that frees its
 * array. Any other initializer is an error, as it would share an array. */
int local_caches (struct ast *decl);

/* Where the name of a cache variable is expected: the declaration of the
 * cache that expr names, or NULL if it is not one. */
const struct ast *cache_variable (const struct scope_table *scopes,
                                  const struct ast *expr);

/* An error if expr names a cache, which at is about to copy: the target or
 * the value of an assignment, or the value of a return. */
int check_cache_copy (const struct scope_table *scopes, const struct ast *expr,
                      const struct ast *at);

/* A warning for each cache passed by value in args, an argument list. */
void check_cache_arguments (const struct scope_table *scopes,
                            const struct ast *args);

/* An error for a typedef of cache (T): caches are known by cache (T) alone. */
int check_cache_declaration (const struct ast *decl);

/* Block macros.
 *
 *   @loop name (T1 p1, T2 p2) : a, b
 *   {
 *     ... {...} ...
 *   }
 *
 * defines a loop; `name (e1, e2) U`, U a statement as after for, then becomes
 *
 *   { T1 p1' = e1; T2 p2' = e2; <the body, with {...} replaced by { U }> }
 *
 * The parameters are typed locals, so each argument is converted and
 * evaluated once. The body's own names, parameters included, are renamed
 * where they are declared and used, so U cannot see them nor they hide U's
 * names, except those listed after `:`, which must be declared where {...}
 * is: those are what U sees. The definition itself emits nothing.
 *
 * The parser keeps the definitions, and checks {...} and return as it goes. */

struct var;

struct loop_inject {
  const char *name;
  const struct var *var;  /* its declaration in the body */
};

struct loop_def {
  const char *name;
  int serial;             /* the loop of struct var: see scope.h */
  const struct token *at; /* the name, for messages */
  struct ast *source;     /* LOOP_SOURCE: every token of the definition */
  struct ast *params;     /* PARAMETER_LIST or PARAMETER_TYPE_LIST, or NULL */
  struct ast *body;       /* COMPOUND_STATEMENT */
  struct ast *injects;    /* LOOP_INJECTS, or NULL */
  struct loop_inject *inject;
  int ninject;
  int placeholders;       /* {...} seen in the body */
  int complete;           /* parsed to its end, and usable */
  struct loop_def *next;
};

/* What goes in the tree for a definition: an empty token with the trivia of
 * all of its tokens, so that no line moves. */
struct ast *loop_stub (const struct ast *source);

/* Checks the parameters of def, once parsed: each is named, and there is no
 * `...`. Returns the number of errors reported. */
int loop_check_params (const struct loop_def *def);

/* Whether two definitions have the same tokens. */
int loop_same (const struct loop_def *a, const struct loop_def *b);

/* name (args) block, a use of def: args is an ARGUMENT_LIST or NULL, block
 * the statement U. given[i] is the declaration that the block saw for the
 * name def->inject[i]: the expansion's own declaration of it takes its place,
 * so that in a definition the two stay one name, renamed together or not.
 * Returns the expansion; on an error, reported and counted in *rejected,
 * returns block. */
struct ast *loop_expand (const struct loop_def *def, struct ast *name,
                         struct ast *args, struct ast *block,
                         const struct var *const *given, int *rejected);

#endif
