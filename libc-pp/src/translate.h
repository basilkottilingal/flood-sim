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

/* Rewrites the scalars that decl, a declaration in a function, declares:
 *
 *   scalar x;          becomes  scalar x __attribute__ ((cleanup (...))) = scalar_new ();
 *   static scalar s;   becomes  static scalar s, *const slot ... = &s;
 *
 * so that x gets an id from scalar_new (), and scalar_free (x) runs however
 * its scope ends, while s is numbered at startup like a global scalar.
 * Pointers, arrays, names with an initializer, and extern, register and
 * thread-local scalars are left alone. The parser calls this as it reduces a
 * declaration in a block. */
void local_scalars (struct ast *decl);

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

#endif
