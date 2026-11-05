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

/* The loop `foreach (v in V) body`, given its tokens and body, rewritten as
 *
 *   for (int __i_v = 0; __i_v < V.n; ++__i_v) { __auto_type v = V.a[__i_v]; body }
 *
 * The body is moved into the result, and the other arguments are dropped. Lines
 * are kept: the result spans as many as the loop did. */
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

struct ast *foreach_to_for (struct ast *kw, struct ast *lparen, struct ast *var,
                            struct ast *in, struct ast *coll, struct ast *rparen,
                            struct ast *body);

#endif
