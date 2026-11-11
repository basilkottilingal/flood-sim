/* The names declared so far, by scope. C cannot be parsed without them, since
 * `T * x;` is a declaration only if T is a type: the lexer asks scope_is_type ()
 * about each identifier, and the parser adds every name it declares, typedef or
 * not, because an ordinary name can hide a typedef in an inner scope.
 *
 * Nothing here reports an error or stops the run. A redeclaration, a conflict
 * of kinds or an unbalanced scope_pop () is accepted as it comes, and judging
 * the code is left to gcc, which compiles the output. All of it lives in the
 * arena, so free_arena () releases it. */
#ifndef SCOPE_H
#define SCOPE_H

struct ast;
struct scope_table;

enum scope_kind {
  SCOPE_VARIABLE,         /* any ordinary name: variable, function,
                           * parameter or enumerator */
  SCOPE_TYPEDEF,
  SCOPE_LOOP,             /* a name defined by @loop */
  SCOPE_STENCIL,          /* a name declared by @foreach */
};

struct var {
  const char *name;
  enum scope_kind kind;
  const struct ast *decl; /* where it was declared, NULL for gcc's builtins */
  unsigned hash;
  int loop;               /* declared in the @loop definition with this
                           * serial, as a parameter or a local; else 0 */
  struct var *chain;      /* next in the hash bucket */
};

struct scope_table *scope_table_new (void);

/* Declares name in the current scope, replacing an earlier declaration of it
 * there. Returns the entry, or NULL if name is NULL. The entry stays valid after
 * its scope closes. */
struct var *scope_add (struct scope_table *st, const char *name,
                       enum scope_kind kind, const struct ast *decl);

/* The declaration of name visible now, innermost scope first, or NULL. */
const struct var *scope_find (const struct scope_table *st, const char *name);
int scope_is_type (const struct scope_table *st, const char *name);

/* Block and for-statement scopes. Popping the file scope does nothing. */
void scope_push (struct scope_table *st);
void scope_pop (struct scope_table *st);

/* Function parameters get a scope of their own, which closes at the `)` and
 * reopens for the body with scope_resume_params (). The scope kept is that of
 * the declarator's first parameter list after its name, so in
 * `int (*f (int a)) (int b) { ... }` the body sees a, not b. */
void scope_push_params (struct scope_table *st);
void scope_pop_params (struct scope_table *st);
void scope_resume_params (struct scope_table *st);

/* For the parser's error recovery: closes a parameter scope whose function
 * will not be parsed, and clears what a declaration left half done. */
void scope_abandon_params (struct scope_table *st);
void scope_error_recovered (struct scope_table *st);
void scope_declarator_named (struct scope_table *st);

/* Is the declaration being parsed a typedef? The parser reports the keyword as
 * it passes, then latches it once the specifiers end, so that a declaration
 * nested in a struct body cannot disturb it. */
void scope_typedef_seen (struct scope_table *st);
void scope_specifiers_done (struct scope_table *st);
void scope_specifiers_reset (struct scope_table *st);
enum scope_kind scope_declaring (const struct scope_table *st);

#endif
