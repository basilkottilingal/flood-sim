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

/* Rewrites the scalars that decl, a declaration in function, declares:
 *
 *   scalar x;          becomes  scalar x __attribute__ ((cleanup (...)))
 *                                 = __cwrap_scalar_new_named ("x");
 *   static scalar s;   becomes  static scalar s; static const struct
 *                                 __cwrap_scalar_slot slot ... = { &s, "f.s" };
 *
 * Each is named after its declaration, a string literal, so nothing is
 * copied or allocated: "x", or "f.x" if it hides another scalar of that name
 * (a global, or a local of an outer block: then "f.x.2", "f.x.3", ...), so
 * that two scalars alive in one scope never share a name; a static one is
 * always "f.s". A vector's are "u.x", "u.y", .... function may be a @loop's
 * name, in its definition. Only these are named: any other scalar, as one
 * from scalar_new (), has the name "".
 *
 * so that x gets an id from scalar_new (), and scalar_free (x) runs however
 * its scope ends, while s is numbered at startup like a global scalar. A
 * vector or tensor gets the same for each of its scalars. Pointers, arrays,
 * names with an initializer, and extern, register and thread-local scalars are
 * left alone. The parser calls this as it reduces a declaration in a block,
 * and uses what it returns: decl, or decl with its slots after it. */
struct ast *local_scalars (struct ast *decl, const char *function);

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

/* An error if a call, fn (args), is scalar_free, vector_free or tensor_free
 * of a global scalar, vector or tensor (at file scope, static or extern), or
 * of a member of one: a global scalar is numbered once for the whole run and
 * never freed, as other scalars' operators may depend on it. Returns the
 * number of errors reported. */
int check_global_free (const struct ast *fn, const struct ast *args);

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

/* Stencils, with -dimension=N: the cells of every leaf block of libb-amr.
 *
 *   foreach ()
 *     s[] = (t[1] + t[-1] + t[0,1] + t[0,-1])/4.;
 *
 * s[] is the cell of scalar s, and s[i], s[i,j], s[i,j,k] the cell at that
 * offset from it: offsets left out are 0. The loop is two loops, the @loops
 * that <block-loops.h> defines:
 *
 *   foreach_leaf () {                          every leaf block
 *     Real *const       __cwrap_block_L_0 = block_data (sblocks, s);
 *     const Real *const __cwrap_block_L_1 = block_data (sblocks, t);
 *     foreach_cell (tree_block_shape (), 0)    every interior cell
 *       body, with s[i,j] as __cwrap_block_L_0[c + (i) + (j) * shape.sy]
 *   }
 *
 * Only the scalars that the body reads or writes get a block, once each, and
 * those it only reads a const one. The names the two loops declare (sblocks,
 * c, shape, i, j, level, ...) are renamed, so the body sees none of them.
 *
 * x, y, z and Delta, unless the body declares them, are the centre and size
 * of the cell, as in Basilisk, from the leaf's i, j, k and level and the
 * cell's ci, cj, ck; only those the body uses are declared. */

/* A stencil loop, declared with @foreach:
 *
 *   @foreach foreach () = foreach_leaf () foreach_cell (tree_block_shape (), 0);
 *
 * the loop over the blocks and the loop over the cells of each, two @loops,
 * with their arguments, which may use the parameters, as in
 * `@foreach name (int l) = ...`. A use, name (args) body, binds each
 * parameter to its argument, once, then expands as foreach () does. The
 * cell loop is to visit interior cells: the ghosts are kept up to date
 * around it, not by it. */
struct stencil_def {
  const char *name;
  const struct token *at;       /* the name, for messages */
  struct ast *source;           /* every token, to compare two definitions */
  struct ast *params;           /* PARAMETER_LIST, or NULL */
  const struct loop_def *leaf, *cell;
  struct ast *leaf_args, *cell_args;  /* ARGUMENT_LIST, or NULL */
  struct stencil_def *next;
};

/* @foreach, parsed: added to *list, unless defined before the same way;
 * an error if differently, or if leaf or cell is not a @loop. Returns what
 * goes in the tree: the definition's trivia. */
struct ast *stencil_define (struct stencil_def **list, struct ast *source,
                            struct ast *name, struct ast *params,
                            const struct loop_def *leaf,
                            struct ast *leaf_name, struct ast *leaf_args,
                            const struct loop_def *cell, struct ast *cell_name,
                            struct ast *cell_args, int *rejected);

struct stencil_loop {
  struct stencil_loop *outer;   /* the loop around this one, or NULL */
  const struct token *at;       /* foreach, for messages */
  const char *name;             /* foreach, or another @foreach */
  struct ast **access;          /* every s[...] of the body, in order */
  int naccess, cap;
  int errors;                   /* reported before it: those after are its */
  struct ast **mover;           /* its foreach_neighbor (), in order */
  int nmover, mcap;
  int in_mover;                 /* in one now */
  const char *mover_name;       /* its name */
};

/* foreach_neighbor (r) statement, in a stencil loop: the statement runs for
 * the point moved to each cell within r of it along every axis, as in
 * Basilisk. s[], s[i,j] and x, y, z, Delta are then those of the moved
 * point. It only reads: writing a cell in it is an error, as it would write
 * a neighbour; its reads are reads at an offset, so the lazy ghost update
 * applies. r, as every offset of s[...], must be an integer constant within
 * BAMR_GHOSTS (<block-loops.h>), the offsets inside adding r: gcc checks it
 * with a static assertion, so no access leaves the block. Outside a stencil
 * loop, or in
 * another foreach_neighbor (), it is an error. */
void stencil_mover_begin (struct stencil_loop *loop, struct ast *kw,
                          int *rejected);

/* Loops over two levels, as libb-amr's foreach_prolongation (pr) and
 * foreach_restriction (p), are told apart by the names their two loops give,
 * as sblocks, c and shape are:
 *
 *   coarse_sblocks (block loop), cc and sgn (cell loop): the point has a
 *     coarse cell, cc in the coarse blocks. coarse (s), coarse (s, i, j)
 *     read it and its neighbours at that level, |i| at most
 *     BAMR_COARSE_REACH; child.x, .y, .z is sgn, the half the point is in.
 *   child_sblocks (block loop) and fc (cell loop): the point has children,
 *     the first at fc in the child's blocks. foreach_child () { ... } moves
 *     the point to each, as foreach_neighbor () does, and only reads.
 *
 * In such a loop, an offset on the point's own level (s[1]) and
 * foreach_neighbor () are errors, and no ghost is tracked: it runs inside
 * the ghost update itself. coarse () only reads.
 *
 * A loop whose block loop gives ghosts_by_caller (libb-amr's foreach_level
 * (l), on a level of the tree, not its leaves) has every rule of foreach (),
 * offsets and their reach checked, a write at an offset an error and a read
 * of a neighbour the loop writes a warning; but no ghost update or mark is
 * written around it: its caller updates the ghosts of that level, and marks
 * what it wrote.
 *
 * A loop whose cell loop gives checkerboard visits the cells of one colour
 * (libb-amr's foreach_level_or_leaf_colour (l, colour), red-black
 * Gauss-Seidel), or every cell in a fixed order. A scalar it writes may be
 * read at an offset whose sum is odd, the other colour: gcc checks it with
 * a static assertion, as the reach, so a read of its own colour is a
 * compile error; and none of the warning above. foreach_neighbor ()
 * reading such a scalar is an error: it reaches both colours. */
struct ast *stencil_coarse (struct stencil_loop *loop, struct ast *fn,
                            struct ast *open, struct ast *args,
                            struct ast *close, int *rejected);
struct ast *stencil_mover_end (struct stencil_loop *loop, struct ast *kw,
                               struct ast *open, struct ast *radius,
                               struct ast *close, struct ast *body);

/* Ghost operators: prolongations, restrictions and boundary conditions,
 * chosen per scalar by an integer id and called through a switch, with no
 * function pointer:
 *
 *   @prolongation name (params) { body }   at file scope: the operator,
 *       static inline void __cwrap_prolongation_name (params) { body }
 *       const int name = id;
 *     ids count from 1, per kind, in the order of definition in the file;
 *     0 is left for "the default". A name defined twice is an error.
 *   @prolongation name;                     extern const int name;
 *   @prolongation (id) (args);              a statement:
 *       switch (id) { case 1: __cwrap_prolongation_a (args); break; ...
 *                     default: __builtin_trap (); }
 *     over every operator of that kind defined before it in the file.
 *
 * and the same with @restriction and @boundary. So the one file that
 * dispatches (libb-amr's block-ghosts.c, or a program's own) must see every
 * definition; any other file attaches an operator by name, an int. */
struct ast *operator_define (struct ast *at, struct ast *name, struct ast *open,
                             struct ast *params, struct ast *close,
                             struct ast *deps, struct ast *body, int *rejected);
struct ast *operator_declare (struct ast *at, struct ast *name, struct ast *semi);
struct ast *operator_dispatch (struct ast *at, struct ast *open, struct ast *id,
                               struct ast *close, struct ast *open2,
                               struct ast *args, struct ast *close2,
                               struct ast *semi, int *rejected);

/* Coupled operators. A second parameter list, after the first, holds the
 * scalars a definition depends on, bound per scalar:
 *
 *   @prolongation lake (const Prolongation * pr) (scalar zb) { ... }
 *       static inline void __cwrap_prolongation_lake (const Prolongation * pr,
 *                                                     scalar zb) { ... }
 *   @prolongation (id, uses) (args);  as above, passing uses [0], uses [1],
 *       ... to an operator that takes them; without uses, an error if one
 *       of its kind does
 *   @prolongation uses (id)           an expression: how many scalars the
 *       operator with that id takes, -1 for no such id
 *   ghost_ops (h, prolongation (lake, zb), boundary (bc_odd))
 *       ghost_ops_set (h, GHOST_PROLONGATION, lake, (const scalar []) {zb}, 1),
 *       ghost_ops_set (h, GHOST_BOUNDARY, bc_odd, (const scalar *) 0, 0)
 *     each scalar depended on a global one (at file scope, static or
 *     extern, or a member of a global vector or tensor), which is never
 *     freed: else an error. */
struct ast *operator_uses (struct ast *at, struct ast *word, struct ast *open,
                           struct ast *id, struct ast *close, int *rejected);
struct ast *ghost_ops_call (struct ast *fn, struct ast *args, int *rejected);

/* Whether `foreach ()` is the stencil loop: -dimension= is given. */
int stencil_enabled (void);
int stencil_dimension (void);       /* -dimension=, or 0 */

/* foreach (, once `)` is read: the loop whose accesses the body adds to.
 * One inside another is an error. */
struct stencil_loop *stencil_loop_begin (struct stencil_loop *outer,
                                         struct ast *kw, int *rejected);

/* base[offsets] as the parser reduces it, offsets NULL for base[]. If base is
 * a scalar, returns the access, added to loop. Otherwise returns NULL for
 * base[offsets], an ordinary subscript, and for base[] reports an error and
 * returns a subscript. A scalar outside foreach () is an error too. */
/* foreach_dimension () statement, with -dimension=: the statement once for
 * each axis, as in Basilisk. The first as written; in the next, x becomes y,
 * y z, z x (in 2 dimensions x and y swap), and so on: the names x, y, z (the
 * cell's coordinates, a member .x), the names ending in _x, _y, _z, and the
 * offsets of s[i, j, k] and coarse (s, i, j, k), permuted the same way:
 * s[1] becomes s[0,1]. In a stencil loop the copies' accesses are the
 * loop's, as the statement's. The copies follow it on its last line: no
 * line moves. In 1 dimension, the statement as it is. */
struct ast *foreach_dimension (struct stencil_loop *loop, struct ast *kw,
                               struct ast *open, struct ast *close,
                               struct ast *body);

struct ast *stencil_access (struct stencil_loop *loop, struct ast *base,
                            struct ast *open, struct ast *offsets,
                            struct ast *close, int *rejected);

/* The loop with its body: the expansion described above. def is the
 * loop's @foreach, NULL if not declared: an error; args its arguments, an
 * ARGUMENT_LIST or NULL. has_ghosts: ghosts_mark_dirty (), ghosts_mark_all_dirty () and
 * ghosts_update_dirty () are declared, which a loop that writes or reads at
 * an offset needs: an error otherwise.
 *
 * Ghosts are kept up to date lazily, as in Basilisk:
 *
 *   { int __cwrap_stale_L[1], __cwrap_nstale_L = 0;      the scalars read at
 *     __cwrap_stale_L[__cwrap_nstale_L++] = (t);          an offset, updated
 *     ghosts_update_dirty (__cwrap_nstale_L, __cwrap_stale_L); }  if marked
 *   foreach_leaf () { ... }
 *   ghosts_mark_dirty (s);                                 those written
 *
 * A scalar A[k] that changes in the body stands for every element of A, if A
 * is an array of a known size, and any other for every scalar. Writing a
 * cell at an offset (s[1] = ...) is an error; reading at an offset a scalar
 * the loop writes, a warning. */
struct ast *stencil_loop_end (struct stencil_loop *loop,
                              const struct stencil_def *def, struct ast *args,
                              int has_ghosts, int has_located, int has_reach,
                              int has_coarse_reach, struct ast *kw,
                              struct ast *body, int *rejected);

/* A stencil loop by axis: N () statement, where N is not declared but N_x,
 * N_y (, N_z) are, @foreach each (libb-amr's foreach_face_x ()): N_x ()
 * statement, then N_y () with the statement rotated as foreach_dimension ()
 * rotates it, and so on, each loop with the rules of any. A loop whose cell
 * loop gives face_axis has x, y, z at the face: half a cell down that axis.
 * Before a loop, the scalars it reads at their cell are updated if stale
 * and not at the centre (ghosts_update_dirty_located (), if declared): a
 * face's own values change as fluxes are made to agree. */
struct ast *stencil_by_axis (struct stencil_loop *loop,
                             const struct stencil_def *const defs[3],
                             struct ast *args, int has_ghosts, int has_located,
                             int has_reach, int has_coarse_reach,
                             struct ast *kw, struct ast *body, int *rejected);

#endif
