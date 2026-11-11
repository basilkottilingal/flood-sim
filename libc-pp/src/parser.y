/* C11 with the GNU extensions that system headers rely on, as seen after
 * `gcc -E`. Every token is kept in the tree, so emission stays lossless.
 *
 * To add grammar:
 *   1. a new keyword: add a rule in lexer.l and a %token here
 *   2. a new node: add its kind to AST_KINDS in ast.h
 *   3. add the rule, building the node with N (KIND, $1, $2, ...)
 *
 * Every rule of more than one symbol needs an action that keeps all of its
 * symbols. Bison's default is $$ = $1, which silently loses the rest and shows
 * up as a mismatch in the byte-identical test. Pass NULL-valued optional
 * symbols to N () as they are: they are skipped. That includes the scope
 * markers, whose %destructor makes bison warn of a value left unused.
 *
 * %expect 0 turns any grammar conflict into a build error.
 *
 * Why two kinds of specifier list and of declarator: once a declaration has a
 * type, a TYPE_NAME must be the declared name (`typedef int T; { int T; }`), so
 * decl_specs_typed stops accepting one and decl_name accepts both. In a
 * parameter, `int f (int (T))` makes T a type, so an identifier in parentheses
 * must be separated from a type: hence notype_declarator. */
%code requires {
#include <stdio.h>

#include "lexer.h"

struct ast;
struct parse_status;

/* Returns the tree, or NULL if the input does not parse, after a warning for
 * each error found. The tree is in the arena, so it lasts until free_arena ().
 * *rejected is set if the input must not go on to gcc: the errors have then
 * been reported. */
struct ast *parse_file (FILE *in, int *rejected);
}

%code {
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "ast.h"
#include "diag.h"
#include "scope.h"
#include "translate.h"

/* Syntax errors to report before giving up on a file. */
#define MAX_PARSE_ERRORS 20

struct parse_status {
  int errors;             /* syntax errors reported */
  int rejected;           /* errors that must keep the file from gcc */
  struct loop_def *loops; /* @loop definitions, newest first */
  struct loop_def *cur_loop;  /* the one being parsed, or NULL */
  int loop_serial;
  struct stencil_loop *stencil;  /* the foreach () being parsed, or NULL */
  struct stencil_def *stencils;  /* @foreach definitions, newest first */
};

/* The action of an error rule, where the parse has found its feet again.
 * Bison reports no new error until three more tokens have parsed, which keeps
 * the leftovers of the part skipped from being reported as errors of their own. */
#define RECOVERED()                                                     \
  do {                                                                  \
    scope_error_recovered (scopes);                                     \
    if (status->errors >= MAX_PARSE_ERRORS)                             \
      YYABORT;                                                          \
  } while (0)

/* A declaration in a block, whose scalars and caches get what cwrap gives
 * them: decl, or decl followed by more items. */
#define LOCAL_DECLARATION(decl)                                         \
  (status->rejected += local_caches (decl), local_scalars (decl))

#define N(KIND, ...)                                                    \
  ast_node (AST_##KIND, (struct ast *[]) { __VA_ARGS__ },               \
            sizeof ((struct ast *[]) { __VA_ARGS__ }) / sizeof (struct ast *))

int yylex (YYSTYPE *lvalp, int *lineno, yyscan_t yyscanner);
static void yyerror (int *lineno, yyscan_t yyscanner, struct ast **root,
                     struct scope_table *scopes, struct parse_status *status,
                     const char *msg);

/* A name now in scope. In a @loop definition, it belongs to the loop, and
 * its token records the declaration, so the expansion can rename it. */
static void bind_name (struct parse_status *status, struct var *v,
                       struct ast *leaf)
{
  if (!v)
    return;
  v->loop = status->cur_loop ? status->cur_loop->serial : 0;
  if (leaf)
    ((struct ast_leaf *)leaf)->tok.binding = v;
}

/* The token of the name a declarator declares, or NULL. */
static struct ast *declared_leaf (const struct ast *d)
{
  if (d->kind == AST_DECL_NAME)
    return d->child;
  for (d = d->child; d; d = d->next)
    if (d->kind != AST_TOKEN) {
      struct ast *leaf = declared_leaf (d);

      if (leaf)
        return leaf;
    }
  return NULL;
}

static void declare (struct scope_table *scopes, struct parse_status *status,
                     const struct ast *declarator, enum scope_kind kind)
{
  bind_name (status, scope_add (scopes, ast_declared_name (declarator), kind,
                                declarator),
             declared_leaf (declarator));
}

static void loop_error (struct parse_status *status, const struct ast *at,
                        const char *format, const char *name)
{
  const struct token *t = ast_first_token (at);

  diag (t->file, t->line, DIAG_ERROR, format, name);
  status->rejected++;
}

static struct loop_def *find_loop (const struct parse_status *status,
                                   const char *name)
{
  for (struct loop_def *d = status->loops; d; d = d->next)
    if (d->complete && strcmp (d->name, name) == 0)
      return d;
  return NULL;
}

/* @loop name: the name is a loop from here, at file scope, so that a use of
 * it in its own body is caught. */
static void loop_begin (struct scope_table *scopes, struct parse_status *status,
                        struct ast *name)
{
  struct loop_def *def = allocate (sizeof (*def));

  def->name = ast_token (name)->text;
  def->at = ast_token (name);
  def->serial = ++status->loop_serial;
  status->cur_loop = def;
  scope_add (scopes, def->name, SCOPE_LOOP, name);
}

/* {...}: in a definition, once, with every name listed after `:` declared by
 * the loop and visible here. */
static struct ast *loop_placeholder (struct scope_table *scopes,
                                     struct parse_status *status,
                                     struct ast *open, struct ast *dots,
                                     struct ast *close)
{
  struct loop_def *def = status->cur_loop;
  struct ast *node = N (LOOP_BODY, open, dots, close);

  if (!def) {
    loop_error (status, open, "{...} outside a @loop definition%s", "");
    return node;
  }
  if (++def->placeholders > 1) {
    loop_error (status, open, "loop '%s' has more than one {...}", def->name);
    return node;
  }
  if (!def->injects)
    return node;
  struct ast *list = def->injects->last;
  int n = 0;

  for (struct ast *k = list->child; k; k = k->next)
    n += k->kind == AST_TOKEN && ast_token (k)->text[0] != ',';
  def->inject = allocate ((size_t)n * sizeof (*def->inject));
  for (struct ast *k = list->child; k; k = k->next) {
    const char *name = ast_token (k)->text;
    const struct var *v;

    if (name[0] == ',')
      continue;
    v = scope_find (scopes, name);
    if (!v || v->loop != def->serial) {
      loop_error (status, k, "'%s' is not declared by the loop where {...} is",
                  name);
      continue;
    }
    def->inject[def->ninject].name = name;
    def->inject[def->ninject++].var = v;
  }
  return node;
}

/* The end of a definition: checked, kept, and replaced in the tree by its
 * trivia. */
static struct ast *loop_end (struct parse_status *status, struct ast *source,
                             struct ast *params, struct ast *body)
{
  struct loop_def *def = status->cur_loop, *old;

  status->cur_loop = NULL;
  def->source = source;
  def->params = params;
  def->body = body;
  status->rejected += loop_check_params (def);
  if (!def->placeholders)
    loop_error (status, source, "loop '%s' has no {...} for the block",
                def->name);
  old = find_loop (status, def->name);
  if (old && !loop_same (old, def)) {
    loop_error (status, source, "loop '%s' is defined again, differently",
                def->name);
    diag (old->at->file, old->at->line, DIAG_NOTE, "'%s' is defined here",
          old->name);
  }
  if (!old) {
    def->complete = 1;
    def->next = status->loops;
    status->loops = def;
  }
  return loop_stub (source);
}

/* name ( ... ): the block will see the names the loop gives it. */
static void loop_open (struct scope_table *scopes, struct parse_status *status,
                       struct ast *name)
{
  const struct loop_def *def = find_loop (status, ast_token (name)->text);

  if (!def) {
    loop_error (status, name, "loop '%s' is used in its own definition",
                ast_token (name)->text);
    return;
  }
  /* In a definition, they are the outer loop's own: see loop_expand (). */
  for (int i = 0; i < def->ninject; i++)
    bind_name (status, scope_add (scopes, def->inject[i].name,
                                  def->inject[i].var->kind,
                                  def->inject[i].var->decl), NULL);
}

/* Expands the use while the names it gave the block are still in scope. */
static struct ast *loop_use (struct scope_table *scopes,
                             struct parse_status *status, struct ast *name,
                             struct ast *args, struct ast *block)
{
  const struct loop_def *def = find_loop (status, ast_token (name)->text);
  const struct var **given;

  if (!def)
    return block;
  given = allocate ((size_t)def->ninject * sizeof (*given) + 1);
  for (int i = 0; i < def->ninject; i++)
    given[i] = scope_find (scopes, def->inject[i].name);
  return loop_expand (def, name, args, block, given, &status->rejected);
}

static const struct stencil_def *find_stencil (const struct parse_status *status,
                                               const char *name)
{
  for (const struct stencil_def *d = status->stencils; d; d = d->next)
    if (strcmp (d->name, name) == 0)
      return d;
  return NULL;
}

/* The end of a stencil loop's body: the loop, expanded. */
static struct ast *stencil_use (struct scope_table *scopes,
                                struct parse_status *status, struct ast *name,
                                struct ast *args, struct ast *body)
{
  struct stencil_loop *loop = status->stencil;

  status->stencil = loop->outer;
  return stencil_loop_end (loop, find_stencil (status, ast_token (name)->text),
                           args,
                           scope_find (scopes, "ghosts_mark_dirty")
                           && scope_find (scopes, "ghosts_mark_all_dirty")
                           && scope_find (scopes, "ghosts_update_dirty"),
                           scope_find (scopes, "BAMR_GHOSTS") != NULL,
                           name, body, &status->rejected);
}

/* A file may have foreach loops or something named foreach, not both. Reports
 * the loop, and returns whether foreach is declared here. */
static int foreach_declared (const struct scope_table *scopes,
                             const struct ast *loop)
{
  const struct var *name = scope_find (scopes, "foreach");
  const struct token *t = ast_token (loop);

  if (!name)
    return 0;
  diag (t->file, t->line, DIAG_ERROR, "foreach loop, but 'foreach' is already "
        "declared");
  if (name->decl) {
    t = ast_first_token (name->decl);
    diag (t->file, t->line, DIAG_NOTE, "'foreach' is declared here; rename it "
          "to use the loop");
  }
  return 1;
}
}

%define api.pure full
%define api.value.type {struct ast *}
%define api.token.prefix {T_}
%define parse.error verbose
%param {int *lineno} {yyscan_t yyscanner}
%parse-param {struct ast **root} {struct scope_table *scopes}
%parse-param {struct parse_status *status}

%expect 0

/* Error recovery pops these markers without the rule that would undo them, so
 * their scopes are undone here instead. */
%destructor { scope_pop (scopes); } open_scope enter_function
%destructor { scope_abandon_params (scopes); } open_params

%token AUTO BREAK CASE CHAR CONST CONTINUE DEFAULT DO DOUBLE ELSE ENUM EXTERN
%token FLOAT FOR GOTO IF INLINE INT LONG REGISTER RESTRICT RETURN SHORT SIGNED
%token SIZEOF STATIC STRUCT SWITCH TYPEDEF UNION UNSIGNED VOID VOLATILE WHILE
%token ALIGNAS ALIGNOF ATOMIC BOOL COMPLEX GENERIC IMAGINARY NORETURN
%token STATIC_ASSERT THREAD_LOCAL
%token ATTRIBUTE EXTENSION ASM TYPEOF INT128 FLOATN REAL IMAG LABEL AUTO_TYPE
%token FOREACH CACHE
%token FOREACH_CELLS "foreach ()"
%token AT_LOOP "@loop" LOOP_NAME
%token AT_FOREACH "@foreach" STENCIL_NAME
%token FOREACH_NEIGHBOR "foreach_neighbor"
%token AT_OPERATOR "@prolongation, @restriction or @boundary"
%token BUILTIN_VA_ARG BUILTIN_OFFSETOF BUILTIN_TYPES_COMPATIBLE
%token IDENTIFIER TYPE_NAME CONSTANT STRING
%token ELLIPSIS SHL_ASSIGN SHR_ASSIGN ADD_ASSIGN SUB_ASSIGN MUL_ASSIGN
%token DIV_ASSIGN MOD_ASSIGN AND_ASSIGN OR_ASSIGN XOR_ASSIGN
%token INC DEC ARROW SHL SHR LE GE EQ NE LOGICAL_AND LOGICAL_OR
%token EOI "end of input"       /* carries the trailing trivia */

/* `if (a) if (b) x; else y;`: the else binds to the nearest if. */
%precedence THEN
%precedence ELSE
/* `_Atomic (int)` is a specifier, not the qualifier followed by a declarator. */
%precedence ATOMIC_QUALIFIER
/* `foreach (` is a loop or a call, and `cache (` a type or a call, never the
 * name in parentheses. */
%precedence KEYWORD_NAME
%precedence '('

%%

start:
    translation_unit EOI          { *root = ast_add ($1, $2); $$ = $1; }
  ;

/* External definitions */

translation_unit:
    %empty                        { $$ = N (TRANSLATION_UNIT); }
  | translation_unit external_declaration
                                  { $$ = ast_add ($1, $2); }
  ;

external_declaration:
    function_definition
  | EXTENSION function_definition { $$ = N (EXTENSION_DEFINITION, $1, $2); }
  | declaration
  | ASM '(' string ')' ';'        { $$ = N (ASM_DEFINITION, $1, $2, $3, $4, $5); }
  | ';'                           { $$ = N (EMPTY_DECLARATION, $1); }
  | loop_definition
  | stencil_definition
  | operator_definition
  /* Error recovery: skip to the end of the declaration, and go on to find more
   * errors. A file with errors is not translated, so these build nothing. */
  | error ';'                     { RECOVERED (); $$ = NULL; }
  | error '}'                     { RECOVERED (); $$ = NULL; }
  ;

function_definition:
    decl_specs_typed specs_done declarator enter_function compound_statement
    leave_function
                                  /* the name is in scope after the definition
                                   * too, as `foreach ()` must know */
                                  { declare (scopes, status, $3, SCOPE_VARIABLE);
                                    $$ = N (FUNCTION_DEFINITION, $1, $3, $4, $5); }
  ;

/* A block macro: see translate.h. The parameters and the body share a scope,
 * and every name declared in it belongs to the loop. */
loop_definition:
    AT_LOOP any_identifier        { loop_begin (scopes, status, $2); $$ = NULL; }
    '(' open_scope loop_params ')' loop_injects compound_statement close_scope
                                  { $$ = loop_end (status,
                                                   N (LOOP_SOURCE, $1, $2, $3, $4, $5,
                                                      $6, $7, $8, $9, $10),
                                                   $6, $9); }
  ;

/* A stencil loop: see translate.h. Its parameters are in scope for the
 * arguments of its two loops. foreach itself is not declared as a name, as
 * a declared foreach would be a function. */
stencil_definition:
    AT_FOREACH stencil_name '(' open_scope loop_params ')' '='
    LOOP_NAME '(' opt_arguments ')' LOOP_NAME '(' opt_arguments ')' ';' close_scope
                                  { $$ = stencil_define (&status->stencils,
                                                         N (LOOP_SOURCE, $1, $2, $3, $5, $6, $7,
                                                            $8, $9, $10, $11, $12, $13, $14,
                                                            $15, $16),
                                                         $2, $5,
                                                         find_loop (status, ast_token ($8)->text),
                                                         $8, $10,
                                                         find_loop (status, ast_token ($12)->text),
                                                         $12, $14, &status->rejected);
                                    if (strcmp (ast_token ($2)->text, "foreach") != 0)
                                      scope_add (scopes, ast_token ($2)->text, SCOPE_STENCIL, $2);
                                    (void) $4; (void) $17; }
  ;

stencil_name:
    any_identifier
  | FOREACH_CELLS
  ;

/* A ghost operator: see translate.h. Its parameters are in scope for its
 * body, as a function's. */
operator_definition:
    AT_OPERATOR any_identifier '(' open_scope loop_params ')' compound_statement
    close_scope
                                  { $$ = operator_define ($1, $2, $3, $5, $6, $7,
                                                          &status->rejected);
                                    (void) $4; (void) $8; }
  | AT_OPERATOR any_identifier ';'
                                  { $$ = operator_declare ($1, $2, $3); }
  ;

opt_arguments:
    %empty                        { $$ = NULL; }
  | argument_list
  ;

loop_params:
    %empty                        { $$ = NULL; }
  | parameter_type_list
  ;

loop_injects:
    %empty                        { $$ = NULL; }
  | ':' loop_inject_list          { $$ = N (LOOP_INJECTS, $1, $2);
                                    status->cur_loop->injects = $$; }
  ;

loop_inject_list:
    any_identifier                { $$ = N (LOOP_INJECT_LIST, $1); }
  | loop_inject_list ',' any_identifier
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

/* Markers: empty rules that run an action at that point of a rule. They reduce
 * before the next token is read, or with one that cannot be a name, so the
 * lexer sees their effect on the very next identifier. */
specs_done:      %empty           { scope_specifiers_done (scopes); $$ = NULL; } ;
enter_function:  %empty           { scope_resume_params (scopes); $$ = NULL; } ;
leave_function:  %empty           { scope_pop (scopes); $$ = NULL; } ;
open_scope:      %empty           { scope_push (scopes); $$ = NULL; } ;
close_scope:     %empty           { scope_pop (scopes); $$ = NULL; } ;
open_params:     %empty           { scope_push_params (scopes); $$ = NULL; } ;
close_params:    %empty           { scope_pop_params (scopes); $$ = NULL; } ;

/* Declarations */

declaration:
    decl_specs_typed ';'          { scope_specifiers_reset (scopes);
                                    $$ = N (DECLARATION, $1, $2); }
  | decl_specs_untyped ';'        { scope_specifiers_reset (scopes);
                                    $$ = N (DECLARATION, $1, $2); }
  | decl_specs_typed specs_done init_declarator_list ';'
                                  { $$ = N (DECLARATION, $1, $3, $4);
                                    status->rejected += check_cache_declaration ($$); }
  | EXTENSION declaration         { $$ = N (EXTENSION_DECLARATION, $1, $2); }
  | static_assert_declaration
  ;

static_assert_declaration:
    STATIC_ASSERT '(' constant_expression ')' ';'
                                  { $$ = N (STATIC_ASSERT, $1, $2, $3, $4, $5); }
  | STATIC_ASSERT '(' constant_expression ',' string ')' ';'
                                  { $$ = N (STATIC_ASSERT, $1, $2, $3, $4, $5, $6, $7); }
  ;

decl_specs_typed:
    type_specifier                { $$ = N (DECL_SPECS, $1); }
  | decl_specs_untyped type_specifier
                                  { $$ = ast_add ($1, $2); }
  | decl_specs_typed decl_spec_other
                                  { $$ = ast_add ($1, $2); }
  | decl_specs_typed type_specifier_keyword
                                  { $$ = ast_add ($1, $2); }
  ;

decl_specs_untyped:
    decl_spec_other               { $$ = N (DECL_SPECS, $1); }
  | decl_specs_untyped decl_spec_other
                                  { $$ = ast_add ($1, $2); }
  ;

decl_spec_other:
    storage_class
  | type_qualifier
  | function_specifier
  | alignment_specifier
  | attribute
  ;

storage_class:
    TYPEDEF                       { scope_typedef_seen (scopes); $$ = $1; }
  | EXTERN
  | STATIC
  | THREAD_LOCAL
  | AUTO
  | REGISTER
  ;

function_specifier:
    INLINE
  | NORETURN
  ;

alignment_specifier:
    ALIGNAS '(' type_name ')'     { $$ = N (ALIGNMENT_SPECIFIER, $1, $2, $3, $4); }
  | ALIGNAS '(' constant_expression ')'
                                  { $$ = N (ALIGNMENT_SPECIFIER, $1, $2, $3, $4); }
  ;

type_specifier:
    type_specifier_keyword
  | TYPE_NAME
  | cache_specifier
  ;

/* cache (T): struct cache_T, which the translation defines. A complete type,
 * like a typedef name, so no other type specifier goes with it. */
cache_specifier:
    CACHE '(' type_name ')'       { $$ = cache_specifier ($1, $2, $3, $4,
                                                          &status->rejected); }
  ;

type_specifier_keyword:
    VOID
  | CHAR
  | SHORT
  | INT
  | LONG
  | FLOAT
  | DOUBLE
  | SIGNED
  | UNSIGNED
  | BOOL
  | COMPLEX
  | IMAGINARY
  | INT128
  | FLOATN
  | AUTO_TYPE
  | struct_specifier
  | enum_specifier
  | typeof_specifier
  | atomic_specifier
  ;

typeof_specifier:
    TYPEOF '(' expression ')'     { $$ = N (TYPEOF_SPECIFIER, $1, $2, $3, $4); }
  | TYPEOF '(' type_name ')'      { $$ = N (TYPEOF_SPECIFIER, $1, $2, $3, $4); }
  ;

atomic_specifier:
    ATOMIC '(' type_name ')'      { $$ = N (ATOMIC_SPECIFIER, $1, $2, $3, $4); }
  ;

type_qualifier:
    CONST
  | RESTRICT
  | VOLATILE
  | ATOMIC %prec ATOMIC_QUALIFIER
  ;

/* Attributes */

attribute:
    ATTRIBUTE '(' '(' attribute_list ')' ')'
                                  { $$ = N (ATTRIBUTE, $1, $2, $3, $4, $5, $6); }
  ;

attribute_list:
    attribute_item                { $$ = N (ATTRIBUTE_LIST, $1); }
  | attribute_list ',' attribute_item
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

attribute_item:
    %empty                        { $$ = NULL; }
  | attribute_name                { $$ = N (ATTRIBUTE_ITEM, $1); }
  | attribute_name '(' ')'        { $$ = N (ATTRIBUTE_ITEM, $1, $2, $3); }
  | attribute_name '(' argument_list ')'
                                  { $$ = N (ATTRIBUTE_ITEM, $1, $2, $3, $4); }
  ;

attribute_name:
    any_identifier
  | CONST
  ;

attributes:
    attribute                     { $$ = N (ATTRIBUTES, $1); }
  | attributes attribute          { $$ = ast_add ($1, $2); }
  ;

opt_attributes:
    %empty                        { $$ = NULL; }
  | attributes
  ;

/* An identifier where a typedef name means nothing special: tags, members,
 * labels, enumerators. FOREACH and CACHE are names here, in decl_name and in
 * primary_expression, as the lexer returns them for every use but a typedef
 * once <cwrap/cache.h> is in. */
any_identifier:
    IDENTIFIER
  | TYPE_NAME
  | FOREACH
  | CACHE
  | LOOP_NAME
  | STENCIL_NAME
  ;

/* Declarators */

init_declarator_list:
    init_declarator               { $$ = N (INIT_DECLARATOR_LIST, $1); }
  | init_declarator_list ',' init_declarator
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  | init_declarator_list ',' attributes init_declarator
                                  { $$ = ast_add (ast_add (ast_add ($1, $2), $3), $4); }
  ;

init_declarator:
    named_declarator
  | named_declarator '=' initializer
                                  { $$ = N (INIT_DECLARATOR, $1, $2, $3); }
  ;

/* The name is in scope from here, before its initializer. */
named_declarator:
    declarator                    { declare (scopes, status, $1, scope_declaring (scopes));
                                    $$ = $1; }
  | declarator declarator_attrs   { declare (scopes, status, $1, scope_declaring (scopes));
                                    $$ = N (ATTRIBUTED_DECLARATOR, $1, $2); }
  ;

declarator_attrs:
    asm_label
  | attributes
  | asm_label attributes          { $$ = N (DECLARATOR_ATTRS, $1, $2); }
  ;

asm_label:
    ASM '(' string ')'            { $$ = N (ASM_LABEL, $1, $2, $3, $4); }
  ;

declarator:
    direct_declarator
  | pointer direct_declarator     { $$ = N (DECLARATOR, $1, $2); }
  ;

direct_declarator:
    decl_name
  | '(' declarator ')'            { $$ = N (PAREN_DECLARATOR, $1, $2, $3); }
  | '(' attributes declarator ')' { $$ = N (PAREN_DECLARATOR, $1, $2, $3, $4); }
  | direct_declarator declarator_suffix
                                  { $$ = N (SUFFIXED_DECLARATOR, $1, $2); }
  ;

decl_name:
    IDENTIFIER                    { scope_declarator_named (scopes);
                                    $$ = N (DECL_NAME, $1); }
  | TYPE_NAME                     { scope_declarator_named (scopes);
                                    $$ = N (DECL_NAME, $1); }
  | FOREACH                       { scope_declarator_named (scopes);
                                    $$ = N (DECL_NAME, $1); }
  | CACHE                         { scope_declarator_named (scopes);
                                    $$ = N (DECL_NAME, $1); }
  ;

declarator_suffix:
    array_suffix
  | function_suffix
  ;

array_suffix:
    '[' ']'                       { $$ = N (ARRAY_SUFFIX, $1, $2); }
  | '[' '*' ']'                   { $$ = N (ARRAY_SUFFIX, $1, $2, $3); }
  | '[' assignment_expression ']' { $$ = N (ARRAY_SUFFIX, $1, $2, $3); }
  | '[' qualifier_list ']'        { $$ = N (ARRAY_SUFFIX, $1, $2, $3); }
  | '[' qualifier_list '*' ']'    { $$ = N (ARRAY_SUFFIX, $1, $2, $3, $4); }
  | '[' qualifier_list assignment_expression ']'
                                  { $$ = N (ARRAY_SUFFIX, $1, $2, $3, $4); }
  | '[' STATIC assignment_expression ']'
                                  { $$ = N (ARRAY_SUFFIX, $1, $2, $3, $4); }
  | '[' STATIC qualifier_list assignment_expression ']'
                                  { $$ = N (ARRAY_SUFFIX, $1, $2, $3, $4, $5); }
  | '[' qualifier_list STATIC assignment_expression ']'
                                  { $$ = N (ARRAY_SUFFIX, $1, $2, $3, $4, $5); }
  ;

function_suffix:
    '(' ')'                       { $$ = N (FUNCTION_SUFFIX, $1, $2); }
  | '(' open_params parameter_type_list close_params ')'
                                  { $$ = N (FUNCTION_SUFFIX, $1, $2, $3, $5); }
  ;

pointer:
    '*'                           { $$ = N (POINTER, $1); }
  | '*' qualifier_list            { $$ = N (POINTER, $1, $2); }
  | '*' pointer                   { $$ = N (POINTER, $1, $2); }
  | '*' qualifier_list pointer    { $$ = N (POINTER, $1, $2, $3); }
  ;

qualifier_list:
    qualifier                     { $$ = N (QUALIFIER_LIST, $1); }
  | qualifier_list qualifier      { $$ = ast_add ($1, $2); }
  ;

qualifier:
    type_qualifier
  | attribute
  ;

parameter_type_list:
    parameter_list
  | parameter_list ',' ELLIPSIS   { $$ = N (PARAMETER_TYPE_LIST, $1, $2, $3); }
  ;

parameter_list:
    parameter_declaration         { $$ = N (PARAMETER_LIST, $1); }
  | parameter_list ',' parameter_declaration
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

parameter_declaration:
    decl_specs_typed              { $$ = N (PARAMETER, $1); }
  | decl_specs_typed abstract_declarator
                                  { $$ = N (PARAMETER, $1, $2); }
  | decl_specs_typed param_declarator opt_attributes
                                  { declare (scopes, status, $2, SCOPE_VARIABLE);
                                    $$ = N (PARAMETER, $1, $2, $3); }
  ;

param_declarator:
    param_direct
  | pointer param_direct          { $$ = N (DECLARATOR, $1, $2); }
  ;

param_direct:
    decl_name
  | '(' notype_declarator ')'     { $$ = N (PAREN_DECLARATOR, $1, $2, $3); }
  | param_direct declarator_suffix
                                  { $$ = N (SUFFIXED_DECLARATOR, $1, $2); }
  ;

notype_declarator:
    notype_direct
  | pointer param_direct          { $$ = N (DECLARATOR, $1, $2); }
  ;

notype_direct:
    IDENTIFIER                    { $$ = N (DECL_NAME, $1); }
  | '(' notype_declarator ')'     { $$ = N (PAREN_DECLARATOR, $1, $2, $3); }
  | notype_direct declarator_suffix
                                  { $$ = N (SUFFIXED_DECLARATOR, $1, $2); }
  ;

type_name:
    spec_quals_typed              { $$ = N (TYPE_NAME, $1); }
  | spec_quals_typed abstract_declarator
                                  { $$ = N (TYPE_NAME, $1, $2); }
  ;

spec_quals_typed:
    type_specifier                { $$ = N (SPEC_QUALS, $1); }
  | spec_quals_untyped type_specifier
                                  { $$ = ast_add ($1, $2); }
  | spec_quals_typed spec_qual_other
                                  { $$ = ast_add ($1, $2); }
  | spec_quals_typed type_specifier_keyword
                                  { $$ = ast_add ($1, $2); }
  ;

spec_quals_untyped:
    spec_qual_other               { $$ = N (SPEC_QUALS, $1); }
  | spec_quals_untyped spec_qual_other
                                  { $$ = ast_add ($1, $2); }
  ;

spec_qual_other:
    type_qualifier
  | alignment_specifier
  | attribute
  ;

abstract_declarator:
    pointer
  | direct_abstract_declarator
  | pointer direct_abstract_declarator
                                  { $$ = N (ABSTRACT_DECLARATOR, $1, $2); }
  ;

direct_abstract_declarator:
    '(' abstract_declarator ')'   { $$ = N (PAREN_DECLARATOR, $1, $2, $3); }
  | array_suffix
  | '(' ')'                       { $$ = N (FUNCTION_SUFFIX, $1, $2); }
  | '(' open_params parameter_type_list close_params ')'
                                  { $$ = N (FUNCTION_SUFFIX, $1, $2, $3, $5); }
  | direct_abstract_declarator declarator_suffix
                                  { $$ = N (SUFFIXED_DECLARATOR, $1, $2); }
  ;

/* Structures, unions and enumerations */

struct_specifier:
    struct_or_union opt_attributes any_identifier '{' opt_struct_declarations '}'
                                  { $$ = N (STRUCT_SPECIFIER, $1, $2, $3, $4, $5, $6); }
  | struct_or_union opt_attributes '{' opt_struct_declarations '}'
                                  { $$ = N (STRUCT_SPECIFIER, $1, $2, $3, $4, $5); }
  | struct_or_union opt_attributes any_identifier
                                  { $$ = N (STRUCT_SPECIFIER, $1, $2, $3); }
  ;

struct_or_union:
    STRUCT
  | UNION
  ;

opt_struct_declarations:
    %empty                        { $$ = NULL; }
  | struct_declaration_list
  ;

struct_declaration_list:
    struct_declaration            { $$ = N (STRUCT_DECLARATION_LIST, $1); }
  | struct_declaration_list struct_declaration
                                  { $$ = ast_add ($1, $2); }
  ;

struct_declaration:
    spec_quals_typed ';'          { $$ = N (STRUCT_DECLARATION, $1, $2); }
  | spec_quals_typed struct_declarator_list ';'
                                  { $$ = N (STRUCT_DECLARATION, $1, $2, $3); }
  | EXTENSION struct_declaration  { $$ = N (STRUCT_DECLARATION, $1, $2); }
  | static_assert_declaration
  | ';'                           { $$ = N (STRUCT_DECLARATION, $1); }
  ;

struct_declarator_list:
    struct_declarator             { $$ = N (STRUCT_DECLARATOR_LIST, $1); }
  | struct_declarator_list ',' struct_declarator
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

struct_declarator:
    declarator opt_attributes     { $$ = N (STRUCT_DECLARATOR, $1, $2); }
  | ':' constant_expression opt_attributes
                                  { $$ = N (STRUCT_DECLARATOR, $1, $2, $3); }
  | declarator ':' constant_expression opt_attributes
                                  { $$ = N (STRUCT_DECLARATOR, $1, $2, $3, $4); }
  ;

enum_specifier:
    ENUM opt_attributes any_identifier '{' enumerator_list '}'
                                  { $$ = N (ENUM_SPECIFIER, $1, $2, $3, $4, $5, $6); }
  | ENUM opt_attributes any_identifier '{' enumerator_list ',' '}'
                                  { $$ = N (ENUM_SPECIFIER, $1, $2, $3, $4, $5, $6, $7); }
  | ENUM opt_attributes '{' enumerator_list '}'
                                  { $$ = N (ENUM_SPECIFIER, $1, $2, $3, $4, $5); }
  | ENUM opt_attributes '{' enumerator_list ',' '}'
                                  { $$ = N (ENUM_SPECIFIER, $1, $2, $3, $4, $5, $6); }
  | ENUM opt_attributes any_identifier
                                  { $$ = N (ENUM_SPECIFIER, $1, $2, $3); }
  ;

enumerator_list:
    enumerator                    { $$ = N (ENUMERATOR_LIST, $1); }
  | enumerator_list ',' enumerator
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

enumerator:
    enumerator_name opt_attributes
                                  { $$ = N (ENUMERATOR, $1, $2); }
  | enumerator_name opt_attributes '=' constant_expression
                                  { $$ = N (ENUMERATOR, $1, $2, $3, $4); }
  ;

enumerator_name:
    any_identifier                { scope_add (scopes, ast_token ($1)->text, SCOPE_VARIABLE, $1);
                                    $$ = $1; }
  ;

/* Initializers */

initializer:
    assignment_expression
  | brace_initializer
  ;

brace_initializer:
    '{' '}'                       { $$ = N (INITIALIZER_LIST_EXPR, $1, $2); }
  | '{' initializer_list '}'      { $$ = N (INITIALIZER_LIST_EXPR, $1, $2, $3); }
  | '{' initializer_list ',' '}'  { $$ = N (INITIALIZER_LIST_EXPR, $1, $2, $3, $4); }
  ;

initializer_list:
    initializer_item              { $$ = N (INITIALIZER_LIST, $1); }
  | initializer_list ',' initializer_item
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

initializer_item:
    initializer
  | designation '=' initializer   { $$ = N (INITIALIZER_ITEM, $1, $2, $3); }
  | designation initializer       { $$ = N (INITIALIZER_ITEM, $1, $2); }
  | IDENTIFIER ':' initializer    { $$ = N (OLD_DESIGNATOR, $1, $2, $3); }
  ;

designation:
    designator                    { $$ = N (DESIGNATION, $1); }
  | designation designator        { $$ = ast_add ($1, $2); }
  ;

designator:
    '[' constant_expression ']'   { $$ = N (DESIGNATOR, $1, $2, $3); }
  | '[' constant_expression ELLIPSIS constant_expression ']'
                                  { $$ = N (DESIGNATOR, $1, $2, $3, $4, $5); }
  | '.' any_identifier            { $$ = N (DESIGNATOR, $1, $2); }
  ;

/* Statements */

statement:
    labeled_statement
  | compound_statement
  | '{' ELLIPSIS '}'              { $$ = loop_placeholder (scopes, status, $1, $2, $3); }
  | expression_statement
  | selection_statement
  | iteration_statement
  | jump_statement
  | asm_statement
  ;

labeled_statement:
    label statement               { $$ = N (LABELED_STATEMENT, $1, $2); }
  | label declaration             { $$ = N (LABELED_STATEMENT, $1, LOCAL_DECLARATION ($2)); }
  ;

label:
    any_identifier ':'            { $$ = N (NAMED_LABEL, $1, $2); }
  | CASE constant_expression ':'  { $$ = N (CASE_LABEL, $1, $2, $3); }
  | CASE constant_expression ELLIPSIS constant_expression ':'
                                  { $$ = N (CASE_LABEL, $1, $2, $3, $4, $5); }
  | DEFAULT ':'                   { $$ = N (DEFAULT_LABEL, $1, $2); }
  ;

/* A label may also end a block, as in `default: }`. */
compound_statement:
    '{' open_scope '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $2, $3); }
  | '{' open_scope block_items '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $2, $3, $4); }
  | '{' open_scope label '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $2, $3, $4); }
  | '{' open_scope block_items label '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $2, $3, $4, $5); }
  | '{' open_scope error '}' close_scope
                                  { RECOVERED (); $$ = N (COMPOUND_STATEMENT, $1, $2, $4); }
  | '{' open_scope block_items error '}' close_scope
                                  { RECOVERED (); $$ = N (COMPOUND_STATEMENT, $1, $2, $3, $5); }
  ;

block_items:
    block_item                    { $$ = N (BLOCK_ITEMS, $1); }
  | block_items block_item        { $$ = ast_add ($1, $2); }
  ;

/* Scalars declared in a block get their ids from scalar_new (), and caches
 * a cleanup: LOCAL_DECLARATION (). */
block_item:
    declaration                   { $$ = LOCAL_DECLARATION ($1); }
  | statement
  | label_declaration
  | error ';'                     { RECOVERED (); $$ = NULL; }
  ;

label_declaration:
    LABEL label_list ';'          { $$ = N (LABEL_DECLARATION, $1, $2, $3); }
  ;

label_list:
    any_identifier                { $$ = N (LABEL_LIST, $1); }
  | label_list ',' any_identifier { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

expression_statement:
    ';'                           { $$ = N (EXPRESSION_STATEMENT, $1); }
  | expression ';'                { $$ = N (EXPRESSION_STATEMENT, $1, $2); }
  ;

selection_statement:
    IF '(' expression ')' statement %prec THEN
                                  { $$ = N (IF_STATEMENT, $1, $2, $3, $4, $5); }
  | IF '(' expression ')' statement ELSE statement
                                  { $$ = N (IF_STATEMENT, $1, $2, $3, $4, $5, $6, $7); }
  | SWITCH '(' expression ')' statement
                                  { $$ = N (SWITCH_STATEMENT, $1, $2, $3, $4, $5); }
  ;

iteration_statement:
    WHILE '(' expression ')' statement
                                  { $$ = N (WHILE_STATEMENT, $1, $2, $3, $4, $5); }
  | DO statement WHILE '(' expression ')' ';'
                                  { $$ = N (DO_STATEMENT, $1, $2, $3, $4, $5, $6, $7); }
  | FOR '(' open_scope for_init expression_statement opt_expression ')'
    statement close_scope
                                  { $$ = N (FOR_STATEMENT, $1, $2, $3, $4, $5, $6, $7, $8); }
  /* The loop variable is declared after the collection, which still means the
   * outer name in `foreach (v in v)`, and before the body is read. The scope
   * opens after the `)`, since until `in` this may still be a call. The value
   * of the first action is the declaration of V if it is a cache, looked up
   * before v can hide it. */
  | FOREACH '(' any_identifier foreach_in IDENTIFIER ')' open_scope
                                  { if (foreach_declared (scopes, $1))
                                      status->rejected++;
                                    $$ = (struct ast *)cache_variable (scopes, $5);
                                    bind_name (status, scope_add (scopes, ast_token ($3)->text,
                                                                  SCOPE_VARIABLE, $3), $3); }
    statement
                                  /* The loop's scope, with v in it, is still
                                   * open here: the place to check the names
                                   * that the translation brings in. */
                                  { $$ = foreach_to_for ($1, $2, $3, $4, $5, $6, $9,
                                                         $8 != NULL); }
    close_scope
                                  { $$ = N (FOR_EACH_STATEMENT, $10, $7, $11); }
  /* foreach (), the stencil loop: see translate.h. The body's s[...] are
   * gathered as they reduce, then the loop is expanded with them. */
  | FOREACH_CELLS '(' ')'         { status->stencil = stencil_loop_begin (status->stencil, $1,
                                                                        &status->rejected);
                                    $$ = NULL; }
    statement                     { $$ = stencil_use (scopes, status, $1, NULL, $5);
                                    (void) $2; (void) $3; (void) $4; }
  /* @prolongation (id) (args);: the operator with that id, called */
  | AT_OPERATOR '(' expression ')' '(' opt_arguments ')' ';'
                                  { $$ = operator_dispatch ($1, $2, $3, $4, $5, $6, $7, $8,
                                                            &status->rejected); }
  /* foreach_neighbor (r): the point of a stencil loop moved to each of its
   * neighbours within r, for the statement */
  | FOREACH_NEIGHBOR '(' expression ')'
                                  { stencil_mover_begin (status->stencil, $1, &status->rejected);
                                    $$ = NULL; }
    statement                     { $$ = stencil_mover_end (status->stencil, $1, $2, $3, $4, $6);
                                    (void) $5; }
  /* a loop declared with @foreach, with arguments or not */
  | STENCIL_NAME '(' ')'          { status->stencil = stencil_loop_begin (status->stencil, $1,
                                                                        &status->rejected);
                                    $$ = NULL; }
    statement                     { $$ = stencil_use (scopes, status, $1, NULL, $5);
                                    (void) $2; (void) $3; (void) $4; }
  | STENCIL_NAME '(' argument_list ')'
                                  { status->stencil = stencil_loop_begin (status->stencil, $1,
                                                                        &status->rejected);
                                    $$ = NULL; }
    statement                     { $$ = stencil_use (scopes, status, $1, $3, $6);
                                    (void) $2; (void) $4; (void) $5; }
  /* A use of a block macro. Its statement, the block, is parsed in a scope
   * with the names the loop gives it. */
  | LOOP_NAME '(' ')' open_scope  { loop_open (scopes, status, $1); $$ = NULL; }
    statement                     { $$ = loop_use (scopes, status, $1, NULL, $6); (void) $5; }
    close_scope                   { $$ = $7; (void) $4; (void) $8; }
  | LOOP_NAME '(' argument_list ')' open_scope
                                  { loop_open (scopes, status, $1); $$ = NULL; }
    statement                     { $$ = loop_use (scopes, status, $1, $3, $7); (void) $6; }
    close_scope                   { $$ = $8; (void) $5; (void) $9; }
  ;

/* in is not a token, so that it stays an ordinary name. */
foreach_in:
    any_identifier                { if (strcmp (ast_token ($1)->text, "in") != 0) {
                                      yyerror (lineno, yyscanner, root, scopes, status, "expected 'in'");
                                      YYERROR;
                                    }
                                    $$ = $1; }
  ;

for_init:
    expression_statement
  | declaration                   { $$ = LOCAL_DECLARATION ($1); }
  ;

opt_expression:
    %empty                        { $$ = NULL; }
  | expression
  ;

jump_statement:
    GOTO any_identifier ';'       { $$ = N (GOTO_STATEMENT, $1, $2, $3); }
  | GOTO '*' expression ';'       { $$ = N (GOTO_STATEMENT, $1, $2, $3, $4); }
  | CONTINUE ';'                  { $$ = N (CONTINUE_STATEMENT, $1, $2); }
  | BREAK ';'                     { $$ = N (BREAK_STATEMENT, $1, $2); }
  | RETURN ';'                    { if (status->cur_loop)
                                      loop_error (status, $1, "return in loop '%s', "
                                                  "which would leave the function "
                                                  "that uses it", status->cur_loop->name);
                                    $$ = N (RETURN_STATEMENT, $1, $2); }
  | RETURN expression ';'         { if (status->cur_loop)
                                      loop_error (status, $1, "return in loop '%s', "
                                                  "which would leave the function "
                                                  "that uses it", status->cur_loop->name);
                                    status->rejected += check_cache_copy (scopes, $2, $1);
                                    $$ = N (RETURN_STATEMENT, $1, $2, $3); }
  ;

asm_statement:
    ASM opt_asm_qualifiers '(' string ')' ';'
                                  { $$ = N (ASM_STATEMENT, $1, $2, $3, $4, $5, $6); }
  | ASM opt_asm_qualifiers '(' string ':' opt_asm_operands ')' ';'
                                  { $$ = N (ASM_STATEMENT, $1, $2, $3, $4, $5, $6, $7, $8); }
  | ASM opt_asm_qualifiers '(' string ':' opt_asm_operands ':' opt_asm_operands ')' ';'
                                  { $$ = N (ASM_STATEMENT, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10); }
  | ASM opt_asm_qualifiers '(' string ':' opt_asm_operands ':' opt_asm_operands ':' opt_asm_clobbers ')' ';'
                                  { $$ = N (ASM_STATEMENT, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12); }
  | ASM opt_asm_qualifiers '(' string ':' opt_asm_operands ':' opt_asm_operands ':' opt_asm_clobbers ':' asm_labels ')' ';'
                                  { $$ = N (ASM_STATEMENT, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14); }
  ;

opt_asm_qualifiers:
    %empty                        { $$ = NULL; }
  | asm_qualifiers
  ;

asm_qualifiers:
    asm_qualifier                 { $$ = N (ASM_QUALIFIERS, $1); }
  | asm_qualifiers asm_qualifier  { $$ = ast_add ($1, $2); }
  ;

asm_qualifier:
    VOLATILE
  | INLINE
  | GOTO
  ;

opt_asm_operands:
    %empty                        { $$ = NULL; }
  | asm_operands
  ;

asm_operands:
    asm_operand                   { $$ = N (ASM_OPERANDS, $1); }
  | asm_operands ',' asm_operand  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

asm_operand:
    string '(' expression ')'     { $$ = N (ASM_OPERAND, $1, $2, $3, $4); }
  | '[' any_identifier ']' string '(' expression ')'
                                  { $$ = N (ASM_OPERAND, $1, $2, $3, $4, $5, $6, $7); }
  ;

opt_asm_clobbers:
    %empty                        { $$ = NULL; }
  | asm_clobbers
  ;

asm_clobbers:
    string                        { $$ = N (ASM_CLOBBERS, $1); }
  | asm_clobbers ',' string       { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

asm_labels:
    any_identifier                { $$ = N (ASM_LABELS, $1); }
  | asm_labels ',' any_identifier { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

/* Expressions, from tightest binding to loosest */

primary_expression:
    IDENTIFIER
  | FOREACH %prec KEYWORD_NAME
  | CACHE %prec KEYWORD_NAME
  | CONSTANT
  | string
  | '(' expression ')'            { $$ = N (PAREN_EXPRESSION, $1, $2, $3); }
  | '(' compound_statement ')'    { $$ = N (STATEMENT_EXPRESSION, $1, $2, $3); }
  | generic_selection
  | BUILTIN_VA_ARG '(' assignment_expression ',' type_name ')'
                                  { $$ = N (BUILTIN_VA_ARG, $1, $2, $3, $4, $5, $6); }
  | BUILTIN_OFFSETOF '(' type_name ',' offsetof_member ')'
                                  { $$ = N (BUILTIN_OFFSETOF, $1, $2, $3, $4, $5, $6); }
  | BUILTIN_TYPES_COMPATIBLE '(' type_name ',' type_name ')'
                                  { $$ = N (BUILTIN_TYPES_COMPATIBLE, $1, $2, $3, $4, $5, $6); }
  ;

string:
    STRING                        { $$ = N (STRING_LITERAL, $1); }
  | string STRING                 { $$ = ast_add ($1, $2); }
  ;

offsetof_member:
    any_identifier                { $$ = N (OFFSETOF_MEMBER, $1); }
  | offsetof_member '.' any_identifier
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  | offsetof_member '[' expression ']'
                                  { $$ = ast_add (ast_add (ast_add ($1, $2), $3), $4); }
  ;

generic_selection:
    GENERIC '(' assignment_expression ',' generic_associations ')'
                                  { $$ = N (GENERIC_SELECTION, $1, $2, $3, $4, $5, $6); }
  ;

generic_associations:
    generic_association           { $$ = N (GENERIC_ASSOCIATIONS, $1); }
  | generic_associations ',' generic_association
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

generic_association:
    type_name ':' assignment_expression
                                  { $$ = N (GENERIC_ASSOCIATION, $1, $2, $3); }
  | DEFAULT ':' assignment_expression
                                  { $$ = N (GENERIC_ASSOCIATION, $1, $2, $3); }
  ;

postfix_expression:
    primary_expression
  /* s[...] and s[], when s is a scalar: a cell, in foreach () */
  | postfix_expression '[' expression ']'
                                  { $$ = stencil_access (status->stencil, $1, $2, $3, $4,
                                                         &status->rejected);
                                    if (!$$)
                                      $$ = N (SUBSCRIPT, $1, $2, $3, $4); }
  | postfix_expression '[' ']'    { $$ = stencil_access (status->stencil, $1, $2, NULL, $3,
                                                         &status->rejected); }
  | postfix_expression '(' ')'    { $$ = N (CALL, $1, $2, $3); }
  | FOREACH '(' ')'               { $$ = N (CALL, $1, $2, $3); }
  | FOREACH '(' argument_list ')' { $$ = N (CALL, $1, $2, $3, $4); }
  | CACHE '(' ')'                 { $$ = N (CALL, $1, $2, $3); }
  | CACHE '(' argument_list ')'   { $$ = N (CALL, $1, $2, $3, $4); }
  | postfix_expression '(' argument_list ')'
                                  { check_cache_arguments (scopes, $3);
                                    $$ = N (CALL, $1, $2, $3, $4); }
  | postfix_expression '.' any_identifier
                                  { $$ = N (MEMBER, $1, $2, $3); }
  | postfix_expression ARROW any_identifier
                                  { $$ = N (MEMBER, $1, $2, $3); }
  | postfix_expression INC        { check_scalar_assignment (scopes, $1, $2);
                                    $$ = N (POSTFIX_OP, $1, $2); }
  | postfix_expression DEC        { check_scalar_assignment (scopes, $1, $2);
                                    $$ = N (POSTFIX_OP, $1, $2); }
  | '(' type_name ')' brace_initializer
                                  { $$ = N (COMPOUND_LITERAL, $1, $2, $3, $4); }
  ;

argument_list:
    assignment_expression         { $$ = N (ARGUMENT_LIST, $1); }
  | argument_list ',' assignment_expression
                                  { $$ = ast_add (ast_add ($1, $2), $3); }
  ;

unary_expression:
    postfix_expression
  | INC unary_expression          { check_scalar_assignment (scopes, $2, $1);
                                    $$ = N (PREFIX_OP, $1, $2); }
  | DEC unary_expression          { check_scalar_assignment (scopes, $2, $1);
                                    $$ = N (PREFIX_OP, $1, $2); }
  | unary_operator cast_expression
                                  { $$ = N (UNARY_OP, $1, $2); }
  | LOGICAL_AND any_identifier    { $$ = N (LABEL_ADDRESS, $1, $2); }
  | SIZEOF unary_expression       { $$ = N (SIZEOF_EXPRESSION, $1, $2); }
  | SIZEOF '(' type_name ')'      { $$ = N (SIZEOF_TYPE, $1, $2, $3, $4); }
  | ALIGNOF unary_expression      { $$ = N (ALIGNOF_EXPRESSION, $1, $2); }
  | ALIGNOF '(' type_name ')'     { $$ = N (ALIGNOF_TYPE, $1, $2, $3, $4); }
  | EXTENSION cast_expression     { $$ = N (EXTENSION_EXPRESSION, $1, $2); }
  ;

unary_operator:
    '&'
  | '*'
  | '+'
  | '-'
  | '~'
  | '!'
  | REAL
  | IMAG
  ;

cast_expression:
    unary_expression
  | '(' type_name ')' cast_expression
                                  { $$ = N (CAST, $1, $2, $3, $4); }
  ;

multiplicative_expression:
    cast_expression
  | multiplicative_expression '*' cast_expression
                                  { $$ = N (MULTIPLICATIVE, $1, $2, $3); }
  | multiplicative_expression '/' cast_expression
                                  { $$ = N (MULTIPLICATIVE, $1, $2, $3); }
  | multiplicative_expression '%' cast_expression
                                  { $$ = N (MULTIPLICATIVE, $1, $2, $3); }
  ;

additive_expression:
    multiplicative_expression
  | additive_expression '+' multiplicative_expression
                                  { $$ = N (ADDITIVE, $1, $2, $3); }
  | additive_expression '-' multiplicative_expression
                                  { $$ = N (ADDITIVE, $1, $2, $3); }
  ;

shift_expression:
    additive_expression
  | shift_expression SHL additive_expression
                                  { $$ = N (SHIFT, $1, $2, $3); }
  | shift_expression SHR additive_expression
                                  { $$ = N (SHIFT, $1, $2, $3); }
  ;

relational_expression:
    shift_expression
  | relational_expression '<' shift_expression
                                  { $$ = N (RELATIONAL, $1, $2, $3); }
  | relational_expression '>' shift_expression
                                  { $$ = N (RELATIONAL, $1, $2, $3); }
  | relational_expression LE shift_expression
                                  { $$ = N (RELATIONAL, $1, $2, $3); }
  | relational_expression GE shift_expression
                                  { $$ = N (RELATIONAL, $1, $2, $3); }
  ;

equality_expression:
    relational_expression
  | equality_expression EQ relational_expression
                                  { $$ = N (EQUALITY, $1, $2, $3); }
  | equality_expression NE relational_expression
                                  { $$ = N (EQUALITY, $1, $2, $3); }
  ;

and_expression:
    equality_expression
  | and_expression '&' equality_expression
                                  { $$ = N (BIT_AND, $1, $2, $3); }
  ;

exclusive_or_expression:
    and_expression
  | exclusive_or_expression '^' and_expression
                                  { $$ = N (BIT_XOR, $1, $2, $3); }
  ;

inclusive_or_expression:
    exclusive_or_expression
  | inclusive_or_expression '|' exclusive_or_expression
                                  { $$ = N (BIT_OR, $1, $2, $3); }
  ;

logical_and_expression:
    inclusive_or_expression
  | logical_and_expression LOGICAL_AND inclusive_or_expression
                                  { $$ = N (LOGICAL_AND, $1, $2, $3); }
  ;

logical_or_expression:
    logical_and_expression
  | logical_or_expression LOGICAL_OR logical_and_expression
                                  { $$ = N (LOGICAL_OR, $1, $2, $3); }
  ;

conditional_expression:
    logical_or_expression
  | logical_or_expression '?' expression ':' conditional_expression
                                  { $$ = N (CONDITIONAL, $1, $2, $3, $4, $5); }
  | logical_or_expression '?' ':' conditional_expression
                                  { $$ = N (CONDITIONAL, $1, $2, $3, $4); }
  ;

assignment_expression:
    conditional_expression
  | unary_expression assignment_operator assignment_expression
                                  { check_scalar_assignment (scopes, $1, $2);
                                    /* one error for w = v */
                                    status->rejected += check_cache_copy (scopes, $1, $2)
                                                        || check_cache_copy (scopes, $3, $2);
                                    $$ = N (ASSIGNMENT, $1, $2, $3); }
  ;

assignment_operator:
    '='
  | MUL_ASSIGN
  | DIV_ASSIGN
  | MOD_ASSIGN
  | ADD_ASSIGN
  | SUB_ASSIGN
  | SHL_ASSIGN
  | SHR_ASSIGN
  | AND_ASSIGN
  | XOR_ASSIGN
  | OR_ASSIGN
  ;

expression:
    assignment_expression
  | expression ',' assignment_expression
                                  { $$ = N (COMMA, $1, $2, $3); }
  ;

constant_expression:
    conditional_expression
  ;

%%

/* The first error says what becomes of the file; the rest follow it. */
static void yyerror (int *lineno, yyscan_t yyscanner, struct ast **root,
                     struct scope_table *scopes, struct parse_status *status,
                     const char *msg)
{
  const struct token *t = lex_last (yyscanner);

  (void)lineno;
  (void)root;
  (void)scopes;
  if (++status->errors > MAX_PARSE_ERRORS)
    return;
  diag (t->file, t->line, DIAG_WARNING, "%s%s near '%s'",
        status->errors == 1 ? "not translated, leaving this file to gcc: " : "",
        msg, t->text);
  if (status->errors == MAX_PARSE_ERRORS)
    diag (t->file, t->line, DIAG_NOTE, "%d errors, not looking for more",
          MAX_PARSE_ERRORS);
}

struct ast *parse_file (FILE *in, int *rejected)
{
  struct scope_table *scopes = scope_table_new ();
  yyscan_t scanner = lex_open (in, scopes);
  struct ast *root = NULL;
  struct parse_status status = { 0, 0, NULL, NULL, 0, NULL, NULL };
  int lineno = 1;
  int failed = yyparse (&lineno, scanner, &root, scopes, &status);

  lex_close (scanner);
  *rejected = status.rejected > 0;
  return failed || status.errors || status.rejected ? NULL : root;
}
