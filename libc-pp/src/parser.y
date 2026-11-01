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
 * symbols to N () as they are: they are skipped.
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

/* Returns the tree, or NULL after a syntax error. The tree is in the arena, so
 * it lasts until free_arena (). */
struct ast *parse_file (FILE *in);
}

%code {
#include <stdlib.h>

#include "allocator.h"
#include "ast.h"
#include "symbol.h"

#define N(KIND, ...)                                                    \
  ast_node (AST_##KIND, (struct ast *[]) { __VA_ARGS__ },               \
            sizeof ((struct ast *[]) { __VA_ARGS__ }) / sizeof (struct ast *))

int yylex (YYSTYPE *lvalp, int *lineno, yyscan_t yyscanner);
static void yyerror (int *lineno, yyscan_t yyscanner, struct ast **root,
                     struct symtab *syms, const char *msg);

static void declare (struct symtab *syms, const struct ast *declarator,
                     int is_type)
{
  const char *name = ast_declared_name (declarator);

  if (name)
    sym_define (syms, name, is_type);
}
}

%define api.pure full
%define api.value.type {struct ast *}
%define api.token.prefix {T_}
%define parse.error verbose
%param {int *lineno} {yyscan_t yyscanner}
%parse-param {struct ast **root} {struct symtab *syms}

%expect 0

%token AUTO BREAK CASE CHAR CONST CONTINUE DEFAULT DO DOUBLE ELSE ENUM EXTERN
%token FLOAT FOR GOTO IF INLINE INT LONG REGISTER RESTRICT RETURN SHORT SIGNED
%token SIZEOF STATIC STRUCT SWITCH TYPEDEF UNION UNSIGNED VOID VOLATILE WHILE
%token ALIGNAS ALIGNOF ATOMIC BOOL COMPLEX GENERIC IMAGINARY NORETURN
%token STATIC_ASSERT THREAD_LOCAL
%token ATTRIBUTE EXTENSION ASM TYPEOF INT128 REAL IMAG LABEL AUTO_TYPE
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
  ;

function_definition:
    decl_specs_typed specs_done declarator enter_function compound_statement
    leave_function
                                  { $$ = N (FUNCTION_DEFINITION, $1, $3, $5); }
  ;

/* Markers: empty rules that run an action at that point of a rule. They reduce
 * before the next token is read, or with one that cannot be a name, so the
 * lexer sees their effect on the very next identifier. */
specs_done:      %empty           { sym_specifiers_done (syms); $$ = NULL; } ;
enter_function:  %empty           { scope_resume_params (syms); $$ = NULL; } ;
leave_function:  %empty           { scope_pop (syms); $$ = NULL; } ;
open_scope:      %empty           { scope_push (syms); $$ = NULL; } ;
close_scope:     %empty           { scope_pop (syms); $$ = NULL; } ;
open_params:     %empty           { scope_push_params (syms); $$ = NULL; } ;
close_params:    %empty           { scope_pop_params (syms); $$ = NULL; } ;

/* Declarations */

declaration:
    decl_specs_typed ';'          { sym_specifiers_reset (syms);
                                    $$ = N (DECLARATION, $1, $2); }
  | decl_specs_untyped ';'        { sym_specifiers_reset (syms);
                                    $$ = N (DECLARATION, $1, $2); }
  | decl_specs_typed specs_done init_declarator_list ';'
                                  { $$ = N (DECLARATION, $1, $3, $4); }
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
    TYPEDEF                       { sym_typedef_seen (syms); $$ = $1; }
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
 * labels, enumerators. */
any_identifier:
    IDENTIFIER
  | TYPE_NAME
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
    declarator                    { declare (syms, $1, sym_in_typedef (syms));
                                    $$ = $1; }
  | declarator declarator_attrs   { declare (syms, $1, sym_in_typedef (syms));
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
    IDENTIFIER                    { sym_declarator_named (syms);
                                    $$ = N (DECL_NAME, $1); }
  | TYPE_NAME                     { sym_declarator_named (syms);
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
                                  { $$ = N (FUNCTION_SUFFIX, $1, $3, $5); }
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
                                  { declare (syms, $2, 0);
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
                                  { $$ = N (FUNCTION_SUFFIX, $1, $3, $5); }
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
    any_identifier                { sym_define (syms, ast_token ($1)->text, 0);
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
  | expression_statement
  | selection_statement
  | iteration_statement
  | jump_statement
  | asm_statement
  ;

labeled_statement:
    label statement               { $$ = N (LABELED_STATEMENT, $1, $2); }
  | label declaration             { $$ = N (LABELED_STATEMENT, $1, $2); }
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
                                  { $$ = N (COMPOUND_STATEMENT, $1, $3); }
  | '{' open_scope block_items '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $3, $4); }
  | '{' open_scope label '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $3, $4); }
  | '{' open_scope block_items label '}' close_scope
                                  { $$ = N (COMPOUND_STATEMENT, $1, $3, $4, $5); }
  ;

block_items:
    block_item                    { $$ = N (BLOCK_ITEMS, $1); }
  | block_items block_item        { $$ = ast_add ($1, $2); }
  ;

block_item:
    declaration
  | statement
  | label_declaration
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
                                  { $$ = N (FOR_STATEMENT, $1, $2, $4, $5, $6, $7, $8); }
  ;

for_init:
    expression_statement
  | declaration
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
  | RETURN ';'                    { $$ = N (RETURN_STATEMENT, $1, $2); }
  | RETURN expression ';'         { $$ = N (RETURN_STATEMENT, $1, $2, $3); }
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
  | postfix_expression '[' expression ']'
                                  { $$ = N (SUBSCRIPT, $1, $2, $3, $4); }
  | postfix_expression '(' ')'    { $$ = N (CALL, $1, $2, $3); }
  | postfix_expression '(' argument_list ')'
                                  { $$ = N (CALL, $1, $2, $3, $4); }
  | postfix_expression '.' any_identifier
                                  { $$ = N (MEMBER, $1, $2, $3); }
  | postfix_expression ARROW any_identifier
                                  { $$ = N (MEMBER, $1, $2, $3); }
  | postfix_expression INC        { $$ = N (POSTFIX_OP, $1, $2); }
  | postfix_expression DEC        { $$ = N (POSTFIX_OP, $1, $2); }
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
  | INC unary_expression          { $$ = N (PREFIX_OP, $1, $2); }
  | DEC unary_expression          { $$ = N (PREFIX_OP, $1, $2); }
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
                                  { $$ = N (ASSIGNMENT, $1, $2, $3); }
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

static void yyerror (int *lineno, yyscan_t yyscanner, struct ast **root,
                     struct symtab *syms, const char *msg)
{
  const struct token *t = lex_last (yyscanner);

  (void)lineno;
  (void)root;
  (void)syms;
  fprintf (stderr, "%s:%d: error: %s near '%s'\n",
           t->file, t->line, msg, t->text);
}

struct ast *parse_file (FILE *in)
{
  struct symtab *syms = symtab_new ();
  yyscan_t scanner = lex_open (in, syms);
  struct ast *root = NULL;
  int lineno = 1;
  int failed = yyparse (&lineno, scanner, &root, syms);

  lex_close (scanner);
  return failed ? NULL : root;
}
