#ifndef AST_H
#define AST_H

#include <stddef.h>
#include <stdio.h>

#include "lexer.h"

/* One kind per construct in parser.y. To add one: list it here, then build it
 * with N (NAME, ...) in the grammar action. */
#define AST_KINDS(X)                                                          \
  X (TOKEN) X (TRANSLATION_UNIT) X (FUNCTION_DEFINITION)                      \
  X (EXTENSION_DEFINITION) X (ASM_DEFINITION) X (EMPTY_DECLARATION)           \
  X (DECLARATION) X (EXTENSION_DECLARATION) X (STATIC_ASSERT)                 \
  X (DECL_SPECS) X (SPEC_QUALS) X (INIT_DECLARATOR_LIST) X (INIT_DECLARATOR)  \
  X (ATTRIBUTED_DECLARATOR) X (ASM_LABEL) X (DECLARATOR_ATTRS)                \
  X (DECLARATOR) X (PAREN_DECLARATOR) X (SUFFIXED_DECLARATOR) X (DECL_NAME)   \
  X (POINTER) X (QUALIFIER_LIST) X (ARRAY_SUFFIX) X (FUNCTION_SUFFIX)         \
  X (PARAMETER_TYPE_LIST) X (PARAMETER_LIST) X (PARAMETER)                    \
  X (ABSTRACT_DECLARATOR) X (TYPE_NAME)                                       \
  X (STRUCT_SPECIFIER) X (STRUCT_DECLARATION_LIST) X (STRUCT_DECLARATION)     \
  X (STRUCT_DECLARATOR_LIST) X (STRUCT_DECLARATOR)                            \
  X (ENUM_SPECIFIER) X (ENUMERATOR_LIST) X (ENUMERATOR)                       \
  X (TYPEOF_SPECIFIER) X (ATOMIC_SPECIFIER) X (ALIGNMENT_SPECIFIER)           \
  X (CACHE_SPECIFIER)                                                         \
  X (ATTRIBUTE) X (ATTRIBUTE_LIST) X (ATTRIBUTE_ITEM) X (ATTRIBUTES)          \
  X (INITIALIZER_LIST_EXPR) X (INITIALIZER_LIST) X (INITIALIZER_ITEM)         \
  X (DESIGNATION) X (DESIGNATOR) X (OLD_DESIGNATOR)                           \
  X (COMPOUND_STATEMENT) X (BLOCK_ITEMS) X (LABEL_DECLARATION)                \
  X (LABEL_LIST) X (LABELED_STATEMENT) X (NAMED_LABEL) X (CASE_LABEL)         \
  X (DEFAULT_LABEL) X (EXPRESSION_STATEMENT) X (IF_STATEMENT)                 \
  X (SWITCH_STATEMENT) X (WHILE_STATEMENT) X (DO_STATEMENT) X (FOR_STATEMENT) \
  X (FOR_EACH_STATEMENT)                                                      \
  X (GOTO_STATEMENT) X (CONTINUE_STATEMENT) X (BREAK_STATEMENT)               \
  X (RETURN_STATEMENT) X (ASM_STATEMENT) X (ASM_QUALIFIERS)                   \
  X (ASM_OPERANDS) X (ASM_OPERAND) X (ASM_CLOBBERS) X (ASM_LABELS)            \
  X (STRING_LITERAL) X (PAREN_EXPRESSION) X (STATEMENT_EXPRESSION)            \
  X (GENERIC_SELECTION) X (GENERIC_ASSOCIATIONS) X (GENERIC_ASSOCIATION)      \
  X (BUILTIN_VA_ARG) X (BUILTIN_OFFSETOF) X (OFFSETOF_MEMBER)                 \
  X (BUILTIN_TYPES_COMPATIBLE) X (SUBSCRIPT) X (CALL) X (ARGUMENT_LIST)       \
  X (MEMBER) X (POSTFIX_OP) X (COMPOUND_LITERAL)                              \
  X (PREFIX_OP) X (UNARY_OP) X (LABEL_ADDRESS) X (SIZEOF_EXPRESSION)          \
  X (SIZEOF_TYPE) X (ALIGNOF_EXPRESSION) X (ALIGNOF_TYPE)                     \
  X (EXTENSION_EXPRESSION) X (CAST) X (MULTIPLICATIVE) X (ADDITIVE) X (SHIFT) \
  X (RELATIONAL) X (EQUALITY) X (BIT_AND) X (BIT_XOR) X (BIT_OR)              \
  X (LOGICAL_AND) X (LOGICAL_OR) X (CONDITIONAL) X (ASSIGNMENT) X (COMMA)    \
  X (LOOP_DEFINITION) X (LOOP_SOURCE) X (LOOP_INJECTS) X (LOOP_INJECT_LIST)  \
  X (LOOP_BODY) X (STENCIL_ACCESS) X (STENCIL_LOOP)

enum ast_kind {
#define X(name) AST_##name,
  AST_KINDS (X)
#undef X
};

/* Children are chained: child is the first, last the most recent, and next
 * links siblings. */
struct ast {
  enum ast_kind kind;
  struct ast *parent;
  struct ast *child, *last;
  struct ast *next;
};

struct ast_leaf {         /* kind == AST_TOKEN */
  struct ast a;
  struct token tok;
};

/* All nodes come from the arena: free_arena () releases a whole tree. */
struct ast *ast_new (enum ast_kind kind);
struct ast *ast_leaf (const struct token *t);

/* New node of the given kind, over kids[0..n). NULL kids are skipped, which
 * lets optional parts of a rule be passed as they are. */
struct ast *ast_node (enum ast_kind kind, struct ast **kids, size_t n);

/* A copy of the tree under n, tokens included, with no parent. */
struct ast *ast_copy (const struct ast *n);

/* Appends kid to parent, unless kid is NULL. Returns parent. */
struct ast *ast_add (struct ast *parent, struct ast *kid);

/* Puts n where old is under its parent, and returns old, now detached. */
struct ast *ast_replace (struct ast *old, struct ast *n);

/* Puts n just before pos, under pos's parent. Returns n. */
struct ast *ast_insert_before (struct ast *pos, struct ast *n);

const struct token *ast_token (const struct ast *n);

/* The first token under n, or n's own if it is a leaf. */
const struct token *ast_first_token (const struct ast *n);

/* Text of the name a declarator declares, or NULL if it has none. */
const char *ast_declared_name (const struct ast *declarator);

/* Writes the tokens under n, each preceded by its trivia. */
void ast_emit (const struct ast *n, FILE *out);

/* Writes the tree under n as an indented outline, for debugging. */
void ast_dump (const struct ast *n, FILE *out);

#endif
