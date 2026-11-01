#ifndef AST_H
#define AST_H

#include <stddef.h>
#include <stdio.h>

#include "lexer.h"

enum ast_kind {
	AST_TOKEN,
	AST_TRANSLATION_UNIT,
	/* grammar node kinds go here */
};

struct ast {
	enum ast_kind kind;
	struct token tok;	/* AST_TOKEN only */
	struct ast **kids;
	size_t nkids, cap;
};

struct ast *ast_new(enum ast_kind kind);
struct ast *ast_leaf(const struct token *t);
void ast_add(struct ast *parent, struct ast *kid);

/* Writes the tokens under n, each preceded by its trivia. */
void ast_emit(const struct ast *n, FILE *out);

#endif
