#include <stdlib.h>

#include "ast.h"

struct ast *ast_new(enum ast_kind kind)
{
	struct ast *n = calloc(1, sizeof(*n));

	if (!n)
		abort();
	n->kind = kind;
	return n;
}

struct ast *ast_leaf(const struct token *t)
{
	struct ast *n = ast_new(AST_TOKEN);

	n->tok = *t;
	return n;
}

void ast_add(struct ast *parent, struct ast *kid)
{
	if (parent->nkids == parent->cap) {
		parent->cap = parent->cap ? parent->cap * 2 : 8;
		parent->kids = realloc(parent->kids,
				       parent->cap * sizeof(*parent->kids));
		if (!parent->kids)
			abort();
	}
	parent->kids[parent->nkids++] = kid;
}

void ast_emit(const struct ast *n, FILE *out)
{
	if (n->kind == AST_TOKEN) {
		fputs(n->tok.trivia, out);
		fputs(n->tok.text, out);
		return;
	}
	for (size_t i = 0; i < n->nkids; i++)
		ast_emit(n->kids[i], out);
}
