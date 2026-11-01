#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "ast.h"
#include "parser.tab.h"
#include "translate.h"

/* The translation itself: rewrite the tree here. */
static void transform(struct ast *root)
{
	(void)root;
}

int translate_file(const char *in_path, const char *out_path)
{
	FILE *in, *out;
	struct ast *root;
	int failed;

	in = fopen(in_path, "r");
	if (!in) {
		fprintf(stderr, "cwrap: %s: %s\n", in_path, strerror(errno));
		return 1;
	}
	root = parse_file(in);
	fclose(in);
	if (!root)
		return 1;

	transform(root);

	out = fopen(out_path, "w");
	if (!out) {
		fprintf(stderr, "cwrap: %s: %s\n", out_path, strerror(errno));
		return 1;
	}
	ast_emit(root, out);
	failed = ferror(out);
	if (fclose(out) != 0 || failed) {
		fprintf(stderr, "cwrap: %s: write failed\n", out_path);
		return 1;
	}
	return 0;
}
