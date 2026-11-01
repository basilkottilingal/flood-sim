/* The only surface the rest of cwrap sees of the lexer. To switch to another
 * lexer generator, implement these three functions and drop lexer.l.
 *
 * Token kinds are bison's: T_* from parser.tab.h, a plain character for
 * single-character punctuators, and 0 at end of input. */
#ifndef LEXER_H
#define LEXER_H

#include <stdio.h>

struct token {
	int kind;
	const char *text;
	const char *trivia;	/* whitespace, linemarkers and #pragma lines
				 * before this token, kept verbatim */
	const char *file;	/* original location, from linemarkers */
	int line;
};

void lex_open(FILE *in);
int lex_next(struct token *t);
void lex_close(void);

#endif
