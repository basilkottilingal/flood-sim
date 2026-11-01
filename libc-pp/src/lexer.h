/* The only surface the rest of cwrap sees of the lexer. To switch to another
 * lexer generator, implement these functions and yylex (), and drop lexer.l.
 *
 * Token kinds are bison's: T_* from parser.tab.h, a plain character for
 * single-character punctuators, and 0 after T_EOI.
 *
 * yylex () creates one leaf per token and hands it to the parser through its
 * semantic value. The line number is an argument of both yyparse () and
 * yylex (). The rest of the lexer's state belongs to the scanner handle. */
#ifndef LEXER_H
#define LEXER_H

#include <stdio.h>

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void *yyscan_t;
#endif

struct symtab;

struct token {
  int kind;
  const char *text;
  const char *trivia;     /* whitespace, linemarkers and #pragma lines
                           * before this token, kept verbatim */
  const char *file;       /* original location, from linemarkers */
  int line;
};

yyscan_t lex_open (FILE *in, struct symtab *syms);
void lex_close (yyscan_t scanner);

/* The token most recently read, for error messages. */
const struct token *lex_last (yyscan_t scanner);

#endif
