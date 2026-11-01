/* Stub grammar: accepts any token sequence as a flat list, so the pipeline runs
 * end to end today. Replace `translation_unit` and `any` with the C grammar;
 * the token set, AST plumbing and parse_file() can stay. */
%code requires {
#include <stdio.h>
struct ast;
struct ast *parse_file(FILE *in);
}

%code {
#include <stdlib.h>

#include "ast.h"
#include "lexer.h"

static struct token last;
static struct ast *eof_leaf;

static int yylex(void);
static void yyerror(struct ast **root, const char *msg);
}

%define api.value.type {struct ast *}
%define api.token.prefix {T_}
%define parse.error verbose
%parse-param {struct ast **root}

%token AUTO BREAK CASE CHAR CONST CONTINUE DEFAULT DO DOUBLE ELSE ENUM EXTERN
%token FLOAT FOR GOTO IF INLINE INT LONG REGISTER RESTRICT RETURN SHORT SIGNED
%token SIZEOF STATIC STRUCT SWITCH TYPEDEF UNION UNSIGNED VOID VOLATILE WHILE
%token ALIGNAS ALIGNOF ATOMIC BOOL COMPLEX GENERIC IMAGINARY NORETURN
%token STATIC_ASSERT THREAD_LOCAL
%token ATTRIBUTE EXTENSION ASM TYPEOF INT128 REAL IMAG LABEL AUTO_TYPE
%token IDENTIFIER CONSTANT STRING
%token ELLIPSIS SHL_ASSIGN SHR_ASSIGN ADD_ASSIGN SUB_ASSIGN MUL_ASSIGN
%token DIV_ASSIGN MOD_ASSIGN AND_ASSIGN OR_ASSIGN XOR_ASSIGN
%token INC DEC ARROW SHL SHR LE GE EQ NE LOGICAL_AND LOGICAL_OR

%%

start:
	translation_unit	{ *root = $1; }
	;

translation_unit:
	%empty			{ $$ = ast_new(AST_TRANSLATION_UNIT); }
	| translation_unit any	{ ast_add($1, $2); $$ = $1; }
	;

any:
	  AUTO | BREAK | CASE | CHAR | CONST | CONTINUE | DEFAULT | DO | DOUBLE
	| ELSE | ENUM | EXTERN | FLOAT | FOR | GOTO | IF | INLINE | INT | LONG
	| REGISTER | RESTRICT | RETURN | SHORT | SIGNED | SIZEOF | STATIC
	| STRUCT | SWITCH | TYPEDEF | UNION | UNSIGNED | VOID | VOLATILE | WHILE
	| ALIGNAS | ALIGNOF | ATOMIC | BOOL | COMPLEX | GENERIC | IMAGINARY
	| NORETURN | STATIC_ASSERT | THREAD_LOCAL
	| ATTRIBUTE | EXTENSION | ASM | TYPEOF | INT128 | REAL | IMAG | LABEL
	| AUTO_TYPE
	| IDENTIFIER | CONSTANT | STRING
	| ELLIPSIS | SHL_ASSIGN | SHR_ASSIGN | ADD_ASSIGN | SUB_ASSIGN
	| MUL_ASSIGN | DIV_ASSIGN | MOD_ASSIGN | AND_ASSIGN | OR_ASSIGN
	| XOR_ASSIGN | INC | DEC | ARROW | SHL | SHR | LE | GE | EQ | NE
	| LOGICAL_AND | LOGICAL_OR
	| '[' | ']' | '(' | ')' | '{' | '}' | '.' | '&' | '*' | '+' | '-'
	| '~' | '!' | '/' | '%' | '<' | '>' | '^' | '|' | '?' | ':' | ';'
	| '=' | ','
	;

%%

/* Every token becomes a leaf, so grammar actions can attach them to the tree
 * and emission stays lossless.
 *
 * Typedef names belong here: when the real grammar arrives, return a distinct
 * TYPE_NAME for an IDENTIFIER the symbol table knows as a typedef. Lexing on
 * demand (not ahead of time) is what lets declarations seen so far decide. */
static int yylex(void)
{
	lex_next(&last);
	yylval = ast_leaf(&last);
	if (last.kind == 0)
		eof_leaf = yylval;
	return last.kind;
}

static void yyerror(struct ast **root, const char *msg)
{
	(void)root;
	fprintf(stderr, "%s:%d: error: %s near '%s'\n",
		last.file, last.line, msg, last.text);
}

struct ast *parse_file(FILE *in)
{
	struct ast *result = NULL;

	lex_open(in);
	eof_leaf = NULL;
	if (yyparse(&result) != 0)
		result = NULL;
	else
		ast_add(result, eof_leaf);	/* carries trailing trivia */
	lex_close();
	return result;
}
