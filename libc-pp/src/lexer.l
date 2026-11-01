/* Placeholder for the table-based lexer: same contract as lexer.h.
 * Input is `gcc -E` output, so there are no comments or macros left. */
%option noyywrap nounput noinput never-interactive

%{
#include <stdlib.h>
#include <string.h>

#include "ast.h"
#include "parser.tab.h"
#include "lexer.h"

#define YY_DECL static int scan(void)
YY_DECL;

static char *trivia;
static size_t trivia_len, trivia_cap;
static const char *cur_file = "<input>";
static int cur_line = 1;

static void trivia_add(const char *s, size_t n)
{
	if (trivia_len + n + 1 > trivia_cap) {
		trivia_cap = (trivia_len + n + 1) * 2;
		trivia = realloc(trivia, trivia_cap);
		if (!trivia)
			abort();
	}
	memcpy(trivia + trivia_len, s, n);
	trivia_len += n;
	trivia[trivia_len] = '\0';
}

/* Linemarker: # LINE "FILE" FLAGS... Anything else here is a #pragma. */
static void directive(const char *s)
{
	char *end;
	long line;

	s += strspn(s, " \t") + 1;
	line = strtol(s, &end, 10);
	if (end == s) {
		cur_line++;
		return;
	}
	cur_line = (int)line;
	s = strchr(end, '"');
	if (s) {
		const char *q = ++s;

		while (*q && *q != '"')
			q += (*q == '\\' && q[1]) ? 2 : 1;
		cur_file = strndup(s, (size_t)(q - s));
	}
}
%}

ID	[A-Za-z_$\x80-\xff][A-Za-z0-9_$\x80-\xff]*
PPNUM	\.?[0-9]([eEpP][+-]|[0-9A-Za-z_.])*
PFX	(u8|u|U|L)

%%

^[ \t]*#[^\n]*\n	{ trivia_add(yytext, yyleng); directive(yytext); }
\n			{ trivia_add(yytext, yyleng); cur_line++; }
[ \t\r\f\v]+		{ trivia_add(yytext, yyleng); }

	/* C11 */
"auto"			return T_AUTO;
"break"			return T_BREAK;
"case"			return T_CASE;
"char"			return T_CHAR;
"const"|"__const"|"__const__"			return T_CONST;
"continue"		return T_CONTINUE;
"default"		return T_DEFAULT;
"do"			return T_DO;
"double"		return T_DOUBLE;
"else"			return T_ELSE;
"enum"			return T_ENUM;
"extern"		return T_EXTERN;
"float"			return T_FLOAT;
"for"			return T_FOR;
"goto"			return T_GOTO;
"if"			return T_IF;
"inline"|"__inline"|"__inline__"		return T_INLINE;
"int"			return T_INT;
"long"			return T_LONG;
"register"		return T_REGISTER;
"restrict"|"__restrict"|"__restrict__"		return T_RESTRICT;
"return"		return T_RETURN;
"short"			return T_SHORT;
"signed"|"__signed"|"__signed__"		return T_SIGNED;
"sizeof"		return T_SIZEOF;
"static"		return T_STATIC;
"struct"		return T_STRUCT;
"switch"		return T_SWITCH;
"typedef"		return T_TYPEDEF;
"union"			return T_UNION;
"unsigned"		return T_UNSIGNED;
"void"			return T_VOID;
"volatile"|"__volatile"|"__volatile__"		return T_VOLATILE;
"while"			return T_WHILE;
"_Alignas"		return T_ALIGNAS;
"_Alignof"|"__alignof"|"__alignof__"		return T_ALIGNOF;
"_Atomic"		return T_ATOMIC;
"_Bool"			return T_BOOL;
"_Complex"|"__complex"|"__complex__"		return T_COMPLEX;
"_Generic"		return T_GENERIC;
"_Imaginary"		return T_IMAGINARY;
"_Noreturn"		return T_NORETURN;
"_Static_assert"	return T_STATIC_ASSERT;
"_Thread_local"|"__thread"			return T_THREAD_LOCAL;

	/* GNU. Plain asm and typeof are keywords in gnu modes only. */
"__attribute"|"__attribute__"			return T_ATTRIBUTE;
"__extension__"		return T_EXTENSION;
"asm"|"__asm"|"__asm__"				return T_ASM;
"typeof"|"__typeof"|"__typeof__"		return T_TYPEOF;
"__int128"		return T_INT128;
"__real"|"__real__"	return T_REAL;
"__imag"|"__imag__"	return T_IMAG;
"__label__"		return T_LABEL;
"__auto_type"		return T_AUTO_TYPE;

{ID}			return T_IDENTIFIER;
{PPNUM}			return T_CONSTANT;
{PFX}?'([^'\\\n]|\\.)+'	return T_CONSTANT;
{PFX}?\"([^"\\\n]|\\.)*\"	return T_STRING;

"..."	return T_ELLIPSIS;
"<<="	return T_SHL_ASSIGN;
">>="	return T_SHR_ASSIGN;
"+="	return T_ADD_ASSIGN;
"-="	return T_SUB_ASSIGN;
"*="	return T_MUL_ASSIGN;
"/="	return T_DIV_ASSIGN;
"%="	return T_MOD_ASSIGN;
"&="	return T_AND_ASSIGN;
"|="	return T_OR_ASSIGN;
"^="	return T_XOR_ASSIGN;
"++"	return T_INC;
"--"	return T_DEC;
"->"	return T_ARROW;
"<<"	return T_SHL;
">>"	return T_SHR;
"<="	return T_LE;
">="	return T_GE;
"=="	return T_EQ;
"!="	return T_NE;
"&&"	return T_LOGICAL_AND;
"||"	return T_LOGICAL_OR;
"<:"	return '[';
":>"	return ']';
"<%"	return '{';
"%>"	return '}';

[\[\](){}.&*+\-~!/%<>^|?:;=,]	return yytext[0];

	/* Unknown to the grammar, so it surfaces as a syntax error. */
.	return yytext[0];

%%

void lex_open(FILE *in)
{
	yyin = in;
	yyrestart(in);
	BEGIN(INITIAL);
	cur_file = "<input>";
	cur_line = 1;
	trivia_len = 0;
}

int lex_next(struct token *t)
{
	trivia_len = 0;
	trivia_add("", 0);
	t->kind = scan();
	t->text = strdup(yytext ? (t->kind ? yytext : "") : "");
	t->trivia = strdup(trivia);
	t->file = cur_file;
	t->line = cur_line;
	return t->kind;
}

void lex_close(void)
{
	yylex_destroy();
}
