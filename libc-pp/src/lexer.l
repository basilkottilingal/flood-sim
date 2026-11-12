/* Placeholder for the table-based lexer: same contract as lexer.h.
 * Input is `gcc -E` output, so there are no comments or macros left.
 *
 * No file-scope state: the line number comes in as an argument, and the rest
 * lives in struct lexer, reached through yyextra.
 *
 * The argument is not called yylineno, which a reentrant scanner defines as a
 * macro that would also break the prototypes in parser.tab.h. */
%option reentrant noyywrap nounput noinput never-interactive
%option extra-type="struct lexer *"

%{
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "ast.h"
#include "parser.tab.h"
#include "lexer.h"
#include "scope.h"
#include "translate.h"

#define YY_DECL int yylex (YYSTYPE *lvalp, int *lineno, yyscan_t yyscanner)
YY_DECL;

struct lexer {
  yyscan_t scanner;
  struct scope_table *scopes;
  char *trivia;           /* scratch: the only part outside the arena */
  size_t trivia_len, trivia_cap;
  const char *file;
  struct token last;
  int eof_seen;
  int cache_on;           /* <cwrap/cache.h> seen: cache and foreach are keywords */
  const struct var *binding;  /* of the name just read, for its token */
};

/* The header that turns on cache (T) and foreach, as a linemarker names it. */
#define CACHE_HEADER "/cwrap/cache.h"


/* Placeholders: punctuators share these strings instead of the arena. */
static const char *const punctuators[256] = {
  ['['] = "[", [']'] = "]", ['('] = "(", [')'] = ")", ['{'] = "{", ['}'] = "}",
  ['.'] = ".", ['&'] = "&", ['*'] = "*", ['+'] = "+", ['-'] = "-", ['~'] = "~",
  ['!'] = "!", ['/'] = "/", ['%'] = "%", ['<'] = "<", ['>'] = ">", ['^'] = "^",
  ['|'] = "|", ['?'] = "?", [':'] = ":", [';'] = ";", ['='] = "=", [','] = ",",
};

static void trivia_add (struct lexer *lx, const char *s, size_t n)
{
  if (lx->trivia_len + n + 1 > lx->trivia_cap) {
    lx->trivia_cap = (lx->trivia_len + n + 1) * 2;
    lx->trivia = realloc (lx->trivia, lx->trivia_cap);
    if (!lx->trivia)
      abort ();
  }
  memcpy (lx->trivia + lx->trivia_len, s, n);
  lx->trivia_len += n;
  lx->trivia[lx->trivia_len] = '\0';
}

/* Linemarker: # LINE "FILE" FLAGS... Anything else here is a #pragma. */
static void directive (struct lexer *lx, int *lineno, const char *s)
{
  char *end;
  long line;

  s += strspn (s, " \t") + 1;
  line = strtol (s, &end, 10);
  if (end == s) {
    (*lineno)++;
    return;
  }
  *lineno = (int)line;
  s = strchr (end, '"');
  if (s) {
    const char *q = ++s;
    int n;

    while (*q && *q != '"')
      q += (*q == '\\' && q[1]) ? 2 : 1;
    n = (int)(q - s);
    if (strncmp (lx->file, s, n) != 0 || lx->file[n])
      lx->file = allocate_str (s, n);
    if (n >= (int)strlen (CACHE_HEADER)
        && strncmp (s + n - strlen (CACHE_HEADER), CACHE_HEADER,
                    strlen (CACHE_HEADER)) == 0)
      lx->cache_on = 1;
  }
}

/* A name's kind, from the declaration of it visible now, which its token
 * keeps: a typedef name, a loop defined by @loop, or an identifier. */
static int identifier_kind (struct lexer *lx, const char *name)
{
  const struct var *v = scope_find (lx->scopes, name);

  lx->binding = v;
  if (v && v->kind == SCOPE_TYPEDEF)
    return T_TYPE_NAME;
  if (v && v->kind == SCOPE_LOOP)
    return T_LOOP_NAME;
  if (v && v->kind == SCOPE_STENCIL)
    return T_STENCIL_NAME;
  /* N, a stencil loop by axis, if N_x is one (stencil_by_axis ()) */
  if (!v && stencil_enabled ()) {
    char axis[256];
    const struct var *w;

    snprintf (axis, sizeof (axis), "%s_x", name);
    w = scope_find (lx->scopes, axis);
    if (w && w->kind == SCOPE_STENCIL)
      return T_STENCIL_NAME;
  }
  return T_IDENTIFIER;
}

/* cache and foreach, once <cwrap/cache.h> is in: a keyword, but for a typedef.
 * The parser takes the keyword as a name too. */
static int cache_word (struct lexer *lx, const char *name, int keyword)
{
  if (!lx->cache_on || scope_is_type (lx->scopes, name))
    return identifier_kind (lx, name);
  lx->binding = scope_find (lx->scopes, name);
  return keyword;
}

/* foreach (): the stencil loop, with -dimension=, unless foreach is declared:
 * then it is what it was before, as is any other foreach. */
static int stencil_word (struct lexer *lx)
{
  if (!stencil_enabled () || scope_find (lx->scopes, "foreach"))
    return cache_word (lx, "foreach", T_FOREACH);
  lx->binding = NULL;
  return T_FOREACH_CELLS;
}

/* foreach_neighbor, which moves the point of a stencil loop: a word of its
 * own with -dimension=, unless declared. */
static int mover_word (struct lexer *lx, const char *name, int keyword)
{
  if (!stencil_enabled () || scope_find (lx->scopes, name))
    return identifier_kind (lx, name);
  lx->binding = NULL;
  return keyword;
}

/* The leaf for one token, which also takes the trivia gathered before it. */
static struct ast *leaf (struct lexer *lx, int kind, const char *text, int len,
                         int line)
{
  const char *placeholder = len == 1 ? punctuators[(unsigned char)*text] : NULL;

  lx->last.kind = kind;
  lx->last.text = placeholder ? placeholder : allocate_str (text, len);
  lx->last.trivia = lx->trivia_len ? allocate_str (lx->trivia, lx->trivia_len)
                                   : "";
  lx->last.file = lx->file;
  lx->last.line = line;
  lx->last.binding = lx->binding;
  lx->binding = NULL;
  lx->trivia_len = 0;
  return ast_leaf (&lx->last);
}

#define AST(KIND)                                                       \
  do {                                                                  \
    int kind_ = (KIND);                                                 \
    *lvalp = leaf (yyextra, kind_, yytext, yyleng, *lineno);            \
    return kind_;                                                       \
  } while (0)
%}

ID      [A-Za-z_$\x80-\xff][A-Za-z0-9_$\x80-\xff]*
PPNUM   \.?[0-9]([eEpP][+-]|[0-9A-Za-z_.])*
PFX     (u8|u|U|L)

%%

^[ \t]*#[^\n]*\n        { trivia_add (yyextra, yytext, yyleng);
                          directive (yyextra, lineno, yytext); }
\n                      { trivia_add (yyextra, yytext, yyleng); (*lineno)++; }
[ \t\r\f\v]+            { trivia_add (yyextra, yytext, yyleng); }

  /* C11 */
"auto"                  { AST (T_AUTO); }
"break"                 { AST (T_BREAK); }
"case"                  { AST (T_CASE); }
"char"                  { AST (T_CHAR); }
"const"|"__const"|"__const__"                   { AST (T_CONST); }
"continue"              { AST (T_CONTINUE); }
"default"               { AST (T_DEFAULT); }
"do"                    { AST (T_DO); }
"double"                { AST (T_DOUBLE); }
"else"                  { AST (T_ELSE); }
"enum"                  { AST (T_ENUM); }
"extern"                { AST (T_EXTERN); }
"float"                 { AST (T_FLOAT); }
"for"                   { AST (T_FOR); }
"goto"                  { AST (T_GOTO); }
"if"                    { AST (T_IF); }
"inline"|"__inline"|"__inline__"                { AST (T_INLINE); }
"int"                   { AST (T_INT); }
"long"                  { AST (T_LONG); }
"register"              { AST (T_REGISTER); }
"restrict"|"__restrict"|"__restrict__"          { AST (T_RESTRICT); }
"return"                { AST (T_RETURN); }
"short"                 { AST (T_SHORT); }
"signed"|"__signed"|"__signed__"                { AST (T_SIGNED); }
"sizeof"                { AST (T_SIZEOF); }
"static"                { AST (T_STATIC); }
"struct"                { AST (T_STRUCT); }
"switch"                { AST (T_SWITCH); }
"typedef"               { AST (T_TYPEDEF); }
"union"                 { AST (T_UNION); }
"unsigned"              { AST (T_UNSIGNED); }
"void"                  { AST (T_VOID); }
"volatile"|"__volatile"|"__volatile__"          { AST (T_VOLATILE); }
"while"                 { AST (T_WHILE); }
"_Alignas"              { AST (T_ALIGNAS); }
"_Alignof"|"__alignof"|"__alignof__"            { AST (T_ALIGNOF); }
"_Atomic"               { AST (T_ATOMIC); }
"_Bool"                 { AST (T_BOOL); }
"_Complex"|"__complex"|"__complex__"            { AST (T_COMPLEX); }
"_Generic"              { AST (T_GENERIC); }
"_Imaginary"            { AST (T_IMAGINARY); }
"_Noreturn"             { AST (T_NORETURN); }
"_Static_assert"        { AST (T_STATIC_ASSERT); }
"_Thread_local"|"__thread"                      { AST (T_THREAD_LOCAL); }

  /* GNU. Plain asm and typeof are keywords in gnu modes only. */
"__attribute"|"__attribute__"                   { AST (T_ATTRIBUTE); }
"__extension__"         { AST (T_EXTENSION); }
"asm"|"__asm"|"__asm__"                         { AST (T_ASM); }
"typeof"|"__typeof"|"__typeof__"                { AST (T_TYPEOF); }
"__int128"              { AST (T_INT128); }
  /* Floating types that gcc knows as keywords, as in `_Complex _Float32`,
   * where a typedef name could not follow another type. */
"_Float16"|"_Float32"|"_Float64"|"_Float128"|"_Float32x"|"_Float64x"|"_Float128x"  { AST (T_FLOATN); }
"_Decimal32"|"_Decimal64"|"_Decimal128"|"__float80"|"__float128"|"__ibm128"        { AST (T_FLOATN); }
"__real"|"__real__"     { AST (T_REAL); }
"__imag"|"__imag__"     { AST (T_IMAG); }
"__label__"             { AST (T_LABEL); }
"__auto_type"           { AST (T_AUTO_TYPE); }

  /* gcc builtins whose arguments include a type, so the grammar must know them */
"__builtin_va_arg"      { AST (T_BUILTIN_VA_ARG); }
"__builtin_offsetof"    { AST (T_BUILTIN_OFFSETOF); }
"__builtin_types_compatible_p"  { AST (T_BUILTIN_TYPES_COMPATIBLE); }

  /* Ordinary names until <cwrap/cache.h> comes in. */
  /* foreach (): the stencil loop. The parentheses are read again after it. */
"foreach"/[ \t\r\f\v\n]*"("[ \t\r\f\v\n]*")"  { AST (stencil_word (yyextra)); }
"foreach_neighbor"      { AST (mover_word (yyextra, yytext, T_FOREACH_NEIGHBOR)); }
"foreach_child"         { AST (mover_word (yyextra, yytext, T_FOREACH_CHILD)); }
"foreach_dimension"     { AST (mover_word (yyextra, yytext, T_FOREACH_DIMENSION)); }
"foreach"               { AST (cache_word (yyextra, yytext, T_FOREACH)); }
"cache"                 { AST (cache_word (yyextra, yytext, T_CACHE)); }

{ID}                    { AST (identifier_kind (yyextra, yytext)); }
{PPNUM}                 { AST (T_CONSTANT); }
{PFX}?'([^'\\\n]|\\.)+' { AST (T_CONSTANT); }
{PFX}?\"([^"\\\n]|\\.)*\"       { AST (T_STRING); }

"..."   { AST (T_ELLIPSIS); }
"<<="   { AST (T_SHL_ASSIGN); }
">>="   { AST (T_SHR_ASSIGN); }
"+="    { AST (T_ADD_ASSIGN); }
"-="    { AST (T_SUB_ASSIGN); }
"*="    { AST (T_MUL_ASSIGN); }
"/="    { AST (T_DIV_ASSIGN); }
"%="    { AST (T_MOD_ASSIGN); }
"&="    { AST (T_AND_ASSIGN); }
"|="    { AST (T_OR_ASSIGN); }
"^="    { AST (T_XOR_ASSIGN); }
"++"    { AST (T_INC); }
"--"    { AST (T_DEC); }
"->"    { AST (T_ARROW); }
"<<"    { AST (T_SHL); }
">>"    { AST (T_SHR); }
"<="    { AST (T_LE); }
">="    { AST (T_GE); }
"=="    { AST (T_EQ); }
"!="    { AST (T_NE); }
"&&"    { AST (T_LOGICAL_AND); }
"||"    { AST (T_LOGICAL_OR); }
"<:"    { AST ('['); }
":>"    { AST (']'); }
"<%"    { AST ('{'); }
"%>"    { AST ('}'); }

[\[\](){}.&*+\-~!/%<>^|?:;=,]   { AST ((unsigned char)yytext[0]); }

  /* A block macro: see translate.h. C has no other use for @. */
"@loop"/[^A-Za-z0-9_$\x80-\xff]  { AST (T_AT_LOOP); }
  /* A ghost operator, defined, declared or dispatched: see translate.h. */
"@prolongation"/[^A-Za-z0-9_$\x80-\xff]  { AST (T_AT_OPERATOR); }
"@restriction"/[^A-Za-z0-9_$\x80-\xff]   { AST (T_AT_OPERATOR); }
"@boundary"/[^A-Za-z0-9_$\x80-\xff]      { AST (T_AT_OPERATOR); }
  /* A stencil loop: see translate.h. */
"@foreach"/[^A-Za-z0-9_$\x80-\xff]  { AST (T_AT_FOREACH); }

  /* Unknown to the grammar, so it surfaces as a syntax error. */
.       { AST ((unsigned char)yytext[0]); }

  /* T_EOI first, for its trailing trivia: then 0 ends the parse. */
<<EOF>> {
          if (yyextra->eof_seen)
            return 0;
          yyextra->eof_seen = 1;
          *lvalp = leaf (yyextra, T_EOI, "", 0, *lineno);
          return T_EOI;
        }

%%

yyscan_t lex_open (FILE *in, struct scope_table *scopes)
{
  struct lexer *lx = allocate (sizeof (*lx));

  lx->scopes = scopes;
  lx->file = "<input>";
  if (yylex_init_extra (lx, &lx->scanner) != 0)
    abort ();
  yyset_in (in, lx->scanner);
  return lx->scanner;
}

void lex_close (yyscan_t scanner)
{
  struct lexer *lx = yyget_extra (scanner);

  yylex_destroy (scanner);
  free (lx->trivia);
}

const struct token *lex_last (yyscan_t scanner)
{
  return &((struct lexer *)yyget_extra (scanner))->last;
}
