%{
/*
.. refer to the combined lexer + grammar rules of json at
.. https://www.json.org/json-en.html 
.. https://www.crockford.com/mckeeman.html
*/

#include "json.tab.h"
#include "ast.h"
#define YY_DECL   int yylex (YYSTYPE * lvalp)
#define AST(TKN)  *lvalp = ast_leaf (yytext, yyleng, ast_token (TKN)); return TKN

%}

%%

true                                            { AST (jtrue); }
false                                           { AST (jfalse); }
null                                            { AST (jnull); }
\"([^"\\\x00-\x1F]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*\"   { AST (jstring); }
\[|\]|:|,|\{|\}                                 { return (yytext [0]); }
-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?  { AST (jnumber); }
[ \t\v\n\f]+                                    { /* skip */ }
.                                               { fprintf (stderr, "\nerror : %s", yytext); }

%%
