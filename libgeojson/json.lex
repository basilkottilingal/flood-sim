%{
/*
.. refer to the combined lexer + grammar rules of json at
.. https://www.json.org/json-en.html 
.. https://www.crockford.com/mckeeman.html
*/

#include "json.tab.h"
#include "ast.h"
#define YY_DECL   int yylex (YYSTYPE * lvalp)
#define AST(TKN)  *lvalp = ast_leaf (yytext); return TKN

%}

%%

true                                            { AST (J_TRUE); }
false                                           { AST (J_FALSE); }
null                                            { AST (J_NULL); }
\"([^"\\\x00-\x1F]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*\"   { AST (J_STRING); }
\[|\]|:|,|\{|\}                                 { AST (yytext [0]); }
-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?  { AST (J_NUMBER); }
[ \t\v\n\f]+                                    { /* skip */ }
.                                               { fprintf (stderr, "\nerror : %s", yytext); }

%%
