%code requires {
/* these line will be visible in json.tab.h */
#include "ast.h"
}
%define api.pure full
%define api.value.type {Ast *}

%{

/*
.. This parser is generated from json.y using [the table based parser 
.. generator] bison. Note that -d option of bison splits the parser into
.. two [json.tab.h and json.tab.c]
.. The header will contain YYSTYPE, yyparse forward declarations.
.. $ bison -d json.y
..
.. Refer to the combined lexer + grammar rules of json at
.. https://www.json.org/json-en.html 
.. https://www.crockford.com/mckeeman.html
*/

#include <stdio.h>
#include <stdlib.h>

#include "lexer.h"
#include "ast.h"

/* forward declaration of functions declared in the lexer */
extern int yylex (YYSTYPE *);
extern int yylineno;
extern int yycolno;

/* error handling */
void yyerror (const char *s);

%}
%token J_TRUE J_FALSE J_NULL J_STRING J_NUMBER

%%

JSON
    : J_VALUE                  {$$ = ast1 ($1);}
    ;

J_VALUE
    : J_OBJECT                 {$$ = ast1 ($1);}
    | J_ARRAY                  {$$ = ast1 ($1);}
    | J_STRING                 {$$ = ast1 ($1);}
    | J_NUMBER                 {$$ = ast1 ($1);}
    | J_TRUE                   {$$ = ast1 ($1);}
    | J_FALSE                  {$$ = ast1 ($1);}
    | J_NULL                   {$$ = ast1 ($1);}
    ;

J_OBJECT
    : '{' '}'                  {$$ = ast2 ($1, $2);}
    | '{' J_MEMBERS '}'        {$$ = ast3 ($1, $2, $3);}
    ;

J_MEMBERS
    : J_PAIR                   {$$ = ast1 ($1);}
    | J_MEMBERS ',' J_PAIR     {$$ = ast3 ($1, $2, $3);}
    ;

J_PAIR
    : J_STRING ':' J_VALUE     {$$ = ast3 ($1, $2, $3);}
    ;

J_ARRAY
    : '[' ']'                  {$$ = ast2 ($1, $2);}
    | '[' J_ELEMENTS ']'       {$$ = ast3 ($1, $2, $3);}
    ;

J_ELEMENTS
    : J_VALUE                  {$$ = ast1 ($1);}
    | J_ELEMENTS ',' J_VALUE   {$$ = ast3 ($1, $2, $3);}
    ;

%%

#include "json.tab.h"

/* Error handling routine called automatically by Yacc on syntax error */
void yyerror (const char *s)
{
  fprintf(stderr, "JSON Syntax Error: %s at line %d\n", s, yylineno);
}

#ifdef JSON_PARSER_TESTING
int main (void)
{
  printf("Parsing JSON input...\n");
  if (yyparse() == 0)
  {
    printf("Result: Valid JSON document!\n");
    return 0;
  } else
  {
    printf("Result: Invalid JSON.\n");
    return 1;
  }
}
#endif
