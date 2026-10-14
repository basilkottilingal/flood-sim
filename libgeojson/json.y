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

#include "ast.h"
#include "lexer.h"

/* forward declaration of functions declared in the lexer */
extern int yylex (YYSTYPE *);
extern int yylineno;
extern int yycolno;

/* error handling */
void yyerror (const char *s);

%}
%token jtrue jfalse jnull jstring jnumber

%%

jjson
    : jvalue                  {$$ = ast_json ($1);}
    ;

jvalue
    : jobject                 {$$ = ast_val ($1); }
    | jarray                  {$$ = ast_val ($1); }
    | jstring                 {$$ = ast_val ($1); }
    | jnumber                 {$$ = ast_val ($1); }
    | jtrue                   {$$ = ast_val ($1); }
    | jfalse                  {$$ = ast_val ($1); }
    | jnull                   {$$ = ast_val ($1); }
    ;

jobject
    : '{' '}'                 {$$ = ast_obj (NULL);}
    | '{' jmembers '}'        {$$ = ast_obj ($2);}
    ;

jmembers
    : jpair                   {$$ = ast_append (NULL, $1);}
    | jmembers ',' jpair      {$$ = ast_append ($1, $3);}
    ;

jpair
    : jstring ':' jvalue      {$$ = ast_pair ($1, $3);}
    ;

jarray
    : '[' ']'                 {$$ = ast_array (NULL);}
    | '[' jelements ']'       {$$ = ast_array ($2);}
    ;

jelements
    : jvalue                  {$$ = ast_append (NULL, $1);}
    | jelements ',' jvalue    {$$ = ast_append ($1, $3);}
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
