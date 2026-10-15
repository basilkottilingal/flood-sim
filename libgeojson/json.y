%code requires {
  /* these lines will be visible in json.tab.h */
  #include "ast.h"

  /* fixme : lexer (github.com/basilkottilingal/lexer) cannot
  .. split the lexer file into a header and .c. So we cannot
  .. use multiple instances of lexer.h (produced using json.l)
  .. Thus we limit the use of lexer.h only in this parser. API to
  .. the lexer is thus called via parser (json.tab.h)
  */

  /* call before yyparse () */
  void yyparse_open_file (const char * file);

  /* can be called to clean the lexer buffer and fclose () the source */
  void yyparse_close_file ();
}

%define api.pure full
%define api.value.type {Ast *}
%parse-param {Ast ** json}
%parse-param {Allocator * alloc}
%lex-param   {Allocator * alloc}

%{

/*
.. This parser is generated from json.y using [the table based parser 
.. generator] bison. Note that -d option of bison splits the parser into
.. two [json.tab.h and json.tab.c]
.. The header will contain YYSTYPE, yyparse forward declarations.
.. You can update the parser using :
.. $ bison -d json.y
.. however it's not recommended to put this in the Makefile. Run it only
.. when you want to update the parser.
..
.. Refer to the combined lexer + grammar rules of json at
.. https://www.json.org/json-en.html 
.. https://www.crockford.com/mckeeman.html
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "ast.h"
#include "lexer.h"

/* forward declaration of functions declared in the lexer */
extern int yylex (Ast **, Allocator *);
extern int yylineno;
extern int yycolno;

/* error handling */
void yyerror (Ast **, Allocator *, const char *s);

%}
%token jtrue jfalse jnull jstring jnumber

%%

jjson
    : jvalue                  {*json = $$ = ast_json (alloc, $1);}
    ;

jvalue
    : jobject                 {$$ = ast_val (alloc, $1); }
    | jarray                  {$$ = ast_val (alloc, $1); }
    | jstring                 {$$ = ast_val (alloc, $1); }
    | jnumber                 {$$ = ast_val (alloc, $1); }
    | jtrue                   {$$ = ast_val (alloc, $1); }
    | jfalse                  {$$ = ast_val (alloc, $1); }
    | jnull                   {$$ = ast_val (alloc, $1); }
    ;

jobject
    : '{' '}'                 {$$ = ast_obj (alloc, NULL);}
    | '{' jmembers '}'        {$$ = ast_obj (alloc, $2);}
    ;

jmembers
    : jpair                   {$$ = ast_append (alloc, NULL, $1);}
    | jmembers ',' jpair      {$$ = ast_append (alloc, $1, $3);}
    ;

jpair
    : jstring ':' jvalue      {$$ = ast_pair (alloc, $1, $3);}
    ;

jarray
    : '[' ']'                 {$$ = ast_array (alloc, NULL);}
    | '[' jelements ']'       {$$ = ast_array (alloc, $2);}
    ;

jelements
    : jvalue                  {$$ = ast_append (alloc, NULL, $1);}
    | jelements ',' jvalue    {$$ = ast_append (alloc, $1, $3);}
    ;

%%

#include "json.tab.h"

/* Error handling routine called automatically by Yacc on syntax error */
void yyerror (Ast ** a, Allocator * alloc, const char *s)
{
  (void) a;
  (void) alloc;
  fprintf(stderr, "JSON Syntax Error: %s at line %d\n", s, yylineno);
}

void yyparse_open_file (const char * f)
{
  /* in absence of calling, lexer will assume stdin as the source */
  lxr_source (f);
}

void yyparse_close_file ()
{
  /*
  .. warning : will clean all the token buffers. recommend
  .. a copy of tokens (may use strdup (yytext)) while creating AST
  */
  lxr_clean ();
}

#ifdef JSON_PARSER_TESTING
int main (void)
{
  printf("Parsing JSON input...\n");
  Ast * json = NULL;
  if (yyparse(&json, allocator ()) == 0)
  {
    printf ("Result: Valid JSON document!\n");
    assert (json);
    return 0;
  } else
  {
    printf("Result: Invalid JSON.\n");
    return 1;
  }
  lxr_clean ();
  if (json)
    ast_free (json);
}
#endif
