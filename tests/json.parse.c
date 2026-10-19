#include "json.tab.h"
#include "ast.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main (int argc, const char ** argv)
{
  /*
  .. $ make json.parse.s && ./json.parse.s ../libgeojson/test.json
  .. $ make json.parse.s && ./json.parse.s < ../libgeojson/test.json
  */

  fprintf(stderr, "Parsing JSON input...\n");
  Ast * json = NULL;
  int y =
    argc > 1 ? yyparse_file (&json, argv [1]) :
               yyparse (&json, allocator ());
  if (!y)
  {
    fprintf (stderr, "Result: Valid JSON document!\n");
    assert (json);
  } else
  {
    fprintf(stderr, "Result: Invalid JSON.\n");
    exit (-1);
  }
  ast_print (json);
  ast_free (json);
}
