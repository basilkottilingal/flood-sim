#include "json.tab.h"
#include "ast.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main (int argc, const char ** argv)
{
  /*
  .. $ make json.parse.s && ./json.parse.s ../libgeojson/test.json
  */

  fprintf(stderr, "Parsing JSON input...\n");
  Ast * json = NULL;
  if (argc > 1)
    yyparse_open_file (argv [1]);
  if (yyparse (&json, allocator ()) == 0)
  {
    fprintf (stderr, "Result: Valid JSON document!\n");
    assert (json);
  } else
  {
    fprintf(stderr, "Result: Invalid JSON.\n");
    exit (-1);
  }
  ast_print (json);
  yyparse_close_file ();
  if (json)
    ast_free (json);
}
