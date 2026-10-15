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

  printf("Parsing JSON input...\n");
  Ast * json = NULL;
  if (argc > 1)
    yyparse_open_file (argv [1]);
  if (yyparse (&json, allocator ()) == 0)
  {
    printf ("Result: Valid JSON document!\n");
    assert (json);
    return 0;
  } else
  {
    printf("Result: Invalid JSON.\n");
    return 1;
  }
  yyparse_close_file ();
  if (json)
    ast_free (json);
}
