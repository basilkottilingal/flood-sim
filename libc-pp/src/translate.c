#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "ast.h"
#include "parser.tab.h"
#include "translate.h"

/* The translation itself: rewrite the tree here. */
static void transform (struct ast *root)
{
  (void)root;
}

/* gcc gets the input as it came. It either accepts code our grammar does not
 * cover, or reports the error better than we could. */
static int copy_file (FILE *in, FILE *out)
{
  char buf[1 << 16];
  size_t n;

  rewind (in);
  while ((n = fread (buf, 1, sizeof (buf), in)) > 0)
    if (fwrite (buf, 1, n, out) != n)
      return 1;
  return ferror (in);
}

int translate_file (const char *in_path, const char *out_path)
{
  FILE *in, *out;
  struct ast *root;
  int failed;

  in = fopen (in_path, "r");
  if (!in) {
    fprintf (stderr, "cwrap: %s: %s\n", in_path, strerror (errno));
    return 1;
  }
  root = parse_file (in);
  if (root) {
    transform (root);
    if (getenv ("CWRAP_DUMP_AST"))
      ast_dump (root, stderr);
  }

  out = fopen (out_path, "w");
  if (!out) {
    fprintf (stderr, "cwrap: %s: %s\n", out_path, strerror (errno));
    fclose (in);
    free_arena ();
    return 1;
  }
  if (root)
    ast_emit (root, out);
  failed = root ? 0 : copy_file (in, out);
  fclose (in);
  free_arena ();
  failed |= ferror (out);
  if (fclose (out) != 0 || failed) {
    fprintf (stderr, "cwrap: %s: write failed\n", out_path);
    return 1;
  }
  return 0;
}
