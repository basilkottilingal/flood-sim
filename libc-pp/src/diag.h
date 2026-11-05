/* Diagnostics, written and colored as gcc writes them:
 *
 *   file:line: error: message with a 'quoted' name
 *
 * The location is bold, the kind colored, and text in single quotes bold.
 * Colors follow gcc: -fdiagnostics-color=auto (the default: only on a
 * terminal), =always or =never, and GCC_COLORS, whose keys error, warning,
 * note, locus and quote are used, and which turns color off if set but empty. */
#ifndef DIAG_H
#define DIAG_H

enum diag_kind { DIAG_ERROR, DIAG_WARNING, DIAG_NOTE };

/* when is the last -fdiagnostics-color option given, or NULL for none. */
void diag_init (const char *when);

/* A message about file and line. */
void diag (const char *file, int line, enum diag_kind kind, const char *fmt, ...)
  __attribute__ ((format (printf, 4, 5)));

/* A message about the run itself: `cwrap: error: ...`. */
void diag_tool (enum diag_kind kind, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));

#endif
