#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "diag.h"

/* SGR codes: gcc's defaults, which GCC_COLORS may change. */
enum { C_ERROR, C_WARNING, C_NOTE, C_LOCUS, C_QUOTE, NCOLORS };

static const char *const color_keys[NCOLORS] = {
  "error", "warning", "note", "locus", "quote",
};
static char colors[NCOLORS][32] = {
  "01;31", "01;35", "01;36", "01", "01",
};
static const char *const kind_names[] = { "error", "warning", "note" };

static int ready, use_color;

/* key=sgr:key=sgr:... into colors. A key not set keeps its default, as an
 * empty value turns that color off. */
static void read_gcc_colors (const char *s)
{
  while (*s) {
    size_t key = strcspn (s, "="), end = strcspn (s, ":");

    if (key < end)
      for (int i = 0; i < NCOLORS; i++)
        if (strlen (color_keys[i]) == key && strncmp (s, color_keys[i], key) == 0
            && end - key - 1 < sizeof (colors[i])) {
          memcpy (colors[i], s + key + 1, end - key - 1);
          colors[i][end - key - 1] = '\0';
        }
    s += end + (s[end] == ':');
  }
}

void diag_init (const char *when)
{
  const char *gcc_colors = getenv ("GCC_COLORS");
  const char *term = getenv ("TERM");

  ready = 1;
  if (when && strcmp (when, "-fno-diagnostics-color") == 0)
    use_color = 0;
  else if (when && (strcmp (when, "-fdiagnostics-color") == 0
                    || strcmp (when, "-fdiagnostics-color=always") == 0))
    use_color = 1;
  else if (when && strcmp (when, "-fdiagnostics-color=never") == 0)
    use_color = 0;
  else
    use_color = isatty (STDERR_FILENO) && term && strcmp (term, "dumb") != 0;
  if (gcc_colors && !*gcc_colors)
    use_color = 0;
  else if (gcc_colors)
    read_gcc_colors (gcc_colors);
}

static void start (int c)
{
  if (use_color && colors[c][0])
    fprintf (stderr, "\033[%sm\033[K", colors[c]);
}

static void stop (int c)
{
  if (use_color && colors[c][0])
    fputs ("\033[m\033[K", stderr);
}

/* The message, with the text between pairs of single quotes in bold. An odd
 * quote, as in an apostrophe, leaves the rest as it is. */
static void message (const char *fmt, va_list ap)
{
  char buf[1024], *s = buf;
  int n = vsnprintf (buf, sizeof (buf), fmt, ap);
  int quotes = 0;

  if (n < 0)
    return;
  for (const char *p = buf; *p; p++)
    quotes += *p == '\'';
  for (char *q; quotes >= 2 && (q = strchr (s, '\'')); quotes -= 2) {
    char *close = strchr (q + 1, '\'');

    fwrite (s, 1, (size_t)(q - s), stderr);
    fputc ('\'', stderr);
    start (C_QUOTE);
    fwrite (q + 1, 1, (size_t)(close - q - 1), stderr);
    stop (C_QUOTE);
    fputc ('\'', stderr);
    s = close + 1;
  }
  fputs (s, stderr);
  if ((size_t)n >= sizeof (buf))
    fputs ("...", stderr);
  fputc ('\n', stderr);
}

static void kind_label (enum diag_kind kind)
{
  start (kind);
  fprintf (stderr, "%s:", kind_names[kind]);
  stop (kind);
  fputc (' ', stderr);
}

void diag (const char *file, int line, enum diag_kind kind, const char *fmt, ...)
{
  va_list ap;

  if (!ready)
    diag_init (NULL);
  start (C_LOCUS);
  fprintf (stderr, "%s:%d:", file, line);
  stop (C_LOCUS);
  fputc (' ', stderr);
  kind_label (kind);
  va_start (ap, fmt);
  message (fmt, ap);
  va_end (ap);
}

void diag_tool (enum diag_kind kind, const char *fmt, ...)
{
  va_list ap;

  if (!ready)
    diag_init (NULL);
  start (C_LOCUS);
  fputs ("cwrap:", stderr);
  stop (C_LOCUS);
  fputc (' ', stderr);
  kind_label (kind);
  va_start (ap, fmt);
  message (fmt, ap);
  va_end (ap);
}
