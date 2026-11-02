#ifndef TRANSLATE_H
#define TRANSLATE_H

enum {
  TRANSLATE_OK,
  TRANSLATE_FAILED,       /* I/O error */
  TRANSLATE_SKIPPED,      /* did not parse, so copied as it is for gcc to judge */
};

/* Translates preprocessed C from in_path to out_path. Returns a TRANSLATE_ value. */
int translate_file (const char *in_path, const char *out_path);

#endif
