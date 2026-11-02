#ifndef TRANSLATE_H
#define TRANSLATE_H

/* Translates preprocessed C from in_path to out_path. Input that does not parse
 * is copied as it is, for gcc to judge. Returns non-zero only on an I/O error. */
int translate_file (const char *in_path, const char *out_path);

#endif
