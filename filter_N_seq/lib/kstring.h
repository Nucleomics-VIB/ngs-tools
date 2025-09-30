/* Copied from klib/kstring.h for minimal build */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

/* Additional kstring functions and macros can be added here if needed by kseq.h */
