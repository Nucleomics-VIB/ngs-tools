/* Copied from klib/kseq.h for minimal build */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

#define KS_SEP_SPACE 0
#define KS_SEP_TAB   1
#define KS_SEP_LINE  2
#define KS_SEP_MAX   2

#define __KS_TYPE(type_t) \
    typedef struct __kstream_t { \
        unsigned char *buf; \
        int begin, end, is_eof; \
        type_t f; \
    } kstream_t;

#define ks_err(ks) ((ks)->end == -1)
#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#define __KS_BASIC(type_t, __bufsize) \
    static inline kstream_t *ks_init(type_t f) { \
        kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t)); \
        ks->f = f; \
        ks->buf = (unsigned char*)malloc(__bufsize); \
        return ks; \
    } \
    static inline void ks_destroy(kstream_t *ks) { \
        if (ks) { \
            free(ks->buf); \
            free(ks); \
        } \
    }

#define __KS_GETC(__read, __bufsize) \
    static inline int ks_getc(kstream_t *ks) { \
        if (ks_err(ks)) return -3; \
        if (ks->is_eof && ks->begin >= ks->end) return -1; \
        if (ks->begin >= ks->end) { \
            ks->begin = 0; \
            ks->end = __read(ks->f, ks->buf, __bufsize); \
            if (ks->end == 0) { ks->is_eof = 1; return -1;} \
            if (ks->end == -1) { ks->is_eof = 1; return -3;} \
        } \
        return (int)ks->buf[ks->begin++]; \
    }

#define __KS_GETUNTIL(__read, __bufsize) \
    static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append) { \
        int gotany = 0; \
        if (dret) *dret = 0; \
        str->l = append? str->l : 0; \
        for (;;) { \
            int i; \
            if (ks_err(ks)) return -3; \
            if (ks->begin >= ks->end) { \
                if (!ks->is_eof) { \
                    ks->begin = 0; \
                    ks->end = __read(ks->f, ks->buf, __bufsize); \
                    if (ks->end == 0) { ks->is_eof = 1; break; } \
                    if (ks->end == -1) { ks->is_eof = 1; return -3; } \
                } else break; \
            } \
            if (delimiter == KS_SEP_LINE) { \
                unsigned char *sep = (unsigned char*)memchr(ks->buf + ks->begin, '\n', ks->end - ks->begin); \
                i = sep != NULL ? sep - ks->buf : ks->end; \
            } else if (delimiter > KS_SEP_MAX) { \
                unsigned char *sep = (unsigned char*)memchr(ks->buf + ks->begin, delimiter, ks->end - ks->begin); \
                i = sep != NULL ? sep - ks->buf : ks->end; \
            } else if (delimiter == KS_SEP_SPACE) { \
                for (i = ks->begin; i < ks->end; ++i) \
                    if (isspace(ks->buf[i])) break; \
            } else if (delimiter == KS_SEP_TAB) { \
                for (i = ks->begin; i < ks->end; ++i) \
                    if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
            } else i = 0; \
            if (str->m - str->l < (size_t)(i - ks->begin + 1)) { \
                str->m = str->l + (i - ks->begin) + 1; \
                str->s = (char*)realloc(str->s, str->m); \
            } \
            gotany = 1; \
            memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
            str->l = str->l + (i - ks->begin); \
            ks->begin = i + 1; \
            if (i < ks->end) { \
                if (dret) *dret = ks->buf[i]; \
                break; \
            } \
        } \
        if (!gotany && ks_eof(ks)) return -1; \
        if (str->s == 0) { \
            str->m = 1; \
            str->s = (char*)calloc(1, 1); \
        } else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
        str->s[str->l] = '\0'; \
        return str->l; \
    } \
    static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) { \
        return ks_getuntil2(ks, delimiter, str, dret, 0); \
    }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
    __KS_TYPE(type_t) \
    __KS_BASIC(type_t, __bufsize) \
    __KS_GETC(__read, __bufsize) \
    __KS_GETUNTIL(__read, __bufsize)

#define __KSEQ_TYPE(type_t) \
    typedef struct { \
        kstring_t name, comment, seq, qual; \
        int last_char; \
        kstream_t *f; \
    } kseq_t;

#define __KSEQ_BASIC(SCOPE, type_t) \
    SCOPE kseq_t *kseq_init(type_t fd) { \
        kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t)); \
        s->f = ks_init(fd); \
        return s; \
    } \
    SCOPE void kseq_destroy(kseq_t *ks) { \
        if (!ks) return; \
        free(ks->name.s); free(ks->comment.s); free(ks->seq.s); free(ks->qual.s); \
        ks_destroy(ks->f); \
        free(ks); \
    }

#define __KSEQ_READ(SCOPE) \
    SCOPE int kseq_read(kseq_t *seq) { \
        int c,r; \
        kstream_t *ks = seq->f; \
        if (seq->last_char == 0) { \
            while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '@'); \
            if (c < 0) return c; \
            seq->last_char = c; \
        } \
        seq->comment.l = seq->seq.l = seq->qual.l = 0; \
        if ((r=ks_getuntil(ks, 0, &seq->name, &c)) < 0) return r; \
        if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); \
        if (seq->seq.s == 0) { \
            seq->seq.m = 256; \
            seq->seq.s = (char*)malloc(seq->seq.m); \
        } \
        while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '+' && c != '@') { \
            if (c == '\n') continue; \
            seq->seq.s[seq->seq.l++] = c; \
            ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); \
        } \
        if (c == '>' || c == '@') seq->last_char = c; \
        if (seq->seq.l + 1 >= seq->seq.m) { \
            seq->seq.m = seq->seq.l + 2; \
            seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
        } \
        seq->seq.s[seq->seq.l] = 0; \
        if (c != '+') return seq->seq.l; \
        if (seq->qual.m < seq->seq.m) { \
            seq->qual.m = seq->seq.m; \
            seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
        } \
        while ((c = ks_getc(ks)) >= 0 && c != '\n'); \
        if (c == -1) return -2; \
        while ((c = ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l)); \
        if (c == -3) return -3; \
        seq->last_char = 0; \
        if (seq->seq.l != seq->qual.l) return -2; \
        return seq->seq.l; \
    }

#define KSEQ_INIT2(SCOPE, type_t, __read) \
    KSTREAM_INIT(type_t, __read, 16384) \
    __KSEQ_TYPE(type_t) \
    __KSEQ_BASIC(SCOPE, type_t) \
    __KSEQ_READ(SCOPE)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read)

#define KSEQ_DECLARE(type_t) \
    __KS_TYPE(type_t) \
    __KSEQ_TYPE(type_t) \
    extern kseq_t *kseq_init(type_t fd); \
    void kseq_destroy(kseq_t *ks); \
    int kseq_read(kseq_t *seq);
