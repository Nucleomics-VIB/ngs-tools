
#define _POSIX_C_SOURCE 200112L

// filter_N_seq.c
// Description: Split a fasta file into two files:
// - one with sequences containing only N/./- (empty)
// - one with all other sequences (noempty).
// Usage: filter_N_seq -i input.fasta
//
// Requires: kseq.h (https://github.com/attractivechaos/klib)
//           zlib (for file reading)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // Required for getopt/optarg on Ubuntu
#include <stdbool.h>
#include <zlib.h>

#include "lib/kseq.h"
KSEQ_INIT(gzFile, gzread)

// Returns true if the sequence contains only N, ., or -
bool is_empty_seq(const char* seq) {
    for (const char* p = seq; *p; ++p) {
        if (*p != 'N' && *p != '.' && *p != '-') {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    char* infile = NULL;
    int opt;
    while ((opt = getopt(argc, argv, "i:")) != -1) {
        switch (opt) {
            case 'i': infile = optarg; break;
            default:
                fprintf(stderr, "Usage: %s -i input.fasta\n", argv[0]);
                return 1;
        }
    }
    if (!infile) {
        fprintf(stderr, "Usage: %s -i input.fasta\n", argv[0]);
        return 1;
    }

    char empty_out[4096], noempty_out[4096];
    snprintf(empty_out, sizeof(empty_out), "%s.empty.fasta", infile);
    snprintf(noempty_out, sizeof(noempty_out), "%s.noempty.fasta", infile);

    FILE* f_empty = fopen(empty_out, "w");
    FILE* f_noempty = fopen(noempty_out, "w");
    if (!f_empty || !f_noempty) {
        fprintf(stderr, "ERROR: Cannot open output files.\n");
        return 2;
    }

    gzFile fp = gzopen(infile, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open input file %s\n", infile);
        return 3;
    }
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        if (is_empty_seq(seq->seq.s)) {
            fprintf(f_empty, ">%s", seq->name.s);
            if (seq->comment.l) fprintf(f_empty, " %s", seq->comment.s);
            fprintf(f_empty, "\n%s\n", seq->seq.s);
        } else {
            fprintf(f_noempty, ">%s", seq->name.s);
            if (seq->comment.l) fprintf(f_noempty, " %s", seq->comment.s);
            fprintf(f_noempty, "\n%s\n", seq->seq.s);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    fclose(f_empty);
    fclose(f_noempty);
    return 0;
}
