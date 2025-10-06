// filter_N_seq.cpp
// Description: Split a fasta file into two files:
// - one with sequences containing only N/./- (empty)
// - one with all other sequences (noempty).
// Usage: filter_N_seq_cpp -i input.fasta
//
// Requires: kseq.h (https://github.com/attractivechaos/klib)
//           zlib (for file reading)

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>
#include <zlib.h>

extern "C" {
#include "lib/kseq.h"
}
KSEQ_INIT(gzFile, gzread)

// Returns true if the sequence contains only N, ., or -
bool is_empty_seq(const std::string& seq) {
    for (char c : seq) {
        if (c != 'N' && c != '.' && c != '-') {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::string infile;
    int opt;
    while ((opt = getopt(argc, argv, "i:")) != -1) {
        switch (opt) {
            case 'i': infile = optarg; break;
            default:
                std::cerr << "Usage: " << argv[0] << " -i input.fasta" << std::endl;
                return 1;
        }
    }
    if (infile.empty()) {
        std::cerr << "Usage: " << argv[0] << " -i input.fasta" << std::endl;
        return 1;
    }

    std::string empty_out = infile + ".empty.fasta";
    std::string noempty_out = infile + ".noempty.fasta";

    std::ofstream f_empty(empty_out);
    std::ofstream f_noempty(noempty_out);
    if (!f_empty.is_open() || !f_noempty.is_open()) {
        std::cerr << "ERROR: Cannot open output files." << std::endl;
        return 2;
    }

    gzFile fp = gzopen(infile.c_str(), "r");
    if (!fp) {
        std::cerr << "ERROR: Cannot open input file " << infile << std::endl;
        return 3;
    }
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        if (is_empty_seq(seq->seq.s)) {
            f_empty << ">" << seq->name.s;
            if (seq->comment.l) f_empty << " " << seq->comment.s;
            f_empty << "\n" << seq->seq.s << "\n";
        } else {
            f_noempty << ">" << seq->name.s;
            if (seq->comment.l) f_noempty << " " << seq->comment.s;
            f_noempty << "\n" << seq->seq.s << "\n";
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    f_empty.close();
    f_noempty.close();
    return 0;
}