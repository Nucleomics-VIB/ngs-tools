// filter_N_seq2.cpp
// Description: Split a fasta file into two files or read from stdin and write to stdout:
// - one with sequences containing only N/./- (empty)
// - one with all other sequences (noempty).
// Usage: filter_N_seq2_cpp -i input.fasta
//        cat input.fasta | filter_N_seq2_cpp
//
// Requires: kseq.h (https://github.com/attractivechaos/klib)
//           zlib (for file reading)

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>

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
    while ((opt = getopt(argc, argv, "i:h")) != -1) {
        switch (opt) {
            case 'i': infile = optarg; break;
            case 'h':
                std::cout << "Usage: " << argv[0] << " [-i input.fasta]" << std::endl;
                std::cout << "  -i FILE    Input FASTA file (default: read from stdin)" << std::endl;
                std::cout << "  -h         Show this help message" << std::endl;
                std::cout << "\nExamples:" << std::endl;
                std::cout << "  " << argv[0] << " -i input.fasta" << std::endl;
                std::cout << "  cat input.fasta | " << argv[0] << std::endl;
                std::cout << "  gzip -dc input.fasta.gz | " << argv[0] << std::endl;
                return 0;
            default:
                std::cerr << "Usage: " << argv[0] << " [-i input.fasta]" << std::endl;
                return 1;
        }
    }

    // Determine if reading from stdin or file
    bool use_stdin = infile.empty();
    gzFile fp;
    
    if (use_stdin) {
        // Read from stdin
        fp = gzdopen(fileno(stdin), "r");
        if (!fp) {
            std::cerr << "ERROR: Cannot open stdin" << std::endl;
            return 3;
        }
    } else {
        // Read from file
        fp = gzopen(infile.c_str(), "r");
        if (!fp) {
            std::cerr << "ERROR: Cannot open input file " << infile << std::endl;
            return 3;
        }
    }

    // Setup output: stdout or files
    std::ofstream f_empty_file;
    std::ofstream f_noempty_file;
    std::ostream *f_empty = &std::cout;
    std::ostream *f_noempty = &std::cout;
    
    if (!use_stdin) {
        // Write to separate files when reading from a file
        std::string empty_out = infile + ".empty.fasta";
        std::string noempty_out = infile + ".noempty.fasta";
        
        f_empty_file.open(empty_out);
        f_noempty_file.open(noempty_out);
        
        if (!f_empty_file.is_open() || !f_noempty_file.is_open()) {
            std::cerr << "ERROR: Cannot open output files." << std::endl;
            gzclose(fp);
            return 2;
        }
        
        f_empty = &f_empty_file;
        f_noempty = &f_noempty_file;
    }

    // Process sequences
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        if (use_stdin) {
            // When reading from stdin, write all to stdout with a tag
            if (is_empty_seq(seq->seq.s)) {
                // Skip empty sequences or optionally mark them
                // Uncomment next line to output empty sequences with a tag
                // std::cout << ">" << seq->name.s << " [EMPTY]\n" << seq->seq.s << "\n";
            } else {
                std::cout << ">" << seq->name.s << "\n" << seq->seq.s << "\n";
            }
        } else {
            // When reading from file, split into two files
            if (is_empty_seq(seq->seq.s)) {
                *f_empty << ">" << seq->name.s << "\n" << seq->seq.s << "\n";
            } else {
                *f_noempty << ">" << seq->name.s << "\n" << seq->seq.s << "\n";
            }
        }
    }
    
    kseq_destroy(seq);
    gzclose(fp);
    
    if (!use_stdin) {
        f_empty_file.close();
        f_noempty_file.close();
    }
    
    return 0;
}
