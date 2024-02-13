#include <iostream>
#include <sstream>
#include <cmath>

// produce metrics for R-plotting
// start analysing at base 61 (adaptor excluded)
// with plot_read_metrics.R -i fastq_metrics_file
// SP@NC; 2024_01_05 (+GPT)
//
// script: fastq2metrics.cpp
// compile with: g++ -o fastq2metrics60 fastq2metrics60.cpp

// Function to directly subtract 33 from ASCII value
int ord(char c) {
    return static_cast<int>(c) - 33;
}

// Function to calculate average quality
double aveQual(const std::string& quals) {
    double sum = 0;
    int len = quals.length();
    for (int i = 60; i < len; i++) { // Start from the 61st base
        int q = ord(quals[i]);
        sum += pow(10, q / -10.0);
    }
    double avg_prob = sum / (len - 60); // Adjust length
    double avg_phred = -10 * log10(avg_prob);
    return avg_phred;
}

// Function to calculate GC content
double calculateGC(const std::string& sequence) {
    int gc_count = 0;
    for (size_t i = 60; i < sequence.length(); i++) { // Start from the 61st base
        char base = sequence[i];
        if (base == 'G' || base == 'C') {
            gc_count++;
        }
    }
    return static_cast<double>(gc_count) / (sequence.length() - 60); // Adjust length
}

// Function to round decimal numbers
std::string round_decimal(double number, int decimals) {
    std::ostringstream stream;
    stream.precision(decimals);
    stream << std::fixed << number;
    return stream.str();
}

int main() {
    std::cout << "readid\tmeanq\tlength\tgc" << std::endl;

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line[0] == '@') {
            std::string header = line;
            std::string seq, sep, qual;
            std::getline(std::cin, seq);
            std::getline(std::cin, sep);
            std::getline(std::cin, qual);

            // Extract the first field from the space-separated header and remove leading '@'
            std::istringstream headerStream(header);
            std::string read_id;
            headerStream >> read_id;
            read_id.erase(0, 1);

            // Call the aveQual function and print the result
            double avg_phred = aveQual(qual);

            // Call the calculateGC function and print the result
            double gc_content = calculateGC(seq);

            // Output tab-separated values with rounded GC content and quality
            std::cout << read_id << "\t" << round_decimal(avg_phred, 2) << "\t"
                      << seq.length() - 60 << "\t" << round_decimal(gc_content, 2) << std::endl;
        }
    }

    return 0;
}
