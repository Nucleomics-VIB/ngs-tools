#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <zlib.h>

// source code: replace_quality_scores_gt.cpp 
// SP@NC, 2023-10-19, version 1.5
// read from a fastq (text or gzipped)
// write to stdout (-o -), or to a .fq or a .fq.gz depending on the presence of -z
// replace all base quality scores greater than value given by -m with that value
// compile me with: g++ -o replace_quality_scores_gt replace_quality_scores_gt.cpp -lz

int main(int argc, char* argv[]) {
    int maxQualityScore = 0;
    const char* inputFilename = nullptr;
    const char* outputFilename = nullptr;
    bool gzipOutput = false;
    bool decompressInput = false;

    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-i" && i + 1 < argc) {
            inputFilename = argv[i + 1];
            i++; // Consume the next argument
        } else if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            outputFilename = argv[i + 1];
            i++; // Consume the next argument
        } else if (std::string(argv[i]) == "-m" && i + 1 < argc) {
            maxQualityScore = std::atoi(argv[i + 1]);
            i++; // Consume the next argument
        } else if (std::string(argv[i]) == "-z") {
            gzipOutput = true;
        } else {
            std::cerr << "Error: Invalid option or missing value: " << argv[i] << std::endl;
            return 1;
        }
    }

    if (inputFilename == nullptr || maxQualityScore <= 0 || maxQualityScore > 94) {
        std::cerr << "Usage: " << argv[0] << " -i <input.fastq.gz> -o <output_prefix (- for stdout)> -m <max_quality_score> [-z]" << std::endl;
        return 1;
    }

    // Check if the input file ends with ".gz" to determine if decompression is needed
    if (std::string(inputFilename).rfind(".gz") == std::string(inputFilename).length() - 3) {
        decompressInput = true;
    }

    gzFile input = decompressInput ? gzopen(inputFilename, "rb") : gzopen(inputFilename, "r");
    gzFile gzOutput = nullptr;
    FILE* output = nullptr;

    if (outputFilename != nullptr) {
        if (std::string(outputFilename) == "-") {
            output = stdout;
        } else if (gzipOutput) {
            std::string outputFileName(outputFilename);
            outputFileName += ".fq.gz";  // Append .gz to the output file name
            gzOutput = gzopen(outputFileName.c_str(), "wb");
        } else {
            std::string outputFileName(outputFilename);
            if (outputFileName.rfind(".fq", outputFileName.length() - 3)) {
                outputFileName += ".fq";  // Append .fq to the output file name
            }
            output = fopen(outputFileName.c_str(), "w");
        }
    }

    if (input == NULL) {
        std::cerr << "Error: Could not open input file: " << inputFilename << std::endl;
        return 1;
    }

    if (output == NULL && gzOutput == NULL) {
        std::cerr << "Error: Could not open output file: " << (outputFilename != nullptr ? outputFilename : "stdout") << std::endl;
        return 1;
    }

    char line[2048]; // Adjust the buffer size as needed
    int lineCount = 0;
    while (gzgets(input, line, sizeof(line))) {
        if (lineCount % 4 == 3) {
            // Quality scores line
            for (int i = 0; line[i] != '\0'; i++) {
                if (line[i] > char(maxQualityScore + 33)) {
                    line[i] = char(maxQualityScore + 33);  // Convert to Phred quality score
                }
            }
        }

        if (gzOutput) {
            gzputs(gzOutput, line);
        } else {
            fputs(line, output);
        }
        lineCount++;
    }

    gzclose(input);
    if (gzOutput) {
        gzclose(gzOutput);
    } else if (outputFilename != nullptr && std::string(outputFilename) != "-") {
        fclose(output);
    }

    return 0;
}
