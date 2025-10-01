# filter_N_seq

`filter_N_seq` is a tool for splitting FASTA sequence files into two files:
- one with sequences made exclusively of undefined bases (N, . or -)
- one with all other sequences (containing at least one defined base)

The main functionality is provided by the Bash script `filter_N_seq.sh`, which efficiently processes input files and separates sequences based on whether they are made up only of undefined bases or not.

This is particularly useful for downstream analyses that require unambiguous sequence data. The script is simple, portable, and easy to use on any Unix-like system.

To further improve performance and flexibility, a C/C++ implementation is being developed as a drop-in replacement for the script, especially for large datasets.

## What does the shell script do?

The `filter_N_seq.sh` script reads a sequence file line by line, and for each sequence, checks if it is made up exclusively of undefined bases (N, . or -). Sequences that are only undefined bases are written to the `.empty.fasta` output, while all other sequences (containing at least one defined base) are written to the `.noempty.fasta` output. This is designed to be fast and memory-efficient, making it suitable for large-scale sequence filtering tasks.

## Project Structure

- `filter_N_seq.sh` — Bash script for filtering N sequences (main tool)
- `filter_N_seq.c` — C implementation of the filter logic (performance improvement, in progress or as reference)
- `filter_N_seq2.cpp` — C++ implementation with stdin/pipe support for use in pipelines
- `klib/` — External C utility library (used for efficient data structures and algorithms in the C/C++ version)
- `klib/test/` — Test files for klib components

## Usage

### To use the Bash script (recommended for most users):

```bash
bash filter_N_seq.sh -i input.fasta
```

### C/C++ version (for large datasets or advanced users):

After compiling, run:

```bash
./filter_N_seq -i input.fasta
```

This will produce:
- `input.fasta.empty.fasta` (sequences with only N/./-)
- `input.fasta.noempty.fasta` (all other sequences)

### Version 2 with pipe support (for use in pipelines):

The `filter_N_seq2.cpp` version is designed for piping and can read from stdin:

```bash
# Pipe from standard input - outputs only non-empty sequences to stdout
cat input.fasta | ./filter_N_seq2_cpp

# Works with compressed files
gzip -dc input.fasta.gz | ./filter_N_seq2_cpp

# Can be used in a pipeline
cat input.fasta | ./filter_N_seq2_cpp | some_other_tool

# Still supports file mode to create both output files
./filter_N_seq2_cpp -i input.fasta
```

**Note:** When using pipes (stdin), version 2 only outputs **non-empty sequences** to stdout. Empty sequences (containing only N/./-)  are **discarded** and not saved. Use the `-i` flag if you need both output files.

## Getting the klib Library

The C/C++ version uses [klib](https://github.com/attractivechaos/klib) for efficient data structures and algorithms. To use it, you need to download the library:

```bash
# Clone klib into the project directory (if not already present)
git clone https://github.com/attractivechaos/klib.git
```

The strictly necessary libraries were copied from this git version and placed in the lib folder to build this project (2025-09-30).

## Building (C++ version)

To compile the C++ version (example):

```bash
# Compile filter_N_seq.cpp with klib (adjust filenames as needed)
mkdir -p bin
g++ -O2 -std=c++11 -I./lib -o bin/filter_N_seq_cpp filter_N_seq.cpp -lz
```

**Note for Mac (Apple Silicon):**
If you are compiling the C version (`filter_N_seq.c`), use `gcc` instead of `g++`:

```bash
mkdir -p bin
gcc -O2 -std=c11 -I./lib -o bin/filter_N_seq_c filter_N_seq.c -lz
```

If you use a Makefile, ensure to add `-I./klib` to your `CXXFLAGS` so the compiler can find the klib headers.

## Note for C++ compatibility with klib/kseq.h

If you want to compile the C++ version, you must patch `klib/kseq.h` to add an explicit cast for C++ compatibility. Edit the following lines in `klib/kseq.h`:

Replace:
```c
unsigned char *sep = memchr(ks->buf + ks->begin, '\n', ks->end - ks->begin);
unsigned char *sep = memchr(ks->buf + ks->begin, delimiter, ks->end - ks->begin);
```
With:
```c
unsigned char *sep = (unsigned char*)memchr(ks->buf + ks->begin, '\n', ks->end - ks->begin);
unsigned char *sep = (unsigned char*)memchr(ks->buf + ks->begin, delimiter, ks->end - ks->begin);
```
This ensures the code compiles with both gcc and g++. Note that the edits were applied to the version copied to the lib folder

## License

See klib licensing info in `lib/README.md` and `lib/LICENSE.txt` for third-party library licensing.

## Authors

- SP@Nucleomics-VIB; 2025-09-30; v1.0
