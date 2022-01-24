#!/bin/bash

# script name: consensusQualFilter.sh
# compute the % UpperCase nucleotides in consensus
# compare to minimal value provided by user
# output only sequences that have a higher percentage of capital nucleotides
# corresponding to better consensus
# report bad sequences to stderr
#
# Stephane Plaisance (VIB-NC) 2021/12/30; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# somehow not doing it completely
# see fastaConsensusFilter.pl for OK variant
exit 1

infile=$1
minpc=${2:-100}

bioawk -c fastx -v minpc="${minpc}" 'BEGIN{min=100*minpc}{
cmd1 = "echo -n "$seq" | grep -o [ATGCN] | tr -d \"\n\" | wc -m";
cmd1 | getline UC;
close(cmd1);
cmd2 = "echo -n "$seq" | grep -o [atgcn] | tr -d \"\n\" | wc -m";
cmd2 | getline LC;
close(cmd2);
cmd3 = "echo \"10000*"UC"/("UC"+"LC")\" | bc";
cmd3 | getline perc;
close(cmd3);
cmd4 = "if [ "perc" -lt "min" ]; then echo bad; fi"
cmd4 | getline test
close(cmd4)
if ( test != "bad" )
  print ">"$name"\n"$seq;
else
  print "## low quality sequence ("sprintf("%.2f%%",perc/100)")\n>"$name"\n"$seq > "/dev/stderr";
}' ${infile}
