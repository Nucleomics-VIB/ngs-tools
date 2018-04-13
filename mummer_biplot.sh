#!/bin/bash

# mummer_biplot.sh: produce a pairwise plot from two fasta sequences
#
# Requirements:
# run on a unix computer installed with mummer3 (mummer apps in $PATH)
# two related fasta references to be compared
#
# Stephane Plaisance (VIB-NC+BITS) 2017/09/22; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
version="1.01, 2018_04_13"

usage='# Usage: mummer_biplot.sh -x <reference assembly> -y <query assembly>
# script version '${version}'
# [optional: -o <result folder>]
# [optional: -c <min-cluster|100>]
# [optional: -t <data type (nucmer|promer; default nucmer)>]
# [optional: -I <min-identity to include in show-coords|95>]
# [optional: -L <min-align length to include in show-coords|100>]
# [optional: -f <output format (png,postscript,x11)|png>]
# [optional: -h <this help text>]'

while getopts "x:y:o:c:p:t:I:L:f:h" opt; do
  case $opt in
    x) assembly1=${OPTARG} ;;
    y) assembly2=${OPTARG} ;;
    p) mummerpath=${OPTARG} ;;
    o) outpathopt=${OPTARG} ;;
    f) format=${OPTARG} ;;
    c) clust=${OPTARG} ;;
    t) datatype=${OPTARG} ;;
    I) minidentityopt=${OPTARG} ;;
    L) minalignopt=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
cluster=${clust:-100}
coordfilter=""
outformat=${format:-"png"}

# filtering options for show-coords
minidentity=${minidentityopt:-95}
coordfilter="${coordfilter}_filtered-I.${minidentity}"

minalign=${minalignopt:-100}
coordfilter="${coordfilter}_filtered-L.${minalign}"

# choose protein or nucleic alignment method
if [ -z "${datatype}" ]; then 
    prog="nucmer"
else 
    prog="${datatype}"
fi

# test if minimal arguments were provided
if [ -z "${assembly1}" ]
then
   echo "# no first assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${assembly1}" ]; then
	echo "${assembly1} file not found!"
	exit 1
fi

if [ -z "${assembly2}" ]
then
	echo "#  no second assembly provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${assembly2}" ]; then
    echo "${assembly2} file not found!";
    exit 1
fi

# check if mummer requirements are present
$( hash ${prog} 2>/dev/null ) || ( echo "# ${prog} not found in PATH (nucmer or promer?)"; exit 1 )
$( hash show-coords 2>/dev/null ) || ( echo "# show-coords not found in PATH"; exit 1 )
$( hash mummerplot 2>/dev/null ) || ( echo "# mummer-plot not found in PATH"; exit 1 )

# labels from filenames
xlabel=$(basename ${assembly1%.f*})
ylabel=$(basename ${assembly2%.f*})

# other parameters or defaults
outpath=${outpathopt:-"mummer_results"}
mkdir -p ${outpath}

result="${outpath}/${prog}-plot-${ylabel%.f*}_vs_${xlabel%.f*}"

# build the command
cmd="${prog} --maxmatch \
	-c ${cluster} \
	-p ${result} \
	${assembly1} ${assembly2}
	> ${ylabel}_vs_${xlabel}_mummer3-log.txt 2>&1"

# show and execute	
echo "# ${cmd}"
eval ${cmd}
 
# after success create alignment file and plot
if [ $? -eq 0 ]; then
	cmd="(show-coords -r -c -l ${result}.delta > ${result}_all_coords.txt && \
		show-coords -r -c -l -I ${minidentity} -L ${minalign} ${result}.delta \
		> ${result}${coordfilter}_coords.txt && \
		mummerplot --fat --filter --layout --${outformat} --large -p ${result} ${result}.delta) \
		>> mummer3-log.txt 2>&1"

#		-r ${xlabel} \
#		-q ${ylabel} \

	echo "# ${cmd}"
	eval ${cmd}
else
    echo "Mummer analysis seems to have failed, please check mummer3-log.txt!"
fi

exit 0

########################################################################################
# man pages for the main executables used above

#   USAGE: nucmer  [options]  <Reference>  <Query>
# 
#   DESCRIPTION:
#     nucmer generates nucleotide alignments between two mutli-FASTA input
#     files. The out.delta output file lists the distance between insertions
#     and deletions that produce maximal scoring alignments between each
#     sequence. The show-* utilities know how to read this format.
# 
#   MANDATORY:
#     Reference       Set the input reference multi-FASTA filename
#     Query           Set the input query multi-FASTA filename
# 
#   OPTIONS:
#     --mum           Use anchor matches that are unique in both the reference
#                     and query
#     --mumcand       Same as --mumreference
#     --mumreference  Use anchor matches that are unique in in the reference
#                     but not necessarily unique in the query (default behavior)
#     --maxmatch      Use all anchor matches regardless of their uniqueness
# 
#     -b|breaklen     Set the distance an alignment extension will attempt to
#                     extend poor scoring regions before giving up (default 200)
#     --[no]banded    Enforce absolute banding of dynamic programming matrix
#                     based on diagdiff parameter EXPERIMENTAL (default no)
#     -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
#     --[no]delta     Toggle the creation of the delta file (default --delta)
#     --depend        Print the dependency information and exit
#     -D|diagdiff     Set the maximum diagonal difference between two adjacent
#                     anchors in a cluster (default 5)
#     -d|diagfactor   Set the maximum diagonal difference between two adjacent
#                     anchors in a cluster as a differential fraction of the gap
#                     length (default 0.12)
#     --[no]extend    Toggle the cluster extension step (default --extend)
#     -f
#     --forward       Use only the forward strand of the Query sequences
#     -g|maxgap       Set the maximum gap between two adjacent matches in a
#                     cluster (default 90)
#     -h
#     --help          Display help information and exit
#     -l|minmatch     Set the minimum length of a single match (default 20)
#     -o
#     --coords        Automatically generate the original NUCmer1.1 coords
#                     output file using the 'show-coords' program
#     --[no]optimize  Toggle alignment score optimization, i.e. if an alignment
#                     extension reaches the end of a sequence, it will backtrack
#                     to optimize the alignment score instead of terminating the
#                     alignment at the end of the sequence (default --optimize)
#     -p|prefix       Set the prefix of the output files (default "out")
#     -r
#     --reverse       Use only the reverse complement of the Query sequences
#     --[no]simplify  Simplify alignments by removing shadowed clusters. Turn
#                     this option off if aligning a sequence to itself to look
#                     for repeats (default --simplify)
#     -V
#     --version       Display the version information and exit

#   USAGE: promer  [options]  <Reference>  <Query>
# 
#   DESCRIPTION:
#     promer generates amino acid alignments between two mutli-FASTA DNA input
#     files. The out.delta output file lists the distance between insertions
#     and deletions that produce maximal scoring alignments between each
#     sequence. The show-* utilities know how to read this format. The DNA
#     input is translated into all 6 reading frames in order to generate the
#     output, but the output coordinates reference the original DNA input.
# 
#   MANDATORY:
#     Reference       Set the input reference multi-FASTA DNA file
#     Query           Set the input query multi-FASTA DNA file
# 
#   OPTIONS:
#     --mum           Use anchor matches that are unique in both the reference
#                     and query
#     --mumcand       Same as --mumreference
#     --mumreference  Use anchor matches that are unique in in the reference
#                     but not necessarily unique in the query (default behavior)
#     --maxmatch      Use all anchor matches regardless of their uniqueness
# 
#     -b|breaklen     Set the distance an alignment extension will attempt to
#                     extend poor scoring regions before giving up, measured in
#                     amino acids (default 60)
#     -c|mincluster   Sets the minimum length of a cluster of matches, measured in
#                     amino acids (default 20)
#     --[no]delta     Toggle the creation of the delta file (default --delta)
#     --depend        Print the dependency information and exit
#     -d|diagfactor   Set the clustering diagonal difference separation factor
#                     (default .11)
#     --[no]extend    Toggle the cluster extension step (default --extend)
#     -g|maxgap       Set the maximum gap between two adjacent matches in a
#                     cluster, measured in amino acids (default 30)
#     -h
#     --help          Display help information and exit.
#     -l|minmatch     Set the minimum length of a single match, measured in amino
#                     acids (default 6)
#     -m|masklen      Set the maximum bookend masking lenth, measured in amino
#                     acids (default 8)
#     -o
#     --coords        Automatically generate the original PROmer1.1 ".coords"
#                     output file using the "show-coords" program
#     --[no]optimize  Toggle alignment score optimization, i.e. if an alignment
#                     extension reaches the end of a sequence, it will backtrack
#                     to optimize the alignment score instead of terminating the
#                     alignment at the end of the sequence (default --optimize)
# 
#     -p|prefix       Set the prefix of the output files (default "out")
#     -V
#     --version       Display the version information and exit
#     -x|matrix       Set the alignment matrix number to 1 [BLOSUM 45], 2 [BLOSUM
#                     62] or 3 [BLOSUM 80] (default 2)

# USAGE: show-coords  [options]  <deltafile>
# 
# -b          Merges overlapping alignments regardless of match dir
#             or frame and does not display any idenitity information.
# -B          Switch output to btab format
# -c          Include percent coverage information in the output
# -d          Display the alignment direction in the additional
#             FRM columns (default for promer)
# -g          Deprecated option. Please use 'delta-filter' instead
# -h          Display help information
# -H          Do not print the output header
# -I float    Set minimum percent identity to display
# -k          Knockout (do not display) alignments that overlap
#             another alignment in a different frame by more than 50%
#             of their length, AND have a smaller percent similarity
#             or are less than 75% of the size of the other alignment
#             (promer only)
# -l          Include the sequence length information in the output
# -L long     Set minimum alignment length to display
# -o          Annotate maximal alignments between two sequences, i.e.
#             overlaps between reference and query sequences
# -q          Sort output lines by query IDs and coordinates
# -r          Sort output lines by reference IDs and coordinates
# -T          Switch output to tab-delimited format

#   USAGE: mummerplot  [options]  <match file>
# 
#   DESCRIPTION:
#     mummerplot generates plots of alignment data produced by mummer, nucmer,
#     promer or show-tiling by using the GNU gnuplot utility. After generating
#     the appropriate scripts and datafiles, mummerplot will attempt to run
#     gnuplot to generate the plot. If this attempt fails, a warning will be
#     output and the resulting .gp and .[frh]plot files will remain so that the
#     user may run gnuplot independently. If the attempt succeeds, either an x11
#     window will be spawned or an additional output file will be generated
#     (.ps or .png depending on the selected terminal). Feel free to edit the
#     resulting gnuplot script (.gp) and rerun gnuplot to change line thinkness,
#     labels, colors, plot size etc.
# 
#   MANDATORY:
#     match file      Set the alignment input to 'match file'
#                     Valid inputs are from mummer, nucmer, promer and
#                     show-tiling (.out, .cluster, .delta and .tiling)
# 
#   OPTIONS:
#     -b|breaklen     Highlight alignments with breakpoints further than
#                     breaklen nucleotides from the nearest sequence end
#     --[no]color     Color plot lines with a percent similarity gradient or
#                     turn off all plot color (default color by match dir)
#                     If the plot is very sparse, edit the .gp script to plot
#                     with 'linespoints' instead of 'lines'
#     -c
#     --[no]coverage  Generate a reference coverage plot (default for .tiling)
#     --depend        Print the dependency information and exit
#     -f
#     --filter        Only display .delta alignments which represent the "best"
#                     hit to any particular spot on either sequence, i.e. a
#                     one-to-one mapping of reference and query subsequences
#     -h
#     --help          Display help information and exit
#     -l
#     --layout        Layout a .delta multiplot in an intelligible fashion,
#                     this option requires the -R -Q options
#     --fat           Layout sequences using fattest alignment only
#     -p|prefix       Set the prefix of the output files (default 'out')
#     -rv             Reverse video for x11 plots
#     -r|IdR          Plot a particular reference sequence ID on the X-axis
#     -q|IdQ          Plot a particular query sequence ID on the Y-axis
#     -R|Rfile        Plot an ordered set of reference sequences from Rfile
#     -Q|Qfile        Plot an ordered set of query sequences from Qfile
#                     Rfile/Qfile Can either be the original DNA multi-FastA
#                     files or lists of sequence IDs, lens and dirs [ /+/-]
#     -r|rport        Specify the port to send reference ID and position on
#                     mouse double click in X11 plot window
#     -q|qport        Specify the port to send query IDs and position on mouse
#                     double click in X11 plot window
#     -s|size         Set the output size to small, medium or large
#                     --small --medium --large (default 'small')
#     -S
#     --SNP           Highlight SNP locations in each alignment
#     -t|terminal     Set the output terminal to x11, postscript or png
#                     --x11 --postscript --png (default 'x11')
#     -t|title        Specify the gnuplot plot title (default none)
#     -x|xrange       Set the xrange for the plot '[min:max]'
#     -y|yrange       Set the yrange for the plot '[min:max]'
#     -V
#     --version       Display the version information and exit