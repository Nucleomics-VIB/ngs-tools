#!/bin/bash

# gepard_plot.sh: create a diagalign DNA plot from two assemblies / references
# get the code from https://github.com/univieCUBE/gepard
# ref: Krumsiek J, Arnold R, Rattei T, Bioinformatics 2007; 23(8): 1026-8. PMID: 17309896
#
# Stephane Plaisance (VIB-NC+BITS) 2017/04/21; v1.1
# added color limits with -l -u -g
#
# visit our Git: https://github.com/Nucleomics-VIB

# check parameters for your system
# Post-install
# create an alias named 'gepard.jar' (pointing to the actual jar file) next to a matrices folder containing edna.mat and other matrices 
# refer with -p to that place
version="1.1, 2017_09_05"

usage='# Usage: gepard_plot.sh -x <reference assembly> -y <draft assembly> -p <path to gepard.jar and matrices>
# script version '${version}'
# [optional: -o <result path/prefix>]
# [optional: -w <word size:10>]
# [optional: -W <window size:0>]
# [optional: -l <lower value% (default:0)>]
# [optional: -u <upper value% (default:100)>]
# [optional: -g <greyscale start value% (default:0)>]
# [optional: -J <java extra parameters (eg -Xmx1G, put between double quotes if it contains spaces)>
# [optional: -h <this help text>]'

while getopts "x:y:w:W:o:p:J:l:u:g:h" opt; do
  case $opt in
    x) reference=${OPTARG} ;;
    y) draftassembly=${OPTARG} ;;
    w) wordopt=${OPTARG} ;;
    W) windowopt=${OPTARG} ;;
    l) lowopt=${OPTARG} ;;
    u) uppopt=${OPTARG} ;;
    g) greyopt=${OPTARG} ;;
    o) outfile=${OPTARG} ;;
    p) gepardpath=${OPTARG} ;;
    J) javaargs=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# test if minimal arguments were provided
if [ -z "${reference}" ]
then
	echo "# no reference assembly provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${reference}" ]; then
    echo "${reference} file not found!";
    exit 1
fi

if [ -z "${draftassembly}" ]
then
   echo "# no draft assembly provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${draftassembly}" ]; then
	echo "${draftassembly} file not found!"
	exit 1
fi

if [ -z "${gepardpath}" ]
then
	echo "# no path to gepard.jar provided!"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${gepardpath}/gepard.jar" ]; then
    echo "gepard.jar file not found at ${gepardpath}!";
    exit 1
fi

# required matrix file
if [ ! -f "${gepardpath}/matrices/edna.mat" ]; then
    echo "edna.mat file not found at ${gepardpath}/matrices!";
    exit 1
fi

# other parameters or defaults
destfile=${outfile:-"gepard-$(basename "${draftassembly}" | cut -d. -f1)_vs_$(basename "${reference}" | cut -d. -f1)"}
optargs="w${wordopt:-10}_W${windowopt:-0}_${lowopt:-0}_${uppopt:-100}_${greyopt:-0}"

# build the command
cmd="java ${javaargs:-""} -cp ${gepardpath}/gepard.jar org.gepard.client.cmdline.CommandLine \
	-seq1 ${reference} \
	-seq2 ${draftassembly} \
	-matrix ${gepardpath}/matrices/edna.mat \
	-word ${wordopt:-10} \
	-window ${windowopt:-0} \
	-lower ${lowopt:-0} \
	-upper ${uppopt:-100} \
	-greyscale ${greyopt:-0} \
	-format png \
	-outfile ${destfile}_${optargs}.png"

# show and execute
echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"

exit 0

# Gepard 1.40 final - command line mode
# 
# Reference:
# Krumsiek J, Arnold R, Rattei T
# Gepard: A rapid and sensitive tool for creating dotplots on genome scale.
# Bioinformatics 2007; 23(8): 1026-8. PMID: 17309896
# 
# Parameters are supplied as -name value
# 
# Required parameters:
#   -seq1:        first sequence file
#   -seq2:        second sequence file
#   -matrix:      substitution matrix file
#   -outfile:     output file name
# 
# Dotplot image parameters:
#   -maxwidth:    maximum width of the generated image (default: 750)
#   -maxheight:   maximum height of the generated image (default: 750)
#   -zoom:        specify a zoom factor for the dotplot
#   note: you can only use maxwidth/maxheight OR zoom
#         when using maxwidth/maxheight the program tries to generate the largest
#         possible dotplot within the given bounds
#   -format:      output format, one of:  'png', 'jpg', 'bmp' (default:PNG)
# 
# Dotplot computation parameters:
#   -secondcomp   use complementary of second sequence
#   -word:        word length for suffix array mode (default: 10)
#   -window:      window size for ordinary dotplot mode (default: 0)
#   if a window value and no word value is specified, word=0 is assumed
# 
# Suffix array parameters:
#   -safile       load suffix array from file instead of calculating it
#   -sasecondseq  the suffix array is for the second sequence
#   if -sasecondseq is NOT specified, the suffix array will be used for first sequence
# 
# Coordinate parameters (absolute values of % values)
#   -from1,-to1   coordinates of first sequence
#   -from2,-to2   coordinates of second sequence
#   if these parameters are not specified the full sequence will be used
# 
# Display parameters:
#   -lower        lower limit for dot intensity (in %)
#   -upper        upper limit for dot intensity (in %)
#   -greyscale    greyscale start value (in %)
# 
# Miscellaneous:
#   -silent       generate no output (except error messages)