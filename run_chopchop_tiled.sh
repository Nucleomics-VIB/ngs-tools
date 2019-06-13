#!/bin/bash
# script name: run_chopchop_tiled.sh
# run chopchop prediction upstream and downstream of a region
# and add guides every Nkb within the region
#
# Stephane Plaisance (VIB-NC) 2019/06/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2019_06_12"

usage='# Usage: run_chopchop_tiled.sh
# -r <region (eg. chr1:16740273-16972964)>
# -s <steps within region (default to 20000)>
# -w <width for prediction at both ends (default to 3000)>
# -m <median DNA fragment size (default 30000)>
# script version '${version}'
# [-h for this help]'

while getopts "r:w:m:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    s) opts=${OPTARG} ;;
    w) optw=${OPTARG} ;;
    m) optm=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# add all logs to a global log file
datetag=$(date +%s)
exec &> >(tee -a "run_chopchop_tiled_${datetag}.log")

# provide the region to be captured (start to end shorter than 20kb)
region=${optr:-"chr1:16740273-16972964"}
chr=${region%%:*}
coord=${region#*:}
start=${coord%-*}
end=${coord#*-}

# define window and step sizes
window=${optw:-3000}
steps=${opts:-30000}

# check region width and compare to median molecule size limit
width=$((${end}-${start}))
median=${optm:-30000}
if [ "${width}" -lt "${median}" ]; then
echo "# The region is too short (${width} bps) for tiled sgRNA pairs, consider running run_chopchop.sh !"
exit 0
fi

# create genome regions in steps of $steps
bedtools makewindows -b <(echo -e ${chr}$'\t'${start}$'\t'${end}) -w ${steps} > tmp.regions.bed

# default primer3 settings
primer3args='PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60'
backbone="AGGCTAGTCCGT"
scoring="DOENCH_2014"

# loop over regions and predict '+' and '-' strands
while read rchr rstart rend; do

    # tile limits
	tile="${rchr}_${rstart}_${rend}"

	# create index name and folders for results
	pfx="run_chopchop_tiled_${datetag}_data/${rchr}_${rstart}_${rend}_results"
	mkdir -p ${pfx}

	# define upstream & downstream
	upstream="${rchr}:$((${rstart}-${window}))-${rstart}"
	downstream="${rchr}:${rend}-$((${rend}+${window}))"
	
	# search 3kb upstream
	echo "# predicting ${window}bps upstream of ${tile}"
	chopchop.py -J \
			-P \
			-T 1 \
			-M NGG \
			--maxMismatches 3 \
			-g 20 \
			--scoringMethod ${scoring} \
			-f NN \
			--backbone ${backbone} \
			--replace5P GG \
			-G hg38 \
			-t WHOLE \
			-n N \
			-R 4 \
			-3 ${primer3args} \
			-A 290 \
			-a 20 \
			--rm1perfOff \
			-o ${pfx}/upstream/ \
			--filterGCmin 40 \
			--filterGCmax 80 \
			--filterSelfCompMax 0 \
			--BED \
			--GenBank \
			-Target ${upstream} \
			> ${pfx}/upstream_results.txt 2> ${pfx}/upstream_python.err

	# filter best hits on + strand   
	gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/+/ && $6==0 && $11>0.3)) print $0}' \
		${pfx}/upstream_results.txt | column -s $'\t' -t \
		> ${pfx}/upstream_results-filtered.txt
	echo "# => found $(($(wc -l ${pfx}/upstream_results.txt|cut -d" " -f1)-1)) \
	candidates of which $(($(wc -l ${pfx}/upstream_results-filtered.txt|cut -d" " -f1)-1)) HQ hits"

	# search 3kb downstream
	echo "# predicting ${window}bps downstream of ${tile}"
	chopchop.py -J \
			-P \
			-T 1 \
			-M NGG \
			--maxMismatches 3 \
			-g 20 \
			--scoringMethod ${scoring} \
			-f NN \
			--backbone ${backbone} \
			--replace5P GG \
			-G hg38 \
			-t WHOLE \
			-n N \
			-R 4 \
			-3 ${primer3args} \
			-A 290 \
			-a 20 \
			--rm1perfOff \
			-o ${pfx}/downstream/ \
			--filterGCmin 40 \
			--filterGCmax 80 \
			--filterSelfCompMax 0 \
			--BED \
			--GenBank \
			-Target ${downstream} \
			> ${pfx}/downstream_results.txt 2> ${pfx}/downstream_python.err

	# filter best hits on + strand   
	gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/-/ && $6==0 && $11>0.3)) print $0}' \
		${pfx}/downstream_results.txt | column -s $'\t' -t \
		> ${pfx}/downstream_results-filtered.txt
	echo "# => found $(($(wc -l ${pfx}/downstream_results.txt|cut -d" " -f1)-1)) \
	candidates of which $(($(wc -l ${pfx}/downstream_results-filtered.txt|cut -d" " -f1)-1)) HQ hits"
	
done < tmp.regions.bed



