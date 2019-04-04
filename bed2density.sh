#!/bin/bash

# create windows from a .genome file 
# count window-overlapping features from a second file
# report as a IGV file

genome=${1:-"GCF_003713225.1_Cara_1.0.genome"}
width=${2:-100}
input=${3:-"group_H.bed"}
stat="count"

# create windows from genome
winfile=${genome}_${width}k_windows.bed
bedtools makewindows -w ${width}000 -g ${genome} > ${winfile}

# intersect data and windows
bedtools intersect -a ${winfile} -b ${input} -wa -wb > ${input%.bed}_${width}k_windows.txt && \
  rm ${winfile}

# create density track
lastcol=7
featname=${width}k_density
echo "#track name=${input%.bed}_${width}k_density autoScale=on graphType=bar color=255,0,0" \
  > ${input%.bed}_${width}k_density.igv
bedtools groupby -i ${input%.bed}_${width}k_windows.txt -g 1,2,3 -c ${lastcol} -o ${stat} | \
  awk -v featname=${featname} 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3,featname,$4}' \
    >> ${input%.bed}_${width}k_density.igv && \
    rm ${input%.bed}_${width}k_windows.txt
