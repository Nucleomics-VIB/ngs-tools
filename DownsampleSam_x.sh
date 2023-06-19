#!/usr/bin/bash

# script: DownsampleSam_x.sh
# picard Downsample mappings for benchmarking
# SP@NC; 2023-06-19; v1.0

version="1.1, 2023-06-19"

usage='# Usage: downsample_x.sh
#    -b <input BAM file>
#    -p <percent data to retain (integer number eg. 70 => 70%)>
# optional
#    -s seed <default to 1>
#    -a accuracy <default to 0.001>
#    -S Strategy <default to "ConstantMemory">
# version: '${version}

# Parse command-line options
while getopts "b:p:s:a:S:" opt; do
  case $opt in
    b) opt_bam=$OPTARG ;;
    p) opt_pc=$OPTARG ;;
    s) opt_seed=$OPTARG ;;
    a) opt_acc=$OPTARG ;;
    S) opt_str=$OPTARG ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Check if input BAM file and pc value are provided
if [ -z "$opt_bam" ]; then
  echo "Input BAM file is required (sorted and indexed)"
  exit 1
fi

if [ -z "$opt_pc" ]; then
  echo "Value for -p (pc) option is required (integer number)"
  exit 1
fi

# variables
outbam=${opt_bam%.bam}_${opt_pc}pc_ds.bam
percent=$(echo "scale=2; ${opt_pc}/100" | bc)

# seed and default settings
sn=${opt_seed:-1}
ac=${opt_acc:-0.0001}
st=${opt_str:-"ConstantMemory"}

# work locally
mkdir -p tmpdir

java -jar $PICARD/picard.jar DownsampleSam \
  -I ${opt_bam} \
  -O ${outbam} \
  -R ${sn} \
  -S ${st} \
  -P ${percent} \
  -A ${ac} \
  -M "${opt_bam%.bam}_${opt_pc}pc_ds_metrics.txt" \
  --CREATE_INDEX TRUE \
  --VALIDATION_STRINGENCY LENIENT \
  --VERBOSITY INFO \
  --TMP_DIR tmpdir \
  > "${opt_bam%.bam}_${opt_pc}pc_downsampling_log.txt" 2>&1

# show summary
echo "# summary"
tail -4 "${opt_bam%.bam}_${opt_pc}pc_ds_metrics.txt" | head -2 | transpose -t | column -t
