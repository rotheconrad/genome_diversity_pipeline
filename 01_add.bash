#!/bin/bash

target=$1
dataset=$2
file1=$3
file2=$4
scr=$(readlink -f "$0" 2>/dev/null)
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset || ! -n $file1 ]] ; then
  echo "Usage: $0 target_folder dataset file1 [file2]
  target_folder    Path to the target folder containing the data
  dataset          Name of the dataset to create
  file1            Path to the FWD (or unpaired) reads, in .fastq.gz format
  file2 (optional) Path to the REV reads, in .fastq.gz format
  "
  exit 0
fi

GDIV_PKG="$pkg" . "$pkg/00_env.bash"

"$pkg/00_build.bash" "$target"
cp "$file1" "$target/01_reads/${dataset}.1.fastq.gz"
[[ -n $file2 ]] && cp "$file2" "$target/01_reads/${dataset}.2.fastq.gz"
cd "$target"

SIZE=$(ls -pl "01_reads/${dataset}.1.fastq.gz" | awk '{print $5}')
if [[ -n $file2 ]] ; then
  SIZE2=$(ls -pl "01_reads/${dataset}.2.fastq.gz" | awk '{print $5}')
  let SIZE=$SIZE+$SIZE2
fi
let SIZE_G=$SIZE/1000000000
let TIME_H=2+$SIZE_G

# Launch next step
qsub "$pkg/00_launcher.pbs" -N "GD02-$dataset" \
  -v "PKG=$pkg,TARGET=$target,DATASET=$dataset,STEP=02_trim" \
  -l nodes=1:ppn=1 -l mem=50g -l "walltime=$TIME_H:00:00" \
  -o "xx_log/${dataset}.01.txt" -j oe

