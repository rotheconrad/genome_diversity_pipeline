#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset ]] ; then
  echo "Usage: $0 target_folder dataset
  target_folder Path to the target folder containing the data
  dataset       Name of the dataset to process
  "
  exit 0
fi

. "$pkg/00_env.bash"
RAM=${RAM:-}
cd "$target"
cmd="bbnorm.sh in='02_trim/${dataset}.1.fastq.gz' \
  out='03_norm/${dataset}.1.fastq'"
if [[ -s "02_trim/${dataset}.2.fastq.gz" ]] ; then
  cmd="$cmd in2='02_trim/${dataset}.2.fastq.gz' \
    out2='03_norm/${dataset}.2.fastq'"
fi
cmd="$cmd target=40 min=2 threads=3 prefilter=t -Xmx${RAM}g"
$cmd
gzip -v 03_norm/${dataset}.[12].fastq

# Launch next step
qsub "$pkg/00_launcher.pbs" -N "GD04-$dataset" \
  -v "PKG=$pkg,TARGET=$target,DATASET=$dataset,STEP=04_asm" \
  -l nodes=1:ppn=12 -l mem="200g" -l walltime="72:00:00"
  -o "xx_log/${dataset}.02.txt" -j eo

