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

GDIV_PKG="$pkg" . "$pkg/00_env.bash"
cd "$target"
cmd="bbduk.sh in='01_reads/${dataset}.1.fastq.gz' \
  out='02_trim/${dataset}.1.fastq'"
if [[ -s "01_reads/${dataset}.2.fastq.gz" ]] ; then
  cmd="$cmd in2='01_reads/${dataset}.2.fastq.gz' \
    out2='02_trim/${dataset}.2.fastq'"
fi
cmd="$cmd ktrim=l qtrim=w,3 trimq=17 minlength=50 ref='$pkg/adapters.fa' tbo \
  tossjunk=t cardinalityout=t -Xmx50g"
echo "RUNNING BBDuk:"
echo "$cmd"
$cmd 2> 02_trim/${dataset}.log
gzip -v 02_trim/${dataset}.[12].fastq
card=$(grep 'Unique 31-mers out:' "02_trim/${dataset}.log" | cut -f 2)
let ram=1+$card*12/1000000000

# Launch next step
qsub "$pkg/00_launcher.pbs" -N "GD03-$dataset" \
  -v "PKG=$pkg,TARGET=$target,DATASET=$dataset,STEP=03_norm,RAM=$ram" \
  -l nodes=1:ppn=5 -l mem="$(($ram+10))g" -l walltime="24:00:00"
  -o "xx_log/${dataset}.02.txt" -j eo

