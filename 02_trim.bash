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
cmd="bbduk.sh in='01_reads/${dataset}.1.fastq.gz' out='02_trim/${dataset}.1.fastq'"
if [[ -s "01_reads/${dataset}.2.fastq.gz" ]] ; then
  cmd="$cmd in2='01_reads/${dataset}.2.fastq.gz' out2='02_trim/${dataset}.2.fastq'"
fi
cmd="$cmd ktrim=l qtrim=w,3 trimq=17 minlength=50 ref='$pkg/adapters.fa' tbo \
  tossjunk=t cardinalityout=t -Xmx50g"
$cmd
gzip -v 02_trim/${dataset}.[12].fastq
#qsub "$pkg/00_launcher.pbs" -N "GD03-$dataset" \
#  -v "PKG=$pkg,TARGET=$target,DATASET=$dataset,STEP=03_norm" \
#  -l nodes=1:ppn=1 -l mem=50g -l walltime=$TIME_H:00:00

