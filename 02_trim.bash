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
cd "$target"

cmd="bbduk.sh in='01_reads/${dataset}.1.fastq.gz' \
  out='02_trim/${dataset}.1.fastq'"
if [[ -s "01_reads/${dataset}.2.fastq.gz" ]] ; then
  cmd="$cmd in2='01_reads/${dataset}.2.fastq.gz' \
    out2='02_trim/${dataset}.2.fastq'"
fi
cmd="$cmd qtrim=w,3 trimq=17 minlength=70 ref='$pkg/adapters.fa' tbo \
  tossjunk=t cardinalityout=t -Xmx50g"
echo "RUNNING BBDuk:"
echo "$cmd"
$cmd 2> 02_trim/${dataset}.log
gzip -v 02_trim/${dataset}.[12].fastq

# Launch next step
"$pkg/00_launcher.bash" . "$dataset" 03
