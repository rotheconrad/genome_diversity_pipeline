#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset ]] ; then
  echo "Usage: $0 target_folder dataset [file(s)]
  target_folder    Path to the target folder containing the data
  dataset          Name of the dataset
  file(s) (opt)    Path to the reads if the dataset doesn't exist yet
  "
  exit 0
fi

. "$pkg/00_env.bash"
"$pkg/00_build.bash" "$target"
cd "$target"

if [[ ! -s "01_reads/${dataset}.1.fastq.gz" ]] ; then
  next_step=01
elif [[ ! -s "02_trim/${dataset}.1.fastq.gz" ]] ; then
  next_step=02
elif [[ ! -s "03_norm/${dataset}.1.fastq.gz" ]] ; then
  next_step=03
# Do no run 035, just a temporary test
#elif [[ ! -s "035_norm/${dataset}.1.fastq.gz" ]] ; then
#  next_step=035
elif [[ ! -s "04_asm/${dataset}-norm.LargeContigs.fna" ]] ; then
  next_step=04
elif [[ ! -s "05_maxbin/${dataset}-norm.d/${dataset}-norm.summary" ]] ; then
  next_step=05
elif [[ ! -s "06_metabat/${dataset}-norm.d/${dataset}-norm.1.fa" ]] ; then
  next_step=06
elif [[ ! -s "07_derep/${dataset}/genomospecies.tsv" ]] ; then
  next_step=07
else
  next_step=XX
fi

case "$next_step" in
  XX) echo "Workflow complete for $dataset" ;;
  01) "$pkg/01_add.bash" $@ ;;
  *)  "$pkg/00_launcher.bash" . "$dataset" "$next_step" ;;
esac

