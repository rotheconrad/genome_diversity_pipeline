#!/bin/bash

target=$1
shift
datasets=$@
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder [datasets]
  target_folder    Path to the target folder containing the data
  datasets (opt)   Name of the dataset(s), by default: all the datasets
  "
  exit 0
fi

. "$pkg/00_env.bash"
"$pkg/00_build.bash" "$target"
cd "$target"

# Use all datasets with logs if no datasets are specified
if [[ ! -n $datasets ]] ; then
  datasets=$(for i in $(ls xx_log) ; do basename $i \
    | perl -pe 's/\..*//' ; done | sort | uniq)
fi

for dataset in $datasets ; do
  # Determine the next step
  if [[ ! -s "01_reads/${dataset}.1.fastq.gz" ]] ; then
    next_step=01
  elif [[ ! -s "02_trim/${dataset}.1.fastq.gz" ]] ; then
    next_step=02
  elif [[ ! -s "03_norm/${dataset}.1.fastq.gz" ]] ; then
    next_step=03
  elif [[ ! -s "04_asm/${dataset}-norm.LargeContigs.fna" ]] ; then
    next_step=04
  elif [[ ! -s "05_maxbin/${dataset}-norm.d/${dataset}-norm.summary" ]] ; then
    next_step=05
  elif [[ ! -s "06_metabat/${dataset}-norm.d/${dataset}-norm.1.fa" ]] ; then
    next_step=06
  elif [[ ! -s "07_derep/${dataset}/genomospecies.tsv" ]] ; then
    next_step=07
  elif [[ ! -s "08_anir/${dataset}/anir-95.tsv" ]] ; then
    next_step=08
  else
    next_step=XX
  fi

  # And launch it
  case "$next_step" in
    XX) echo "Workflow complete for $dataset" ;;
    01) echo "Empty dataset $dataset, run 01_add.bash first" ;;
    *)  "$pkg/00_launcher.bash" . "$dataset" "$next_step" ;;
  esac
done
