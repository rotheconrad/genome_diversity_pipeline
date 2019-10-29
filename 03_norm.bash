#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset ]] ; then
  echo "Usage: $0 target_folder dataset"
  exit 0
fi

. "$pkg/00_env.bash"
cd "$target"
cmd="bbnorm.sh in='02_trim/${dataset}.1.fastq.gz' out='03_norm/${dataset}.1.fastq'"
if [[ -s "02_trim/${dataset}.2.fastq.gz" ]] ; then
  cmd="$cmd in2='02_trim/${dataset}.2.fastq.gz' out2='03_norm/${dataset}.2.fastq'"
fi
cmd="$cmd target=40 min=2"
$cmd
gzip -v 03_norm/${dataset}.[12].fastq

