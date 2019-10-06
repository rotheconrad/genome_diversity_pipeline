#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder [dataset]"
  exit 0
fi

GDIV_PKG="$pkg" . "$pkg/00_env.bash"
cd "$target"
if [[ -n $dataset ]] ; then
  cmd="bbduk.sh in='01_reads/${dataset}.1.fastq.gz' out='02_trim/${dataset}.1.fastq'"
  if [[ -s "01_reads/${dataset}.2.fastq.gz" ]] ; then
    cmd="$cmd in2='01_reads/${dataset}.2.fastq.gz' out2='02_trim/${dataset}.2.fastq'"
  fi
  cmd="$cmd ktrim=l qtrim=rl minlength=50 ref='$pkg/adapters.fa' tbo tpe"
  $cmd
  gzip -v 02_trim/${dataset}.[12].fastq
else
  for i in 01_reads/*.1.fastq.gz ; do
    d=$(basename "$i" .1.fastq.gz)
    [[ -s "02_trim/${d}.1.fastq.gz" ]] && continue
    echo -e "\033[31m==[ $d ]\e[0m"
    "$scr" . "$d"
  done
fi

