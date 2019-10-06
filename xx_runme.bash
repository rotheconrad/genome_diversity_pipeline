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
"$pkg/01_build.bash" "$target"
[[ -s "$target/02_trim/${dataset}.1.fastq.gz" ]] \
  || "$pkg/02_trim.bash" "$target" "$dataset"

