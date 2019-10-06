#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder [dataset]"
  exit 0
fi

. "$pkg/00_env.bash"
cd "$target"
if [[ -n $dataset ]] ; then
  for asm in trim norm ; do
    dir="05_maxbin/${dataset}-${asm}.d"
    out="$dir/${dataset}-${asm}"
    mkdir -p "$dir"
    #[[ -d "$dir" ]] && continue
    "$HOME/shared3/apps/MaxBin-2.2.7/run_MaxBin.pl" \
      -contig "04_asm/${dataset}-${asm}.LargeContigs.fna" \
      -reads  "02_trim/${dataset}.coupled.fa" \
      -out    "$out" \
      -thread 4 \
      -preserve_intermediate
  done
else
  for i in 04_asm/*-trim.LargeContigs.fna ; do
    d="$(basename "$i" -trim.LargeContigs.fna)"
    [[ -d "05_maxbin/${d}-norm.d" ]] && continue
    echo -e "\033[31m==[ $d ]\e[0m"
    "$scr" . "$d"
  done
fi

