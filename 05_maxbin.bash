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

