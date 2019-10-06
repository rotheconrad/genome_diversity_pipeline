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
    dir="06_metabat/${dataset}-${asm}.d"
    out="$dir/${dataset}-${asm}"
    mkdir -p "$dir"
    #[[ -d "$dir" ]] && continue
    # sam -> bam
    if [[ ! -s "${out}.bam" ]] ; then
      samtools view -S -b \
        "05_maxbin/${dataset}-${asm}.d/${dataset}-${asm}.sam0" \
        | samtools sort -l 9 - "$out" \
        && rm "05_maxbin/${dataset}-${asm}.d/${dataset}-${asm}.sam0"
    fi
    # depths
    "$HOME/shared3/apps/metabat/jgi_summarize_bam_contig_depths" \
      --outputDepth "${out}.abd" "${out}.bam"
    # binning
    "$HOME/shared3/apps/metabat/metabat2" \
      -i "04_asm/${dataset}-${asm}.LargeContigs.fna" \
      -a "${out}.abd" \
      -o "$out" \
      -t 4
  done
else
  for i in 04_asm/*-trim.LargeContigs.fna ; do
    d="$(basename "$i" -trim.LargeContigs.fna)"
    [[ -d "06_metabat/${d}-norm.d" ]] && continue
    echo -e "\033[31m==[ $d ]\e[0m"
    "$scr" . "$d"
  done
fi

