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

dir="08_anir/$dataset"
mkdir -p "$dir"
for genome in cat 07_derep/${dataset}/*.LargeContigs.fna ; do
  name=$(basename "$genome" .LargeContigs.fna)
  reads="single"
  reads_file="02_trim/${dataset}.single.fa"
  if [[ ! -e "$reads_file" ]] ; then
    reads="interleaved"
    reads_file="02_trim/${dataset}.coupled.fa"
  fi

  # Run ANIr at different thresholds. Note that only the first
  # will run bowtie, all subsequent calls will simply reread the
  # SAM file produced by the first one.
  for identity in 90 95 97.5 ; do
    anir.rb -g "$genome" -r "$reads_file" --r-type "$reads" \
      -m "$dir/${genome}.sam" -t 12 -a fix -i "$identity" \
      -L "$dir/${genome}.identity.txt" \
      > "$dir/${genome}.anir-${identity}.txt"
  done

  # Compress to BAM and sort it
  samtools view -b "$dir/${genome}.sam" -@ 12 \
    | samtools sort -@ 12 -o "$dir/${genome}.bam" -
done

# Launch next step
"$pkg/00_launcher.bash" . "$dataset" 09
