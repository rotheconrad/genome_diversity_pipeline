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
for genome in 07_derep/${dataset}/representatives/*.LargeContigs.fna ; do
  name=$(basename "$genome" .LargeContigs.fna)
  reads="single"
  reads_file="02_trim/${dataset}.1.fastq.gz"
  if [[ -e "02_trim/${dataset}.2.fastq.gz" ]] ; then
    reads="coupled"
    reads_file="${reads_file},02_trim/${dataset}.2.fastq.gz"
  fi

  # Run ANIr at different thresholds. Note that only the first
  # will run bowtie, all subsequent calls will simply reread the
  # SAM file produced by the first one.
  for identity in 90 95 97.5 ; do
    anir.rb -g "$genome" -r "$reads_file" --r-type "$reads" --r-format fastq \
      -m "$dir/${name}.sam" -t 12 -a fix -i "$identity" \
      -L "$dir/${name}.identity.txt" \
      > "$dir/${name}.anir-${identity}.txt"
  done

  # Compress to BAM and sort it
  samtools view -b "$dir/${name}.sam" -@ 12 \
    | samtools sort -@ 12 -o "$dir/${name}.bam" -
done

for i in $dir/*.anir-95.txt ; do
  echo -e "$(basename "$i" .anir-95.txt)\t$(grep ANIr "$i")"
done > $dir/anir-95.tsv

# Launch next step
"$pkg/00_launcher.bash" . "$dataset" 09
