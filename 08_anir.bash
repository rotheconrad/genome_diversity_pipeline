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

# Compile database
for genome in 07_derep/${dataset}/representatives/*.LargeContigs.fna ; do
  name=$(basename "$genome" .LargeContigs.fna)
  perl -pe 's/^>/>'$name':/' < "$genome"
done > 07_derep/${dataset}/representatives.fna

# Determine query settings
reads="single"
reads_file="02_trim/${dataset}.1.fastq.gz"
if [[ -e "02_trim/${dataset}.2.fastq.gz" ]] ; then
  reads="coupled"
  reads_file="${reads_file},02_trim/${dataset}.2.fastq.gz"
fi

# Run base ANIr to generate SAM file
anir.rb -g "07_derep/${dataset}/representatives.fna" \
  -r "$reads_file" --r-type "$reads" --r-format fastq \
  -m "$dir/map.sam" -t 12 -a fix -i 90

# Compress to BAM and sort it
samtools view -b "$dir/map.sam" -@ 12 \
  | samtools sort -@ 12 -o "$dir/map.bam" -
rm "$dir/map.sam"

# Run ANIr for each genome
for genome in 07_derep/${dataset}/representatives/*.LargeContigs.fna ; do
  # Run ANIr at different identity thresholds
  for identity in 97.5 95 90 ; do
    anir.rb -g "$genome" -m "$dir/map.bam" --m-format bam \
      -t 12 -a fix -i "$identity" -L "$dir/${name}.identity.txt" \
      --tab "$dir/${name}.anir-${identity}.tsv"
  done
done

for i in $dir/*.anir-95.txt ; do
  echo -e "$(basename "$i" .anir-95.txt)\t$(tail -n 1 "$i")"
done > $dir/anir-95.tsv

# Launch next step
"$pkg/00_launcher.bash" . "$dataset" 09
