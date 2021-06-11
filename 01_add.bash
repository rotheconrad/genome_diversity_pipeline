#!/bin/bash

# pass in directories and files input from the command line
target=$1
dataset=$2
file1=$3
file2=$4
# get the link to this script
scr=$(readlink -f "$0" 2>/dev/null)
# get the directory this script is stored in
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset || ! -n $file1 ]] ; then
  echo "Usage: $0 target_folder dataset file1 [file2]
  target_folder    Path to the target folder containing the data
  dataset          Name of the dataset to create
  file1            Path to the FWD (or unpaired) reads, in .fastq.gz format
  file2 (optional) Path to the REV reads, in .fastq.gz format
  "
  exit 0
fi

if [[ "$(echo "$dataset" | perl -pe 's/[^A-Z0-9_]//i')" != "$dataset" ]] ; then
  echo "Please use only alpha-numerics and underscores in the dataset names" >&2
  exit 1
fi

. "$pkg/00_env.bash"
"$pkg/00_build.bash" "$target"

cp "$file1" "$target/01_reads/${dataset}.1.fastq.gz"
[[ -n $file2 ]] && cp "$file2" "$target/01_reads/${dataset}.2.fastq.gz"

# Launch next step
"$pkg/00_launcher.bash" "$target" "$dataset" 02

