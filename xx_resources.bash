#!/bin/bash

target=$1
shift
datasets=$@
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder [datasets]
  target_folder    Path to the target folder containing the data
  datasets (opt)   Name of the dataset(s), by default: all the datasets
  "
  exit 0
fi

. "$pkg/00_sizes.bash"
"$pkg/00_build.bash" "$target"
cd "$target"

# Use all datasets with logs if no datasets are specified
if [[ ! -n $datasets ]] ; then
  datasets=$(for i in $(ls xx_log) ; do basename $i \
    | perl -pe 's/\..*//' ; done | sort | uniq)
fi

# Generate tab-delimited summary
if true; then
  echo -e \
    "dataset\tstep\tram_req\tram_use\twall_req\twall_use\tcpu_req\tcpu_use"
  for dataset in $datasets ; do
    for i in 02 03 04 05 06 07 ; do
      echo -e "$dataset\t$i\t$(
        "$pkg/scripts/xx_resources.rb" "xx_log/${dataset}.${i}.txt")"
    done
  done
fi > xx_summary/resources.tsv
echo "File updated: xx_summary/resources.tsv"
"$pkg/scripts/xx_resources.R" xx_summary/resources.tsv xx_summary/resources.pdf
echo "File updated: xx_summary/resources.pdf"

