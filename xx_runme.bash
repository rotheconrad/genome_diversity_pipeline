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

. "$pkg/00_env.bash"
. "$pkg/00_sizes.bash"
"$pkg/00_build.bash" "$target"
cd "$target"

# Use all datasets with logs if no datasets are specified
if [[ ! -n $datasets ]] ; then
  datasets=$(for i in $(ls xx_log) ; do basename $i \
    | perl -pe 's/\..*//' ; done | sort | uniq)
fi

for dataset in $datasets ; do
  echo $dataset

  # Determine the next step
  next_step=$(next_step "$dataset")

  # And launch it
  case "$next_step" in
    XX) echo "Workflow complete for $dataset" ;;
    01) echo "Empty dataset $dataset, run 01_add.bash first" ;;
    *)  "$pkg/00_launcher.bash" . "$dataset" "$next_step" ;;
  esac
done
