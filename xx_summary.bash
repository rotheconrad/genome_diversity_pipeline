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
  # Determine the next step
  next_step=$(next_step "$dataset")
  if [[ "$next_step" == "XX" ]] ; then
    next_step="Complete"
  else
    next_step="Next step: $next_step"
  fi

  echo -e "$(tput setaf 2)${dataset}:$(tput sgr 0) $next_step"
done
