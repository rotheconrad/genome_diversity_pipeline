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

# Everything below will go to the project xx_summary/sizes.tsv file
if true ; then

  # Header
  echo -ne "dataset\tS01\tS02\tC02\tS03\tS04\tN05\tN06\tN07"
  for i in 02 03 04 05 06 07 ; do echo -ne "\tMEM$i\tWALL$i\tCPU$i" ; done
  echo

  # Body
  for dataset in $datasets ; do
    # Sizes, Cardinality, Number of files
    S01=$(size_01 "$dataset")
    S02=$(size_02 "$dataset")
    C02=$(cardinality "$dataset")
    S03=$(size_03 "$datset")
    S04=$(size_04 "$dataset")
    N05=$(ls 05_maxbin/${dataset}-*.d/*.fasta 2>/dev/null | wc -l)
    N06=$(ls 06_metabat/${dataset}-*.d/*.fa 2>/dev/null | wc -l)
    N07=$(cat "07_derep/${dataset}/genomospecies.tsv" 2>/dev/null | wc -l)
    echo -en "$dataset\t$S01\t$S02\t$C02\t$S03\t$S04\t$N05\t$N06\t$N07"
    # Running stats (CPU time, walltime, RAM)
    for i in 02 03 04 05 06 07 ; do
      echo -ne "\t$(run_stats "$dataset" "$i")"
    done
    echo
  done
fi > xx_summary/sizes.tsv 

echo "File updated: xx_summary/sizes.tsv"

