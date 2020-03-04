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

# Header
echo -ne "dataset\tS01\tS02\tC02\tS03\tS035\tS04\tN05\tN06\tN07"
for i in 02 03 035 04 05 06 07 ; do echo -ne "\tCPU$i\tWALL$i\tMEM$i" ; done
echo

# Body
for dataset in $datasets ; do
  # Sizes, Cardinality, Number of files
  S01=$(file_size "01_reads/${dataset}.1.fastq.gz")
  [[ -e "01_reads/${dataset}.2.fastq.gz" ]] \
    && let S01=$S01+$(file_size "01_reads/${dataset}.2.fastq.gz")
  S02=$(file_size "02_trim/${dataset}.1.fastq.gz")
  [[ -e "02_trim/${dataset}.2.fastq.gz" ]] \
    && let S02=$S02+$(file_size "02_trim/${dataset}.2.fastq.gz")
  C02=$(cardinality "$dataset")
  S03=$(file_size "03_norm/${dataset}.1.fastq.gz")
  [[ -e "03_norm/${dataset}.2.fastq.gz" ]] \
    && let S03=$S03+$(file_size "03_norm/${dataset}.2.fastq.gz")
  S035=$(file_size "035_norm/${dataset}.1.fastq.gz")
  [[ -e "035_norm/${dataset}.2.fastq.gz" ]] \
    && let S035=$S035+$(file_size "035_norm/${dataset}.2.fastq.gz")
  S04=$(file_size "04_asm/${dataset}-norm.LargeContigs.fna")
  let S04=$S04+$(file_size "04_asm/${dataset}-trim.LargeContigs.fna")
  N05=$(ls 05_maxbin/${dataset}-*.d/*.fasta 2>/dev/null | wc -l)
  N06=$(ls 06_metabat/${dataset}-*.d/*.fa 2>/dev/null | wc -l)
  N07=$(cat "07_derep/${dataset}/genomospecies.tsv" 2>/dev/null | wc -l)
  echo -en "$dataset\t$S01\t$S02\t$C02\t$S03\t$S035\t$S04\t$N05\t$N06\t$N07"
  # Running stats (CPU time, walltime, RAM)
  for i in 02 03 035 04 05 06 07 ; do
    echo -ne "\t$(run_stats "$dataset" "$i")"
  done
  echo
done

