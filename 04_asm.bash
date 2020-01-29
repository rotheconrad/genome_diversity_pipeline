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
# fq -> fa
for i in 0[23]_*/"$dataset".[12].fastq.gz ; do
  gzip -c -d "$i" | awk -f "$(which FastQ.toFastA.awk)" > "${i}.fa_tmp"
done
# paired
for i in 0[23]_*/"$dataset".2.fastq.gz.fa_tmp ; do
  [[ -s "$i" ]] || continue
  d=$(dirname "$i")
  FastA.interpose.pl "$d/${dataset}.coupled.fa" \
    "$d/${dataset}".[12].fastq.gz.fa_tmp
  rm "$d/${dataset}".[12].fastq.gz.fa_tmp
done
# single
for i in 0[23]_*/"$dataset".1.fastq.gz.fa_tmp ; do
  [[ -s "$i" ]] || continue
  d=$(dirname "$i")
  mv "$i" "$d/${dataset}.single.fa"
done
# assemble
for i in 0[23]*_*/"$dataset".*.fa ; do
  var=$(echo "$(dirname "$i")" | perl -pe 's/^\d+_//')
  base="04_asm/${dataset}-${var}"
  dir="${base}.d"
  rd="r"
  [[ "$i" == *.single.fa.gz ]] && rd="l"
  idba_ud -o "$dir" -$rd "$i" --num_threads 24 --maxk 120
  # link result
  if [[ -s "$dir/scaffold.fa" ]] ; then
    ln "$dir/scaffold.fa" "${base}.AllContigs.fna"
  else
    ln "$dir/contig.fa" "${base}.AllContigs.fna"
  fi
  # filter by length
  FastA.length.pl "${base}.AllContigs.fna" | awk '$2 >= 1000 { print $1 }' \
    | FastA.filter.pl /dev/stdin "${base}.AllContigs.fna" \
    > "${base}.LargeContigs.fna"
  # cleanup
  rm "${base}.AllContigs.fna"
  rm -r "$dir"
done

# Launch next step
qsub "$pkg/00_launcher.pbs" -N "GD05-$dataset" \
  -v "PKG=$pkg,TARGET=$target,DATASET=$dataset,STEP=05_maxbin" \
  -l nodes=1:ppn=12 -l mem="240g" -l walltime="90:00:00" \
  -o "xx_log/${dataset}.05.txt" -j oe


