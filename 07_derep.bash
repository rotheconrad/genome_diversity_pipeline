#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder [dataset]"
  exit 0
fi

. "$pkg/00_env.bash"
cd "$target"
if [[ -n $dataset ]] ; then
  dir="07_derep/$dataset"
  miga new -P "$dir" -t genomes -m aai_p=diamond,ani_p=fastani
  miga add -P "$dir" -t popgenome -m run_mytaxa_scan=false -i assembly \
    -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
    --prefix maxbin_ -v 05_maxbin/"$dataset"-*/*.fasta
  miga add -P "$dir" -t popgenome -m run_mytaxa_scan=false -i assembly \
    -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
    --prefix metabat_ -v 06_metabat/"$dataset"-*/*.fa
  $HOME/shared3/miga-conf/start-node "$dir"
else
  for i in 04_asm/*-trim.LargeContigs.fna ; do
    d="$(basename "$i" -trim.LargeContigs.fna)"
    [[ -d "07_derep/$d" ]] && continue
    echo -e "\033[31m==[ $d ]\e[0m"
    "$scr" . "$d"
  done
fi

