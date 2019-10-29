#!/bin/bash

target=$1
dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target || ! -n $dataset ]] ; then
  echo "Usage: $0 target_folder dataset"
  exit 0
fi

. "$pkg/00_env.bash"
cd "$target"
dir="07_derep/$dataset"
miga new -P "$dir" -t genomes -m aai_p=diamond,ani_p=fastani
miga add -P "$dir" -t popgenome -m run_mytaxa_scan=false -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix maxbin_ -v 05_maxbin/"$dataset"-*/*.fasta
miga add -P "$dir" -t popgenome -m run_mytaxa_scan=false -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix metabat_ -v 06_metabat/"$dataset"-*/*.fa
$HOME/shared3/miga-conf/start-node "$dir"

