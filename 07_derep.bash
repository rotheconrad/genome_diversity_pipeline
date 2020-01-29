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
dir="07_derep/$dataset"
miga new -P "$dir" -t genomes
miga add -P "$dir" -t popgenome -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix maxbin_ -v 05_maxbin/"$dataset"-*/*.fasta
miga add -P "$dir" -t popgenome -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix metabat_ -v 06_metabat/"$dataset"-*/*.fa
miga derep_wf -o "$dir" --fast -j 12 -t 1 -v \
  --daemon "$HOME/shared3/miga-conf/daemon_bash.json"


