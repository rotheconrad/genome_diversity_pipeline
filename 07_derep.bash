#!/bin/bash

# Change `$dataset` variable to `$gd_dataset` to avoid conflicts with
# MiGA's `$DATASET`

gd_target=$1
gd_dataset=$2
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $gd_target || ! -n $gd_dataset ]] ; then
  echo "Usage: $0 target_folder dataset
  target_folder Path to the target folder containing the data
  dataset       Name of the dataset to process
  "
  exit 0
fi

# We don't need the genome_diversity_pipeline environment for this one
#. "$pkg/00_env.bash"
cd "$gd_target"

dir="07_derep/$gd_dataset"

miga new -P "$dir" -t genomes
miga add -P "$dir" -t popgenome -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix maxbin_ -v 05_maxbin/"$gd_dataset"-*/*.fasta
miga add -P "$dir" -t popgenome -i assembly \
  -R '^(?:.*\/)?(.+?)(?i:\.f[nastq]+)?$' \
  --prefix metabat_ -v 06_metabat/"$gd_dataset"-*/*.fa

miga derep_wf -o "$dir" --fast -j 12 -t 1 -v 

# load miga environment to run this ruby script.
eval "$(miga env)"
"$pkg/scripts/07_01_gsp_qual.rb" "$dir" > "$dir/method_qual.tsv"

# Launch next step
"$pkg/00_launcher.bash" . "$gd_dataset" 08
