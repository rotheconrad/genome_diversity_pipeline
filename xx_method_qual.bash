#!/bin/bash

target=$1
shift
scr=$(readlink -f "$0")
pkg=$(dirname "$scr")

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder
  target_folder    Path to the target folder containing the data
  "
  exit 0
fi

. "$pkg/00_sizes.bash"
"$pkg/00_build.bash" "$target"
cd "$target"

"$pkg/scripts/07_02_gsp_qual.R" xx_summary/method_qual.pdf \
  07_derep/*/method_qual.tsv
echo "File updated: xx_summary/method_qual.pdf"

