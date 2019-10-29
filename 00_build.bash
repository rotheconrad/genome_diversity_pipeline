#!/bin/bash

target=$1

if [[ ! -n $target ]] ; then
  echo "Usage: $0 target_folder
  target_folder: Path to the target folder to create
  "
  exit 0
fi

mkdir -p "$target"
cd "$target"
mkdir -p 01_reads 02_trim 03_norm 04_asm 05_maxbin 06_metabat 07_derep xx_log

