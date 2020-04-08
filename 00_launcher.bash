#!/bin/bash

function launch_step_02 {
  local dataset=$1

  # Determine time by input size
  SIZE=$(file_size "01_reads/${dataset}.1.fastq.gz")
  if [[ -n $file2 ]] ; then
    SIZE2=$(file_size "01_reads/${dataset}.2.fastq.gz")
    let SIZE=$SIZE+$SIZE2
  fi
  let SIZE_G=$SIZE/1000000000
  let TIME_H=2+$SIZE_G

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD02-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=02_trim" \
    -l nodes=1:ppn=1 -l mem=50g -l "walltime=${TIME_H}:00:00" \
    -o "xx_log/${dataset}.02.txt" -j oe
}

function launch_step_03 {
  local dataset=$1

  # Determine time and RAM by cardinality
  C02=$(cardinality "$dataset")
  RAM_G=$((1+$C02*12/1000000000))
  TIME_H=$((12+$C02*2/1000000000))

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD03-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=03_norm,RAM=$RAM_G" \
    -l nodes=1:ppn=5 -l mem="$(($RAM_G+10))g" -l walltime="${TIME_H}:00:00" \
    -o "xx_log/${dataset}.03.txt" -j oe
}

function launch_step_035 {
  local dataset=$1

  # Determine time and RAM by cardinality
  C02=$(cardinality "$dataset")
  RAM_G=$((1+$C02*12/1000000000))
  TIME_H=$((12+$C02*2/1000000000))

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD035-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=035_norm,RAM=$RAM_G" \
    -l nodes=1:ppn=5 -l mem="$(($RAM_G+10))g" -l walltime="${TIME_H}:00:00" \
    -o "xx_log/${dataset}.035.txt" -j oe
}

function launch_step_04 {
  local dataset=$1

  # Determine time and RAM by cardinality
  C02=$(cardinality "$dataset")
  TIME_H=$((12+$C02*13/1000000000))
  RAM_G=$((100+$C02*45/1000000000))
  [[ $RAM_G -gt 500 ]] && RAM_G+500 # <- That's the maximum we have

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD04-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=04_asm" \
    -l nodes=1:ppn=24 -l mem="${RAM_G}g" -l walltime="${TIME_H}:00:00" \
    -o "xx_log/${dataset}.04.txt" -j oe
}

function launch_step_05 {
  local dataset=$1

  # Determine time by read size
  S02=$(file_size "02_trim/${dataset}.1.fastq.gz")
  TIME_H=$((6+$S02/1000000000))

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD05-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=05_maxbin" \
    -l nodes=1:ppn=12 -l mem="24g" -l walltime="${TIME_H}:00:00" \
    -o "xx_log/${dataset}.05.txt" -j oe
}

function launch_step_06 {
  local dataset=$1

  # Determine time by read size
  S02=$(file_size "02_trim/${dataset}.1.fastq.gz")
  TIME_H=$((6+$S02/1000000000))

  # Launch
  qsub "$pkg/00_launcher.pbs" -N "GD06-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=06_metabat" \
    -l nodes=1:ppn=12 -l mem="12g" -l walltime="${TIME_H}:00:00" \
    -o "xx_log/${dataset}.06.txt" -j oe
}

function launch_step_07 {
  local dataset=$1

  qsub "$pkg/00_launcher.pbs" -N "GD07-$dataset" \
    -v "PKG=$pkg,TARGET=$PWD,DATASET=$dataset,STEP=07_derep" \
    -l nodes=1:ppn=12 -l mem="12g" -l walltime="24:00:00" \
    -o "xx_log/${dataset}.07.txt" -j oe
}

function launch_step_08 {
  echo "STEP 08: NOT YET IMPLEMENTED" >&2
}

# Source code
scr=$(readlink -f "$0" 2>/dev/null)
export pkg=$(dirname "$scr")
. "$pkg/00_sizes.bash"

# Input variables
target="$1"
dataset="$2"
step="$3"

# Engage!
cd "$target"
"launch_step_$step" "$dataset"

