#!/bin/bash

function cardinality {
  grep 'Unique 31-mers out:' "02_trim/${1}.log" | cut -f 2
}

function file_size {
  if [[ -e "$1" ]] ; then
    ls -pl "$1" | awk '{print $5}'
  else
    echo 0
  fi
}

function time_to_minutes {
  local t=$1
  echo "$t" | tr ':' ' ' | awk '{ print $1*60+$2+$3/60 }'
}

function run_stats {
  local dataset=$1
  local step=$2
  file="xx_log/${dataset}.${step}.txt"
  if [[ -s "$file" ]] ; then
    rsrc=$(grep '^Rsrc Used:' "$file" | head -n 1)
    cput=$(time_to_minutes "$(echo "$rsrc" | perl -pe 's/.*cput=(.+?),.*/$1/')")
    wallt=$(time_to_minutes "$(echo "$rsrc"| perl -pe 's/.*walltime=(.+)/$1/')")
    mem=$(echo "$rsrc" | perl -pe 's/.*mem=(.+?)kb,.*/$1/')
  fi
  echo -e "$cput\t$wallt\t$mem"
}

