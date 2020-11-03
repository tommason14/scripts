#!/usr/bin/env bash

[[ $# -eq 0 || "$1" == "-h" ]] &&
  echo "Plot 'gmx energy' output with gnuplot." && 
  echo "Syntax: $(basename $0) file.xvg" &&
  exit 1

[[ $USER =~ (tommason|tmas0023) ]] && 
  cmd="plot '-' using 1:2 with lines notitle" ||
  cmd="set terminal x11; plot '-' using 1:2 with lines notitle"

tail -n +24 $1 | gnuplot --persist -e "$cmd"
