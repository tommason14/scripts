#!/bin/sh

[ $# -eq 0 ] || [ "$1" = "-h" ] && 
  echo "Unwraps trajectory." &&
  echo "Syntax: $(basename $0) file.xtc" && 
  echo "Assumes that there is a tpr file with the same name in the same directory." && exit 1

file="$1"
tpr="${file%.xtc}.tpr"
output="${file%.xtc}_unwrapped.xtc"

command -v gmx && gmx="gmx" || gmx="gmx_mpi"

echo "System" | "$gmx" trjconv -f "$file" -s "$tpr" -pbc mol -o "$output"
