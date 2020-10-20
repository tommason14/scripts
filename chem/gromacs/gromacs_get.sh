#!/usr/bin/env bash

[[ $# -ne 2 || "$1" == "-h" ]] &&
  echo "Extract property from edr file using 'gmx energy'." &&
  echo "Syntax: $(basename $0) file.edr property" &&
  echo "i.e. $(basename $0) npt.edr Total-Energy" &&
  exit 1

file="$1"
property="$2"
[[ "$property" == *"-"* ]] &&
  ext="$(echo "$property" | sed 's/-/_/g' | tr '[A-Z]' '[a-z]').xvg" ||
  ext="$(echo "$property" | tr '[A-Z]' '[a-z]').xvg"
output="${file%.edr}_$ext"

command -v gmx && gmx="gmx" || gmx="gmx_mpi"

echo "$property" | "$gmx" -f "$file" -o "$output"
