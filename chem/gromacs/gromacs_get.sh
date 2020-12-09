#!/usr/bin/env bash

[[ $# -ne 2 || "$1" == "-h" ]] &&
  echo "Extract property from edr file using 'gmx energy'." &&
  echo "Syntax: $(basename $0) file.edr property" &&
  echo "i.e. $(basename $0) npt.edr Total-Energy" &&
  exit 1

file="$1"
property="$2"
ext="$(echo "$property" | sed 's/-/_/g' | tr '[A-Z]' '[a-z]').xvg"
output="${file%.edr}_$ext"

if command -v gmx; then
  gmx="gmx" 
elif module list 2>&1 | grep -Fq gromacs; then 
  gmx="gmx_mpi"
else 
  module load gromacs/5.1.4
  gmx="gmx_mpi"
fi

echo "$property" | "$gmx" energy -f "$file" -o "$output"
