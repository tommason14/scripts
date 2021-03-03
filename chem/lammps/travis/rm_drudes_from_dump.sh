#!/usr/bin/env bash

[[ $# -eq 0 ]] &&
echo "Remove drude particles from a lammps trajectory file
and adjusts atom numbers accordingly.
Produces a file named no_drude_filename.

Syntax: $(basename $0) filename" && exit 1

# print first timestep by taking lines between the first two 
# occurrences of TIMESTEP
num_atoms_no_drudes=$(
  awk -v "f=0" '/TIMESTEP/{f=1;++d}d!=2{print}f&&d==2{print;exit}' $1 |
  grep -v 'ITEM' | grep '^\s*[A-Z]' | grep -v '^\s*D' | wc -l)  

orig_atoms=$(sed '/BOX BOUNDS/q' $1 | tail -n 2 | head -1)

sed "s|^\s*$orig_atoms.*|$num_atoms_no_drudes|" $1 | grep -v "^\s*D" > no_drude_"$1"
