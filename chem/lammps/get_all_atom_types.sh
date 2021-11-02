#!/usr/bin/env bash
[[ $# -ne 1 || ! -f $1 ]] && 
echo "Print out atom types for each atom in a LAMMPS datafile. 
Reads atom types from the Masses section.
Syntax: $(basename $0) datafile" && exit 1

fname="$1"
declare -A types

while read num element
do
  types[$num]=$element
done < <(sed -n '/Masses/,/^[A-Z]/p' "$fname" | grep '^\s*[0-9]' | awk '{print $1,$NF}')

# one long sed command to replace numbers with types:
# sed 's/^1$/C3/;s/^2$/HA/;'
addtypes="sed "
for t in "${!types[@]}"
do
  addtypes+="s/^$t$/${types[$t]}/;"
done

sed -n '/Atoms/,/^[A-Z]/p' "$fname" | grep '^\s*[0-9]' | sort -k1n | awk '{print $3}' | $addtypes | tr '\n' ' ' | sed 's/ $//'
