#!/usr/bin/env bash

[[ $# -ne 2 || $1 == '-h' ]] && echo "Syntax: $(basename $0) lammps_datafile molecule_name" && exit 1
[[ ! -f $1 ]] && echo "$1 not found" && exit 1

sed -n '/Atoms/,/Bonds/p' $1 |
grep "$2" |
grep -v 'DP$' |
awk '{print $3}' |
sort -nu |
tr '\n' ' ' |
sed 's/ $//'
