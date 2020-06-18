#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] && echo "Syntax: $(basename $0)
lammps_dump" && exit 1

# all x components to average over, then find xhi - xlo and average
sed -n '/ITEM: BOX/{n;p}' $1 | awk '{sum+=($2 - $1)} END {print sum/NR}'
