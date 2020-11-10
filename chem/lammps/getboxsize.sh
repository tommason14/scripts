#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Finds box size by averaging last 30 timesteps printed, cube-rooting the volume.
Syntax: $(basename $0) lammps.out" && exit 1

# Find volume column
col=$(grep Step $1 | tail -1 | tr ' ' '\n' | nl | grep Volume | awk '{print $1}')

grep_lammps_data.sh $1 | tail -30 | awk -F"," -v col=$col '{vol+=$col} END {print (vol/NR)**(1/3)}'
