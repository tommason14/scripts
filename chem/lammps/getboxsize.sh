#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Finds box size by averaging last 30 timesteps printed, cube-rooting the volume.
Assumes volume is the second to last column of the thermo output.
Syntax: $(basename $0) lammps.out" && exit 1

grep '^\s\+[0-9]\+\s\+[0-9]\+' $1 | tail -30 | awk '{vol+=$(NF-1)} END {print (vol/NR)**(1/3)}'
