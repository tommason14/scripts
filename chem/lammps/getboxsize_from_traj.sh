#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Finds box size by averaging last 30 timesteps printed into the lammps trajetory output.
Syntax: $(basename $0) traj.lmp" && exit 1

sed -n '/ITEM: BOX/{n;p}' $1 | awk '{sum+=($2 - $1)} END {print sum/NR}'
