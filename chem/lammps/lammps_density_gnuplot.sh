#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Plots density against timesteps using gnuplot.
Syntax: $(basename $0) lammps.out" && exit 1

# Find step and density columns
step=$(grep Step $1 | tail -1 | tr ' ' '\n' | nl | grep Step | awk '{print $1}')
dens=$(grep Step $1 | tail -1 | tr ' ' '\n' | nl | grep Density | awk '{print $1}')

grep_lammps_data.sh $1 |
  sed '1d' |
  awk -F"," -v step=$step -v dens=$dens '{print $step,$dens}' |
  gnuplot --persist -e "set xlabel 'Timestep'; set ylabel 'Density'; plot '-' using 1:2 with lines notitle"
