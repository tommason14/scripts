#!/bin/sh
cp $(dirname $0)/plot_density.R .
grep_lammps_data.sh lammps data.csv && Rscript plot_density.R && preview plot.png
