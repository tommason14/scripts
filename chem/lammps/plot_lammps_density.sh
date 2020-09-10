#!/bin/sh

[ ! -f "data.csv" ] &&
echo "This script reads a csv named data.csv in the current directory, but data.csv doesn't exist!
Run grep_lammps_data.sh on one or more lammps.out files and save to data.csv" && exit 1

cp $(dirname $0)/plot_density.R .
Rscript plot_density.R 2>/dev/null && preview plot.png
rm Rplots.pdf
