#!/bin/sh

# Grep data from all log.lammps* files or lammps*.out.
# Sort by timestep, the first value of the thermo_style command
# Print to first argument, or data.csv otherwise.

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

[ $# -eq 0 ] && echo "\
Pass in log or lammps to take data from log.lammps* or lammps*.out
Syntax: $(basename $0) [log|lammps] <csvfile>" && exit 1

choice=$1
output=${2-data.csv}
echo "Writing lammps data to $output"

# data- assumes Step is the first value of thermo_style
ls "$choice"* |
  xargs $sed -n '/Step/,/Loop/p' |
  grep -v 'Step\|Loop' |
  sort -nk 1 |
  grep '^\s\+[0-9]\+\s\+-*[0-9]\+' |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > data.tmp

# find the header from first output
ls "$choice"* |
  head -1 |
  xargs grep Step |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > header.tmp

cat header.tmp data.tmp > $output

rm header.tmp data.tmp
