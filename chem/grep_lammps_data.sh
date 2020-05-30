#!/bin/sh

# Grep data from all log.lammps* files.
# Sort by timestep, the first value of the thermo_style command
# Print to first argument, or data.csv otherwise.

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

output=${1-data.csv}
echo "Writing lammps data to $output"

# data- assumes Step is the first value of thermo_style
ls log.lammps* | 
  xargs $sed -n '/Step/,/Loop/p' |
  grep -v 'Step\|Loop' |
  sort -nk 1 |
  grep '^\s\+[0-9]\+\s\+-*[0-9]\+' |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > data.tmp

# find the header from first output
ls log.lammps* |
  head -1 |
  xargs grep Step |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > header.tmp

cat header.tmp data.tmp > $output

rm header.tmp data.tmp
