#!/bin/sh

# Grep data from all lammps*.out.
# Sort by timestep, the first value of the thermo_style command

[ "$1" = "-h" ] || [ $# -eq 0 ] &&
echo "Returns data from lammps.out files in the order given as arguments.
Accounts for restarts and new runs, so the timesteps will always be consecutive.
Syntax: $(basename $0) [files]" && exit 1

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

# Assumes Step is the first value of thermo_style.
# If log is not the first file, and rerun was started with a 
# data file, the Step number starts at 0, when we want it to start at the
# next timestep. So append to data.tmp sequentially.

# find the header from first output
# need the tail because minimisation
# also prints out thermo data but with fewer columns
grep "Step" "$1" | tail -1 | $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > header.tmp

[ -f data.tmp ] && rm data.tmp
for f in "$@"
do
  # extract from "Setting up Verlet run..." onwards because 
  # minimisation is finished by that point
  lines="$(sed -n '/Setting up Verlet run .../,/Loop/p' "$f" |
    grep '^\s\+[0-9]\+\s\+-*[0-9]\+' |
    $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g')"
  # There might not be a line contatining 'Loop' line if the run hasn't finished,
  # in that case just read until the end of the file
  [[ $lines == "" ]] &&
    lines=$(sed -n '/Setting up Verlet run .../, $ p' "$f" | grep '^[0-9]' | $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g')
  # if more than one file
  [ -f data.tmp ] && [ $(cat data.tmp | wc -l) -gt 0 ] && (
  # take last two lines and find step increment
  laststep="$(tail -n 1 data.tmp | cut -d, -f1 | $sed 's/^\s*//;s/\s*$//')"
  # bc fails here if all on one line, not sure why
  increment="$(tail -2 data.tmp | cut -d, -f1 | tr '\n' ' ' | $sed 's/\s*$//;s/ / - /')"
  increment="$(echo $increment | bc | $sed 's/^-//')"
  # If run was from a restart, sometimes there will be a few timesteps
  # duplicated, so instead of taking the first timestep of the new run, take
  # the last step of the old run and use that as the start of the new run,
  # then increment by the $increment * line number of the new data
  grep -q "Reading restart file" $f &&
    echo "$lines" | awk -F"," -v last="$laststep" -v inc="$increment" '$1=last+(inc*NR)' OFS="," >> data.tmp ||
    # adding on a new run that starts from 0, so add on last timestep
    echo "$lines" | awk -F"," -v last="$laststep" -v inc="$increment" '$1+=last+inc' OFS="," >> data.tmp 
  ) || printf "%s\n" "$lines" > data.tmp # write data from first file to data.tmp
done

# if SHAKE algorithm used, the output includes SHAKE statistics- so remove these by checking
# the length of each line wrt to first line
numcols=$(cat header.tmp | tr ',' '\n' | wc -l)
cat header.tmp data.tmp | awk -F"," -v cols=$numcols 'NF==cols'
rm header.tmp data.tmp
