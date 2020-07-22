#!/bin/sh

# Grep data from all lammps*.out.
# Sort by timestep, the first value of the thermo_style command

[ "$1" = "-h" ] || [ $# -eq 0 ] &&
echo "Returns data from lammps.out files in the order given as arguments.
Accounts for restarts and new runs, so the timesteps will always be consecutive.
Syntax: $(basename $0) [files]" && exit 1

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

# data- assumes Step is the first value of thermo_style.
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
  lines="$(sed -n '/Setting up Verlet run .../,/Loop/p' $f |
    grep -v 'Step\|Loop' |
    grep '^\s\+[0-9]\+\s\+-*[0-9]\+' |
    $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g')"
  # if more than one file
  [ -f data.tmp ] && [ $(cat data.tmp | wc -l) -gt 0 ] && (
  # take last two lines and find step increment
  laststep="$(tail -n 1 data.tmp | cut -d, -f1 | $sed 's/^\s*//;s/\s*$//')"
  # bc fails here if all on one line, not sure why
  increment="$(tail -2 data.tmp | cut -d, -f1 | tr '\n' ' ' | $sed 's/\s*$//;s/ / - /')"
  increment="$(echo $increment | bc | $sed 's/^-//')"
  # If run was from a restart, timesteps continue from last run.
  # check first line of new data- first timestep will be the same as the last
  # of the previous run
  first_of_new_run="$(echo $lines | head -1 | cut -d, -f1 | $sed 's/^\s*//;s/\s*$//')" 
  [ "$first_of_new_run" = "$laststep" ] &&
    # continue from restart
    # density of 0th step is the same as the end of the last run, so drop the first line 
    echo "$lines" | tail -n +2 | awk -F"," -v inc="$increment" '$1+=inc' OFS="," >> data.tmp ||
    # adding on a new run that starts from 0, so add on last timestep
    echo "$lines" | awk -F"," -v last="$laststep" -v inc="$increment" '$1+=last+inc' OFS="," >> data.tmp 
  ) || printf "%s\n" "$lines" > data.tmp
done

cat header.tmp data.tmp 
rm header.tmp data.tmp
