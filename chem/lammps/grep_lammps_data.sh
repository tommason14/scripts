#!/bin/sh

# Grep data from all log.lammps* files or lammps*.out.
# Sort by timestep, the first value of the thermo_style command
# Print to first argument, or data.csv otherwise.

[ $# -eq 0 ] || [ $1 == "lammps.out" ] || [ $1 == "log.lammps" ] && 
echo "Pass in log or lammps to take data from log.lammps* or lammps*.out
Syntax: $(basename $0) [log|lammps] <csvfile>" && exit 1

printf "Reminder that if you have multiple log files, 
renumber the first lammps.out or log.lammps to lammps1.out or log.lammps1 ...
Continue? [Y] "
read option
[ ! "$option" = "Y" ] && [ ! "$option" = "y" ] && [ ! "$option" = "" ] && exit 1

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

choice=$1
output=${2-data.csv}
echo "Writing lammps data to $output"

# data- assumes Step is the first value of thermo_style.
# If log is not the first file, and rerun was started with a 
# data file, the Step number starts at 0, when we want it to start at the
# next timestep. So append to data.tmp sequentially.
[ -f data.tmp ] && rm data.tmp
for f in $(ls "$choice"* | sort -n)
do
  lines="$(sed -n '/Step/,/Loop/p' $f |
    grep -v 'Step\|Loop' |
    grep '^\s\+[0-9]\+\s\+-*[0-9]\+' |
    $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g')"
  [ -f data.tmp ] && [ $(cat data.tmp | wc -l) -gt 0 ] && (
  # take last two lines and find step increment
  laststep="$(tail -n 1 data.tmp | awk -F "," '{print $1}'| $sed 's/^\s*//;s/\s*$//')"
  eqn="$(tail -n 2 data.tmp |
    awk -F"," '{print $1}' |
    $sed '/^\s*$/d' |
    tr '\n' '-'|
    $sed 's/\-$//')"
    increment="$(echo "$eqn" | bc | $sed 's/^-//')"
    echo "$lines" |
    awk -F"," -v last="$laststep" -v inc="$increment" '$1+=last+inc' OFS="," >> data.tmp
  ) || printf "%s\n" "$lines" > data.tmp
done

# find the header from first output
ls "$choice"* |
  sort -n |
  head -1 |
  xargs grep Step |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > header.tmp
cat header.tmp data.tmp > $output

rm header.tmp data.tmp
