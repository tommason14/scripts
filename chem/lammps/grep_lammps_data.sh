#!/bin/sh

# Grep data from all lammps*.out.
# Sort by timestep, the first value of the thermo_style command

[ "$1" = "-h" ] && echo "Returns data from lammps*out files." && exit 1

if [ $(find . -maxdepth 1 -name "lammps*out" | wc -l) -gt 1 ]
then
printf "You have multiple log files, make sure you renumber the first lammps.out to lammps1.out.
Continue? [Y] "
read option
[ ! "$option" = "Y" ] && [ ! "$option" = "y" ] && [ ! "$option" = "" ] && exit 1
fi

[ $USER == "tommason" ] || [ $USER == "tmas0023" ] && sed="gsed" || sed="sed"

# data- assumes Step is the first value of thermo_style.
# If log is not the first file, and rerun was started with a 
# data file, the Step number starts at 0, when we want it to start at the
# next timestep. So append to data.tmp sequentially.
[ -f data.tmp ] && rm data.tmp
for f in $(find . -maxdepth 1 -name "lammps*out" | sort -n)
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
first=$(find . -maxdepth 1 -name "lammps*out" | sort -n | head -1)
ls "$first" |
  sort -n |
  head -1 |
  xargs grep Step |
  $sed 's/^\s\+//;s/\s\+$//;s/\s\+/,/g' > header.tmp
cat header.tmp data.tmp 

rm header.tmp data.tmp
