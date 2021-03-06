#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Plots option against timesteps using gnuplot.
Enter column number when prompted.
NB: Assumes that the timestep is the first value in the thermo_style command.

If gnuplot gives a warning of 'slow font initialisation', run 'fc-cache' and
try again.

Should work on multiple files with the same thermo_style output
Syntax: $(basename $0) lammps.out [lammps.out2 ...]" && exit 1

for f in "$@"
do
  [[ ! -f $f ]] && echo "$f not found" && exit 1
done

# remove Step from choice for user, then add one to their choice to find correct column
opts=$(grep Step $1 | tail -1 | tr ' ' '\n' | tail -n +2 | nl)

# read stdin so that "echo density | lammps_plot.sh lammps.out" works
if [[ -p /dev/stdin ]]
then
  option="$(cat /dev/stdin)"
  option=$(printf "%s\n" "${opts[@]}" | grep -i $option | awk '{print $1}')
  [[ $option == "" ]] &&
    echo -e "Option not found. Possible choices are:\n${opts[@]}" &&
    exit 1
else
  while read num opt; do
    echo $num $opt
  done < <(printf "%s\n" "${opts[@]}" | xargs -n2)
  printf "Option (enter number): "
  read option
fi

label=$(grep Step $1 | tail -1 | awk -v choice=$((option + 1)) '{print $choice}' | sed 's/_/-/') 
# gnuplot treats _ as subscript so replace them with hyphens

# gadi gnuplot qt terminal always shows "slow font initialisation warning", so change to x11 terminal.
# This might show stretched text as though the font is still too slow to render, but at least the location
# of the text is correct
[[ $HOSTNAME =~ gadi ]] && extra='set terminal x11;' || extra=''

grep_lammps_data.sh "$@" |
  sed '1d' |
  awk -F"," -v step=1 -v choice=$((option + 1)) '{print $step,$choice}' |
  gnuplot --persist -e "$extra set xlabel 'Timestep'; set ylabel '$label'; plot '-' using 1:2 with lines notitle"
