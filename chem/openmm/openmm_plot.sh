#!/usr/bin/env bash

[[ $# -eq 0 || $1 == '-h' ]] &&
echo "Plots option against timesteps using gnuplot.
Enter column number when prompted.
NB: Assumes that the timestep is the first value set in the StateDataReporter command.

If gnuplot gives a warning of 'slow font initialisation', run 'fc-cache' and
try again.

Works with both stdout and csv formatted output

Syntax: $(basename $0) output" && exit 1

[[ ! -f $1 ]] && echo "$1 not found" && exit 1

plot_with_csv(){
  # remove Step from choice for user, then add one to their choice to find correct column
  opts=$(head -1 $1 | tr ',' '\n' | tail -n +2 | nl)

  # read stdin so that "echo density | openmm_plot.sh output.csv" works
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

  label=$(grep Step $1 | tail -1 | awk -F"," -v choice=$((option + 1)) '{print $choice}' | sed 's/_/-/g;s/"//g') 
  # gnuplot treats _ as subscript so replace them with hyphens, and remove the quotation marks needed when writing csv headers with spaces

  # gadi gnuplot qt terminal always shows "slow font initialisation warning", so change to x11 terminal.
  # This might show stretched text as though the font is still too slow to render, but at least the location
  # of the text is correct
  [[ $HOSTNAME =~ gadi ]] && extra='set terminal x11;' || extra=''

  tail -n +2 $1 |  
    awk -F"," -v step=1 -v choice=$((option + 1)) '{print $step,$choice}' |
    gnuplot --persist -e "$extra set xlabel 'Timestep'; set ylabel '$label'; plot '-' using 1:2 with lines notitle"
}

plot_with_stdout(){
  # remove all lines before Step
  sed -n '/Step/,$p' $1 > cleaned.tmp

  # remove Step from choice for user, then add one to their choice to find correct column
  opts=$(head -1 cleaned.tmp | sed -E 's/\s+/\n/g' | tail -n +2 | nl)

  # read stdin so that "echo density | openmm_plot.sh output.csv" works
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

  label=$(grep Step cleaned.tmp | tail -1 | awk -v choice=$((option + 1)) '{print $choice}' | sed 's/_/-/g;s/"//g') 
  # gnuplot treats _ as subscript so replace them with hyphens, and remove the quotation marks needed when writing csv headers with spaces

  # gadi gnuplot qt terminal always shows "slow font initialisation warning", so change to x11 terminal.
  # This might show stretched text as though the font is still too slow to render, but at least the location
  # of the text is correct
  [[ $HOSTNAME =~ gadi ]] && extra='set terminal x11;' || extra=''

  tail -n +2 cleaned.tmp |  
    awk -v step=1 -v choice=$((option + 1)) '{print $step,$choice}' |
    gnuplot --persist -e "$extra set xlabel 'Timestep'; set ylabel '$label'; plot '-' using 1:2 with lines notitle"
  rm cleaned.tmp
}


if grep 'Step' $1 | grep -q ','
then
  plot_with_csv $1
else
  plot_with_stdout $1
fi
