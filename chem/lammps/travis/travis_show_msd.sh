#!/usr/bin/env bash

# remove colour codes and print diffusion info from travis log file
# obtained by travis ... > file
# IMPORTANT: sed expression is probably bsd only, but seems to work okay on Linux too 

sed -E "s/"$'\E'"\[([0-9]{1,3}((;[0-9]{1,3})*)?)?[m|K]//g" "${1:-/dev/stdin}" | sed -n '/Saving result/,/assuming/p' |
while read -r line 
do 
  [[ $line =~ msd ]] && echo $line | awk -F'_' '{printf "%s ",$2}'
  [[ $line =~ 'm^2/s' ]] && echo $line | awk '{printf "%.6E m^2/s\n", $(NF-1)}'
done
