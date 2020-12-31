#!/bin/sh

# remove colour codes and print diffusion info from travis log file
# obtained by travis ... > file
# IMPORTANT: sed expression is mac os (bsd) only

sed -E "s/"$'\E'"\[([0-9]{1,3}((;[0-9]{1,3})*)?)?[m|K]//g" "${1:-/dev/stdin}" | sed -n '/Saving result/,/assuming/p'
