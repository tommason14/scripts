#!/bin/sh

[ $# -eq 0 ] || [ "$1" == "-h" ] &&
  echo "Plot 'gmx energy' output with gnuplot." && 
  echo "Syntax: $(basename $0) file.xvg" &&
  exit 1

tail -n +24 $1 | gnuplot --persist -e "set terminal x11; plot '-' with lines notitle"
