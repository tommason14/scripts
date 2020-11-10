#!/bin/sh

[ $# -ne 1 ] || [ "$1" = "-h" ] && echo "Syntax: $(basename $0) deform.txt" && exit 1

echo "Strain Stress_ave" > average_deformation.txt
sed '1d' $1 | awk '{ave=(($2+$3+$4)/3); print $1,ave}' >> average_deformation.txt
