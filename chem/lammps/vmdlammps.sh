#!/usr/bin/env bash

# Opens a lammps data file in VMD

[ $# -eq 0 ] && echo "Syntax: $(basename $0) data.lmps" && exit 1

data="$1"
[ ! -f "$data" ] && echo "Error: $data doesn't exist" && exit 1

# Open directly if possible. If a datafile is entered, use topotools
topo(){
[ -f tmpvmd ] && rm tmpvmd
cat << EOF >> tmpvmd
package require topotools
topo readlammpsdata "$data"
EOF

vmd -e tmpvmd

rm tmpvmd
}

grep -Fq 'xlo xhi' "$data" && topo && exit 1 ||
[[ $data =~ lmp$ || $data =~ lmps$ ]] && vmd -lammpstrj "$data" || vmd "$data"

