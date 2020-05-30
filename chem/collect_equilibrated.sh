#!/bin/sh

# Assumes a structure such as 
# .
# ├── struct1
# │   ├── opt.dat
# │   ├── opt.inp
# │   ├── opt.job
# │   ├── opt.log
# │   └── spec
# │       └── opt_equil.xyz
# ├── struct2
# │   ├── opt.dat
# │   ├── opt.inp
# │   ├── opt.job
# │   ├── opt.log
# │   └── spec
# │       └── opt_equil.xyz

# so run chem_assist -e first

newdir="equilibrated"
[ ! -d $newdir ] && mkdir $newdir

find_xyz(){
    struct="$1"
    xyz=$(find . -path "*/$struct/*spec*xyz")
    if [ ! "$xyz" = "" ]
    then 
      echo "$struct -> $xyz" && cp "$xyz" "$newdir/$struct.xyz"
    else
      echo "No equilibrated structure for $struct"
    fi  
}

structs=$(find . -maxdepth 1 -type d | grep -v "^\.$" | awk -F'/' '{print $2}' | grep -v "$newdir")
while read -r struct
do
find_xyz $struct
done <<< "$structs"
