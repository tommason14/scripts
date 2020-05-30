#!/bin/sh

# Opens a lammps data file in VMD

[ $# -eq 0 ] && echo "Syntax: $(basename $0) data.lmps" && exit 1

data="$1"
[ ! -f "$data" ] && echo "Error: $data doesn't exist" && exit 1

[ -f tmpdvmd ] && rm tmpvmd
cat << EOF >> tmpvmd
package require topotools
topo readlammpsdata $data
EOF

if [ "$USER" = "tommason" ]
then
  vmd="/Applications/VMD 1.9.4a38.app/Contents/vmd/vmd_MACOSXX86_64"
elif [ "$USER" = "tmas0023" ]
then
  vmd="/Applications/VMD 1.9.3.app/Contents/vmd/vmd_MACOSXX86"
else
  vmd="vmd" # module loaded
fi

"$vmd" -e tmpvmd

rm tmpvmd
