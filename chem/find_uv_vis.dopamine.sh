#!/bin/sh

orcajob() {
  [ -f "${1%.log}.gbw" ] && return 0 || return 1 
# return 1 is counted as an
# error, doesn't pass if statement
}

uv="$HOME/Documents/repos/dopamine/dopamine_analysis/elucidation_of_structure_in_c2mim_ac/uv-vis-spectra"
for f in $(find . -name "spec.log")
do 
if orcajob "$f" 
then
  name=$(echo "$f" | awk -F '/' '{print $2}') # parent dir name is molecule
  echo "Checking for $name"
  if [ -f "$uv/$name.log" ] 
  then
    echo "$name in \$uv" 
  else
    cp "$f" "$uv/$name.log"
  fi
fi
done
