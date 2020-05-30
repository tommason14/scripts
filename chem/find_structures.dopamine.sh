#!/bin/sh

structures="$HOME/Documents/repos/dopamine/dopamine_analysis/elucidation_of_structure_in_c2mim_ac/equilibrated-structures"
for f in $(ls -1 *xyz)
do 
if [ -f "$structures/$f" ] 
then
  echo "$f in \$structures" 
else
  cp "$f" $structures/
fi
done
