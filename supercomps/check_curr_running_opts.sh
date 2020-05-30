#!/bin/sh

[ $USER = "tm3124" ] && echo "Error: Not implemented for PBS yet" && exit 1

for f in $(squeue -u $USER | grep $USER | awk '{print $1}')
do 
cd $(scontrol show jobid $f | grep WorkDir | awk -F '=' '{print $2}')
pwd
if [[ -f opt-freq.log ]] 
then 
  cat opt-freq.log | plotgauss
else
  if grep -iq 'FMO' opt.log 
  then
  cat opt.log | plotrmsgamess
  else
  cat opt.log | plotmp2
  fi
fi
done
