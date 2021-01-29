#!/usr/bin/env bash

[[ -t 0 ]] && echo "Pipe in a list of job files in subdirectories to change.
This script changes the job name to the topmost directory.
For example:
$ find . -name spec.job | $(basename $0)
takes ./dirname/spec.job and sets #SBATCH --job-name=dirname" && exit 1

while read file
do
  name=$(echo $file | cut -d/ -f2)
  sed -i "s/job-name=.*/job-name=$name/" $file
done < /dev/stdin
