#!/usr/bin/env bash

if ! command -v formchk &>/dev/null
then
  echo "Error: formchk command not found - load a gaussian module if possible"
  exit 1
fi

cwd=$(pwd); for f in $(find . -name spec.chk); do dir=$(dirname $f); [[ ! -f $dir/spec.fchk ]] && echo "$dir" && cd $dir && formchk spec.chk && cd $cwd; done
