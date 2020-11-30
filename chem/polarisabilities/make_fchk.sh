#!/bin/sh
module load gaussian/g16c01

cwd=$(pwd); for f in $(find . -name spec.chk); do dir=$(dirname $f); [[ ! -f $dir/spec.fchk ]] && echo "$dir" && cd $dir && formchk spec.chk && cd $cwd; done
