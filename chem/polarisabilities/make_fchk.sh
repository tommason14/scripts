#!/usr/bin/env bash


cwd=$(pwd); for f in $(find . -name spec.chk); do dir=$(dirname $f); [[ ! -f $dir/spec.fchk ]] && echo "$dir" && cd $dir && formchk spec.chk && cd $cwd; done
