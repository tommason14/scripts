#!/usr/bin/env bash
thisdir=$(dirname $0)
[[ $USER == "tm314" ]] && job='gdma_gadi.job' || job='gdma_mon.job'
for f in $(find . -name spec.fchk); do dir=$(dirname $f); [[ ! -d $dir/gdma ]] && mkdir $dir/gdma && cp $thisdir/gdma.in $dir/gdma/gdma.in && cp $thisdir/$job $dir/gdma/gdma.job ; done
