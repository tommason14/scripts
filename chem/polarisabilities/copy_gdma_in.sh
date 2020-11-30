#!/bin/sh
thisdir=$(dirname $0)
for f in $(find . -name spec.fchk); do dir=$(dirname $f); [[ ! -d $dir/gdma ]] && mkdir $dir/gdma && cp $thisdir/gdma.in $dir/gdma/gdma.in && cp $thisdir/gdma.job $dir/gdma/gdma.job ; done
