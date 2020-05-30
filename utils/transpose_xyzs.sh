#!/usr/bin/env sh

for f in $(ls *xyz)
do
echo $f
tail -n+3 $f |\
awk '{print $1}' |\
sed "s/^/'/;s/$/'/" |\
tr '\n' ', ' | sed 's/^/[/;s/,$/]/'
echo
done
