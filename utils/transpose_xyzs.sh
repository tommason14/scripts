#!/usr/bin/env sh

# take an xyz and print out each atomic symbol in a list of strings
for f in $(ls *xyz)
do
echo $f
tail -n+3 $f |\
awk '{print $1}' |\
sed "s/^/'/;s/$/'/" |\
tr '\n' ', ' | sed 's/^/[/;s/,$/]/'
echo
done
