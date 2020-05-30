#!/bin/sh

files=$(ls -1 *xyz)
columns=("Files" ${@})
num_cols="${#columns[@]}"
header="$(printf ",%s" ${columns[@]})"
echo $header | cut -d ',' -f 2-
for f in $files
do
printf "$f"
echo $(yes "," | head -n $(( $num_cols - 1 )))
done
