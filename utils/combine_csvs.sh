#!/bin/sh

if [ $# -eq 0 ] || [ "$1" = "-h" ]
then 
cat << ENDHELP
Pass in csv files that you would like to combine.
i.e. combine_csvs file1 file2
Note: Assumes that header lines are the same in all csvs.
Result is passed to STDOUT
ENDHELP
exit 
fi

header=$(cat "$@" | head -1)
echo "$header"
cat "$@" | grep -v "$header" | sort -t ',' -k 1 
