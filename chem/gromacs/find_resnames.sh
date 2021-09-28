#!/usr/bin/env bash

[[ $# -ne 1 || $1 =~ '-h' ]] && echo "Syntax: find_resnames.sh file.gro" && exit 1
[[ -f $1 ]] || { echo "Error: $1 not found"; exit 1; }

check_columns(){
# mstools can write 1c2c1i or 100 dhi depending on the length of the resname, 
# so check if column 1 contains letters
awk '{if($1 ~ /[A-z]/){print $1;} else {print $2}}'
}

rm_resids(){
sed 's/^[0-9]\+//'
}

sed 1d $1 | sed '1d;$d' | check_columns | rm_resids | sort | uniq | paste -s -d ' '
