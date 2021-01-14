#!/bin/sh

[ $# -ne 1 ] && echo "Error: pass in filename" && exit 1
[ ! -f $1 ] && echo "$1 doesn't exist" && exit 1

echo "Monomers -> $(grep iFrag $1 | sed '/jFrag/d' | wc -l)" 
echo "Dimers -> $(grep jFrag $1 | sed '/kFrag/d' | wc -l)" 
echo "Trimers -> $(grep kFrag $1 | wc -l)" 
