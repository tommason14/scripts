#!/bin/sh

if [ $# -ne 1 ]
then
echo "Error: pass in filename"
exit 1
fi

echo "Monomers -> $(grep iFrag $1 | sed '/jFrag/d' | wc -l)" 
echo "Dimers -> $(grep jFrag $1 | sed '/kFrag/d' | wc -l)" 
echo "Trimers -> $(grep kFrag $1 | wc -l)" 
