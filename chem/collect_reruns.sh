#!/usr/bin/env bash

# Assumes a directory structure of
# .
# ├── struct1
# ├── struct2

printf "Reminder to run autochem -e. Continue? [Y]"
read option
if [[ $option =~ (y|Y) ]] || [[ $option = "" ]]   
then
  for f in $(find . -type d -maxdepth 1 | sed '/^\.$/d;/equilibrated/d;/reruns/d;s/.\///')
  do
    spec=$(find "$f" -path "*spec*xyz")
    if [[ "$spec" == "" ]]
    then
      echo "No equilibrium found for $f"
      rerun=$(find "$f" -path "*rerun*xyz" | sort -r | head -1)
      [[ ! -d reruns ]] && mkdir reruns
      cp $rerun "reruns/$f.xyz"
      echo "Copied $rerun to reruns/$f.xyz"
    fi
  done
fi
