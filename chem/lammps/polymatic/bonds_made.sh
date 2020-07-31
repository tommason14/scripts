#!/usr/bin/env bash

[[ $# -eq 0 ]] &&
echo "Counts the number of bonds the modified Polymatic code makes.
Syntax $(basename $0) polymatic.out" && exit 1

# [[ $USER =~ (tmas0023|tommason) ]] && grep='ggrep' || grep='grep'
# $grep -o "([a-z].*)$" "$1" | sed 's/[\(\)]//g' | sort | uniq -c

# cross-compatibility
grep "Pair" "$1" | sed 's/^\s\+//' | cut -d ' ' -f 5- | sed 's/[\(\)]//g' | sort | uniq -c
