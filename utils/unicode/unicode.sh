#!/usr/bin/env bash

sed 's/^U\+/\\u/' "$(dirname $0)/unicode.txt" | while read -r line; do echo -e $line ; done | fzf --color=bg+:-1 --preview "echo -e {} | sed 's/, /\n/' | sed 's/, /\n/'" | awk -F',' '{print $1}' | tr -d '\n' | pbcopy
