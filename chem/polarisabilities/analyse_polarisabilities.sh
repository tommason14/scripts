#!/usr/bin/env bash
field="${1:-0.0008}" # by default, 0.0008 au field strength assumed
echo -e "yes\nyes\n$field\nyes\n" | python3 $(dirname $0)/analyse.py | sed -n '/Total polarizability:/,/-------/p' > polarisability.out
