#!/usr/bin/env bash
field="${1:-0.0008}" # by default, 0.0008 au field strength assumed
# numpy linear algebra deprecation warning is sent to /dev/null
echo -e "yes\nyes\n$field\nyes\n" | python3 $(dirname $0)/analyse.py 2>/dev/null | sed -n '/Total polarizability:/,/-------/p' > polarisability.out
