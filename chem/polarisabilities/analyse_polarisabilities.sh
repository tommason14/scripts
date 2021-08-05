#!/usr/bin/env bash
field="${1:-0.0008}" # by default, 0.0008 au field strength assumed
# numpy linear algebra deprecation warning is sent to /dev/null
[[ ! -f connected.in ]] &&
echo "Run $(dirname $0)/make_connected.in.py to create a file with bonding information" &&
echo "and place this in the current directory." &&
exit 1
echo -e "$field\n" | python3 $(dirname $0)/analyse_polarisabilities.py 2>/dev/null | sed -n '/Total polarizability:/,/-------/p' > polarisability.out
