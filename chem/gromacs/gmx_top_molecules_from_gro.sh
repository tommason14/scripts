 #!/usr/bin/env bash

[[ $# -ne 1 || $1 =~ '-h' ]] && echo "Syntax: $(dirname $0) file.gro" && exit 1
[[ -f $1 ]] || { echo "Error: $1 not found"; exit 1; }

: '
Search through a gromacs gro file and print out all residues in the order they are included.
This works for cases where the system is built with gmx insert-molecules which may have repeating units 
such as:
pol 1
SOl 124
NA  24
pol 1
SOL 124
NA 24
etc...
'

check_columns(){
# mstools can write 1c2c1i or 100 dhi depending on the length of the resname,
# so check if column 1 contains letters, otherwise merge columns 1 and 2
awk '{if($1 ~ /[A-z]/){print $1;} else {print $1$2}}'
}

[[ $USER =~ (tmas0023|tommason) ]] && sed='gsed' || sed='sed'

sed 1d $1 | sed '1d;$d' | check_columns  | uniq | $sed 's/[0-9]\+//' | uniq -c | awk '{print $2," ",$1}' | column -t
