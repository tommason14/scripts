#!/usr/bin/env bash

[[ $# -ne 1 ]] && echo "Syntax: $(basename $0) logfile" && exit 1

module list | grep -q orca || module load orca/4.2.1

##########################################################
#  expect script to run through orca_plot automatically  #
##########################################################

cat << EOF > tmp.expect
#!/usr/bin/env expect

set file [lindex \$argv 0]
set orbital [lindex \$argv 1]

spawn orca_plot \$file -i

# orbital

expect "Enter a number:"
send "2\r"
expect "Enter MO:"
send "\$orbital\r"

# Cube output

expect "Enter a number:"
send "5\r"
expect "Enter Format:"
send "7\r"

# Generate file

expect "Enter a number:"
send "10\r"

# Exit

expect "Enter a number:"
send "11\r"
EOF
chmod +x tmp.expect

#################################
#  Parse log file for orbitals  #
#################################

log="$1"
gbw="${log%.*}.gbw"

# Orca prints the orbital occupation along with the orbtial energies, so can use them to find the 
# homo-1, homo, lumo and lumo+1. Then generate cube files for those orbitals

lumo=$(sed -n '/ORBITAL ENERGIES/,/MOLECULAR ORBITALS/p' "$log" | grep "^\s*[0-9]\+\s\+0.0000" | head -1 | awk '{print $1}')
orbs="$(($lumo - 2)) $(($lumo - 1)) $lumo $(($lumo + 1))"
for orb in ${orbs[@]}
do
  echo "Generating cube for orbital $orb"
  ./tmp.expect $gbw $orb >/dev/null
done

###########################################
#  rename cube files to indicate orbital  #
###########################################

cubes=($(find . -maxdepth 1 -name "*cube" | sort -n))
mv "${cubes[0]}" "homo-minus-1.cube"
mv "${cubes[1]}" "homo.cube"
mv "${cubes[2]}" "lumo.cube"
mv "${cubes[3]}" "lumo-plus-1.cube"

rm tmp.expect
