#!/usr/bin/env bash

make_pl0(){
# Modify original fmo3 file
pol="free-state-polarisation"
[[ ! -d $pol ]] && mkdir $pol
inp="$pol/spec.inp"
job="$pol/spec.job"
cp spec.inp $inp
cp spec.job $job

# PIEDA can only be run with FMO2
sed -i 's/NBODY=3/NBODY=2/' $inp
# so remove RITRIM
sed -i '/RITRIM/d' $inp

# add FMOPRP info from Dmitri's examples in tools/fmo/published
# and the efmo0(1) info from the fmo0/spec.dat
efmo=$(sed -n '/efmo0(1)/,/\$end/p' fmo0/spec.dat | grep -v 'efmo\|\$end')
# include newlines by literally typing the string on multiple lines
new_fmoprp=" \$FMOPRP MAXIT=1 IREST=2 MODORB=3 IPIEDA=2 MODPAR=77
 EFMO0(1)=
$efmo
 \$END"
# print before FMOPRP, add new section then print after and including MP2
sed '/$FMOPRP/Q' $inp > "$pol/tmpinp"
echo "$new_fmoprp" >> "$pol/tmpinp"
sed -n -e '/$MP2/,$p' $inp >> "$pol/tmpinp"
mv "$pol/tmpinp" $inp 

# need F30 file from fmo0 run
cp fmo0/tempfiles/spec.F30 $pol/spec.F40.000

sed -i "s|rungms.monarch|~/gamess16-srs/rungms.monarch.copyF40|" $job
}


cwd="$(pwd)"
for f in $(find . -maxdepth 2 -name "spec.log")
  do
    dir="$(dirname "$f")"
    echo "Entering $dir"
    cd "$dir"
    make_pl0
    cd "$cwd"
  done
