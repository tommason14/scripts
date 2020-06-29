#!/usr/bin/env bash


make_pieda(){
# Modify original fmo3 file
pol="free-state-polarisation"
pieda="full-pieda"
[[ ! -d $pieda ]] && mkdir $pieda
inp="$pieda/spec.inp"
cp spec.inp $inp

# PIEDA can only be run with FMO2
sed -i 's/NBODY=3/NBODY=2/' $inp

# add FMOPRP info from Dmitri's examples in tools/fmo/published
# and efmo0(1) info from the fmo0/spec.dat
# and epl0ds(1) and eint0(1) from free-state-polarisation/spec.dat
efmo=$(sed -n '/efmo0(1)/,/ecorr(1)/p' fmo0/spec.dat | grep -v '^\s*e')
unpol_densities=$(sed -n '/epl0ds/,/eint0/p' free-state-polarisation/spec.dat | grep -v '^\s*e')
ints=$(grep '^\s*eint0' free-state-polarisation/spec.dat |
  awk '{OFS="   "; $1=""; print $0}')
# include newlines by literally typing the string on multiple lines
new_fmoprp=" \$FMOPRP MODORB=3 IPIEDA=2
 EFMO0(1)=
$efmo
 EPL0DS(1)=
$unpol_densities
 EINT0(1)=
$ints
 \$END"
# print before FMOPRP, add new section then print after and including MP2
sed '/$FMOPRP/q' $inp > tmpinp
echo "$new_fmoprp" >> tmpinp
sed -n -e '/$MP2/,$p' $inp >> tmpinp
mv tmpinp $inp
# new job script - assumes 4 logical nodes, 2 physical nodes
cat << ENDJOB > $pieda/spec.job
#!/bin/sh
#PBS -P k96
#PBS -l mem=768gb
#PBS -l ncpus=192
#PBS -l jobfs=500gb
#PBS -l walltime=1:00:00
#PBS -l wd

module load gamess/19-srs

rungms.gadi.ln spec.inp \$PBS_NCPUS > spec.log
ENDJOB
}

cwd="$(pwd)"
for f in $(find . -maxdepth 2 -name "spec.log")
  do
    dir="$(dirname "$f")"
    echo "Entering $dir"
    cd "$dir"
    make_pieda
    cd "$cwd"
  done
