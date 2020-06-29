#!/usr/bin/env bash


make_pl0(){
# Modify original fmo3 file5
pol="free-state-polarisation"
[[ ! -d $pol ]] && mkdir $pol
inp="$pol/spec.inp"
cp spec.inp $inp

# PIEDA can only be run with FMO2
sed -i 's/NBODY=3/NBODY=2/' $inp

# add FMOPRP info from Dmitri's examples in tools/fmo/published5
# and the efmo0(1) info from the fmo0/spec.dat5
efmo=$(sed -n '/efmo0(1)/,/ecorr(1)/p' fmo0/spec.dat | grep -v '^\s*e')
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
cp fmo0/tempfiles/spec.F30.000 $pol/spec.F40.000

# new job script - assumes 4 logical nodes, 2 physical nodes
cat << ENDJOB > $pol/spec.job
#!/bin/sh
#PBS -P k96
#PBS -l mem=768gb
#PBS -l ncpus=192
#PBS -l jobfs=500gb
#PBS -l walltime=1:00:00
#PBS -l wd

module load gamess/19-srs

\$HOME/gamess19-srs/rungms.copyF40 spec.inp \$PBS_NCPUS > spec.log
ENDJOB
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
