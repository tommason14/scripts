#!/usr/bin/env bash

# Take in a labelled xyz of a single molecule along with
# a GAMESS geodesic charge log file of the same xyz. Also
# need a 'stock' force field that will not contain the
# correct partial charges, and a pack.inp file in the
# same directory

[ $# -lt 6 ] &&
  echo "Syntax: $(basename "$0") -c geodesic.log -x molecule.xyz -f gaff.ff [-s mon|m3]" &&
  echo "By default, a job script for monarch is created." &&
  echo "Also assumes that there is a pack.inp file for packmol to use." &&
  exit 1

while [ $# -gt 0 ]
do
  option="$1"
  case "$option" in
    -c) charges="$2" 
        shift
        shift
        ;;
    -x) mon_xyz="$2"
        shift
        shift
        ;;
    -f) ff="$2"
        shift
        shift
        ;;
    -s) sc="$2"
        shift
        shift
        ;;
  esac
done

scriptlocation="$(dirname "$0")"
mon_data="${mon_xyz%.*}.data"
mod_ff="gaff_modified.ff"

if [ "$USER" = "tommason" ] || [ "$USER" = "tmas0023" ] 
then
  sed="gsed"
  grep="ggrep"
else
  sed="sed"
  grep="grep"
fi

# make 'monomer' data file with wrong charges - creates $mon_data
"$scriptlocation"/lmp_gaff.py "$mon_xyz" "$ff"

# add partial charges to $ff- creates $mod_ff
"$scriptlocation"/average_partial_charges.py -c "$charges" -d "$mon_data" -f "$ff"

# make the large system
packmol < pack.inp

# now create data file of the large system, using the modified force field
"$scriptlocation"/lmp_gaff.py system.xyz "$mod_ff"

"$scriptlocation"/change_molecule_id.py system.data 

# add molecule IDs to the large system file...otherwise they will all be 1
# currently only works for one molecule in pack.inp!
# num_mols=$(grep number pack.inp | awk '{print $2}')
# num_atoms=$(head -1 system.xyz)
# repeat=$(( $num_atoms / $num_mols ))
# mol_id=$(for i in $(seq $num_mols); do for j in $(seq $repeat); do echo $i; done; done)

# # replace columns in atom sections
# $sed '/Atoms/q' system.data > tmpb4atoms
# $grep -E '^([0-9]+\s+){3}(-?[0-9]+\.[0-9]+\s+){4}#\s?[A-z0-9]' system.data > tmpatoms
# paste <(awk '{print $1}' tmpatoms) <(printf '%s\n' $mol_id) <(awk '{for(i=3; i <= NF; ++i) printf "%s ", $i; print ""}' tmpatoms) > tmptmpatoms
# # format atoms section
# paste <(sed 's/#.*//' tmptmpatoms | column -t) <($grep -Po '\#.*' tmptmpatoms) | $sed 's/\s\+\#/   \#/' > tmpatoms
# sed -ne '/Bonds/,$ p' system.data > tmpafteratoms
# cat tmpb4atoms <(echo) tmpatoms <(echo) tmpafteratoms > system.data
# rm tmpb4atoms tmpatoms tmptmpatoms tmpafteratoms

[ ! -d rundir ] && mkdir rundir
mv system.data rundir/

[ ! -d intermediates ] && mkdir intermediates
mv "$mon_data" "$mod_ff" system.xyz intermediates/

##############################
#  Inputs and modifications  #
##############################

cat << ENDINP >> rundir/in.lmp
# GAFF NPT simulation

units           real
boundary        p p p
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

atom_style      full
bond_style      harmonic
angle_style     harmonic
dihedral_style  fourier
improper_style  cvff

pair_style      lj/charmm/coul/long 9.0 10.0 10.0
kspace_style    pppm 0.0001
pair_modify     mix arithmetic
special_bonds   amber

read_data       system.data

thermo_style    custom step cpu etotal ke pe evdwl ecoul elong temp press vol density
thermo          1000

dump d1 all     custom 1000 traj.lmp element xu yu zu
dump_modify d1  element elements
dump_modify d1  sort id

minimize        0.0 1.0e-8 1000 100000

fix             SHAKE all shake 0.0001 20 0 b bonds
velocity        all create 298 298 dist gaussian
timestep        1
fix 8 all       npt temp 298 298 200 iso 1 1 1000 tchain 3 pchain 3 mtk yes
restart         1000 sim.restart1 sim.restart2
run             200000
undump d1
unfix 8

write_data output.data
ENDINP

# check if $sc is not set, and run monarch by default
if [ -z "$sc" ] || [ "$sc" = "mon" ]
then
cat <<\ENDMONJOB >> rundir/job.slm
#!/bin/bash
#SBATCH --qos=partner
#SBATCH --time=03:00:00
#SBATCH --ntasks=16
#SBATCH --tasks-per-node=16
#SBATCH --mem=16G
#SBATCH --partition=comp,short

module load openmpi
export lammps='/mnt/lustre/projects/p2015120004/apps/clammps/build/lmp_mpi'

srun --export=all -n $SLURM_NTASKS $lammps -in in.lmp > lammps.out
find . -empty -delete
ENDMONJOB
elif [ "$sc" = "m3" ] 
then
cat <<\ENDMASJOB >> rundir/job.slm
#!/bin/bash
#SBATCH --account=sn29
#SBATCH --time=3:00:00
#SBATCH --ntasks=16
#SBATCH --tasks-per-node=16
#SBATCH --mem=16G

module load openmpi/1.10.7-mlx
export lammps='/projects/sn29/apps/clammps/build/lmp_mpi'

srun --export=all -n $SLURM_NTASKS $lammps -in input_file.lmp > lammps.out
find . -empty -delete
ENDMASJOB
fi

#######################
#  Find bonds to fix  #
#######################

bonds=$(
  $sed -n '/Bond Coeffs/,/Angle Coeffs/p' rundir/system.data |
  $grep '# H\|-H' |
  awk '{print $1}' |
  tr '\n' ' '
)

$sed -i "s/b bonds/b $bonds/" rundir/in.lmp

#####################################
#  Sort elements when writing data  #
#####################################

elements=$(
  $sed -n '/Pair Coeffs/,/Bond Coeffs/p' rundir/system.data |
  $grep '#' |
  awk '{print $NF}' |
  $grep -Po "^[A-Z]" |
  tr '\n' ' '
)

$sed -i "s/elements/$elements/" rundir/in.lmp
