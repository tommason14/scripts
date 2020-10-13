#!/usr/bin/env bash 

GAFF="../create_gaff_jobs"

packinp=$(ls pack.inp 2> /dev/null | wc -l)
gaff=$(ls gaff.ff 2> /dev/null | wc -l)
polym=$(ls polym.in 2> /dev/null | wc -l)

[[ $packinp -eq 0 || $gaff -eq 0 || $polym -eq 0 ]] &&
echo "Make sure you have the following files in this dir:" &&
echo "- pack.inp" &&
echo "- gaff.ff" &&
echo "- polym.in" &&
echo "Make sure that pack.inp outputs a file named pack.xyz after" &&
echo "packmol has finished" &&
exit 1

# need additional params for polymatic
mkdir small
cp pack.inp small/pack.inp
cp gaff.ff small/
# just 1 molecule for each different structure
$sed -i 's/number.*/number 1/' small/pack.inp
ls *xyz | xargs -I{} cp {} small/
cp polym.in small/
cd small
# clear output
{ packmol < pack.inp
lmp_gaff.py pack.xyz gaff.ff # creates pack.lmps
} > /dev/null
echo "Finding parameters that Polymatic needs. Could take about a minute..."
add_additional_params.py -l pack.lmps -f gaff.ff -p polym.in
# additional params now in small/pack.lmps
cd ..

## desired system ##
echo "Creating desired system"
{
packmol < pack.inp # generate pack.sh
lmp_gaff.py pack.xyz gaff.ff
} > /dev/null
# resulting pack.lmps assumes all atoms are in the same molecule, and charges
# are incorrect, so fix these
change_molecule_id.py pack.lmps
add_correct_charges.py pack.lmps # asks user to assign _spec.log files to each molecule
# Change params in desired lammps datafile

echo "Adding additional parameters to datafile"
newtypes="$(grep 'types' small/pack.lmps)"
newparams="$(sed -n '/Masses/,/Atoms/p' small/pack.lmps | grep -v '^\s*Atoms')"
sed '/impropers/q' pack.lmps > tmp.start
grep 'lo.*hi' pack.lmps > tmp.boxsize
sed -n '/Atoms/,//p' pack.lmps > tmp.end

cat tmp.start <(echo "$newtypes") tmp.boxsize <(printf "\n$newparams\n\n") tmp.end > tmp.final
mv tmp.final pack.lmps
# Clean up
rm tmp* 
rm -rf small
polymatic_types.py pack.lmps > types.txt
generate_molecule_id_file.py
