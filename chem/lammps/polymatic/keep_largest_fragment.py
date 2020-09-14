#!/usr/bin/env python3
"""
In cases where Polymatic doesn't form a bond with every single monomer,
use this script to keep the largest polymer chain.
In reality, Polymatic will work but maybe one or two monomers (out of hundreds)
will not form a bond- use this script to remove those odd monomers
"""

from collections import Counter
import sys
import subprocess as sp
import re

if len(sys.argv) < 2 or sys.argv[1] == '-h':
    sys.exit('Syntax: polymatic_types.py <lammps_datafile>')
datafile = sys.argv[1]
start, end = datafile.rsplit('.')
output = '.'.join([start, 'largest', end])

found_masses = False
found_pairs = False
found_bond_coeffs = False
found_angle_coeffs = False
found_dihedral_coeffs = False
found_improper_coeffs = False
found_atoms = False
found_bonds = False
found_angles = False
found_dihedrals = False
found_impropers = False

is_data = lambda line: len(line.split()) > 0

# Coeffs
masses = {}
pair_coeffs = {}
bond_coeffs = {}
angle_coeffs = {}
dihedral_coeffs = {}
improper_coeffs = {}

# Molecular properties
atoms = {'id': [], 'mol': [], 'type': [], 'rest': []}
# grep data, then remove IDs of small molecules, and rm bonds/angles etc... containing atoms of those same IDs.
bonds = {'atoms': [], 'type': []}
angles = {'atoms': [], 'type': []}
dihedrals = {'atoms': [], 'type': []}
impropers = {'atoms': [], 'type': []}

with open(datafile) as f:
    for line in f:
        if is_data(line):
            if found_masses and re.search('^\s*[A-Z]', line):
                found_masses = False
            if found_pairs and re.search('^\s*[A-Z]', line):
                found_pairs = False
            if found_bond_coeffs and re.search('^\s*[A-Z]', line):
                found_bond_coeffs = False
            if found_angle_coeffs and re.search('^\s*[A-Z]', line):
                found_angle_coeffs = False
            if found_dihedral_coeffs and re.search('^\s*[A-Z]', line):
                found_dihedral_coeffs = False
            if found_improper_coeffs and re.search('^\s*[A-Z]', line):
                found_improper_coeffs = False
            if found_atoms and re.search('^\s*[A-Z]', line):
                found_atoms = False
            if found_bonds and re.search('^\s*[A-Z]', line):
                found_bonds = False
            if found_angles and re.search('^\s*[A-Z]', line):
                found_angles = False
            if found_dihedrals and re.search('^\s*[A-Z]', line):
                found_dihedrals = False
            if found_impropers and re.search('^\s*[A-Z]', line):
                found_impropers = False
            # Coeffs
            if "Masses" in line:
                found_masses = True
                continue
            if "Pair Coeffs" in line:
                found_pairs = True
                continue
            if "Bond Coeffs" in line:
                found_bond_coeffs = True
                continue
            if "Angle Coeffs" in line:
                found_angle_coeffs = True
                continue
            if "Dihedral Coeffs" in line:
                found_dihedral_coeffs = True
                continue
            if 'Improper Coeffs' in line:
                found_improper_coeffs = True
                continue
            # Molecular properties
            if 'Atoms' in line:
                found_atoms = True
                continue
            if 'Bonds' in line:
                found_bonds = True
                continue
            if "Angles" in line:
                found_angles = True
                continue
            if "Dihedrals" in line:
                found_dihedrals = True
                continue
            if 'Impropers' in line:
                found_impropers = True
                continue
            # Coeffs
            if found_masses:
                if len(line.split()) > 1:
                    count, mass = line.split()
                    masses[count] = mass
            if found_pairs:
                if len(line.split()) > 2:
                    count, *coeffs = line.split()
                    pair_coeffs[count] = coeffs
            if found_bond_coeffs:
                if len(line.split()) > 2:
                    count, *coeffs = line.split()
                    bond_coeffs[count] = coeffs
            if found_angle_coeffs:
                if len(line.split()) > 2:
                    count, *coeffs = line.split()
                    angle_coeffs[count] = coeffs
            if found_dihedral_coeffs:
                if len(line.split()) > 2:
                    count, *coeffs = line.split()
                    dihedral_coeffs[count] = coeffs
            if found_improper_coeffs:
                if len(line.split()) > 2:
                    count, *coeffs = line.split()
                    improper_coeffs[count] = coeffs
            if found_atoms:
                if len(line.split()) > 2:
                    atom_id, mol, atom_type, *rest = line.split()
                    atoms['id'].append(atom_id)
                    atoms['mol'].append(int(mol))
                    atoms['type'].append(atom_type)
                    atoms['rest'].append(rest)
            if found_bonds:
                if len(line.split()) > 2:
                    _, bond_type, *atom_list = line.split()
                    bonds['atoms'].append(atom_list)
                    bonds['type'].append(bond_type)
            if found_angles:
                if len(line.split()) > 2:
                    _, angle_type, *atom_list = line.split()
                    angles['atoms'].append(atom_list)
                    angles['type'].append(angle_type)
            if found_dihedrals:
                if len(line.split()) > 2:
                    _, dihedral_type, *atom_list = line.split()
                    dihedrals['atoms'].append(atom_list)
                    dihedrals['type'].append(dihedral_type)
            if found_impropers:
                if len(line.split()) > 2:
                    _, improper_type, *atom_list = line.split()
                    impropers['atoms'].append(atom_list)
                    impropers['type'].append(improper_type)

# now select fragment with the most atoms
molcounts = Counter(atoms['mol'])
mol_to_keep = molcounts.most_common(1)[0][0]

updated_atom_id = {}  # old: new
new_atoms = {'id': [], 'mol': [], 'type': [], 'rest': []}
# drop fragments
for ind, old_id in enumerate(atoms['id'], 1):
    # loop through atom ids, check if molecule is wanted,
    # and keep if not, noting down the new id of the atom in
    # question (index of the loop).
    # need to update the atom ids of each bond, dih, imp in new datafile
    index_lookup = atoms['id'].index(old_id)
    if atoms['mol'][index_lookup] == mol_to_keep:
        new_atoms['id'].append(ind)
        new_atoms['mol'].append('1')
        new_atoms['type'].append(atoms['type'][index_lookup])
        new_atoms['rest'].append(atoms['rest'][index_lookup])
        updated_atom_id[old_id] = ind

new_bonds = {'atoms': [], 'type': []}
new_angles = {'atoms': [], 'type': []}
new_dihedrals = {'atoms': [], 'type': []}
new_impropers = {'atoms': [], 'type': []}
# check bonds for old ids, carry coeffs over and assign to new id

for ind, atoms in enumerate(bonds['atoms']):
    if any(atom in updated_atom_id for atom in atoms):
        new_bonds['atoms'].append([updated_atom_id[atom] for atom in atoms])
        new_bonds['type'].append(bonds['type'][ind])

for ind, atoms in enumerate(angles['atoms']):
    if any(atom in updated_atom_id for atom in atoms):
        new_angles['atoms'].append([updated_atom_id[atom] for atom in atoms])
        new_angles['type'].append(angles['type'][ind])

for ind, atoms in enumerate(dihedrals['atoms']):
    if any(atom in updated_atom_id for atom in atoms):
        new_dihedrals['atoms'].append(
            [updated_atom_id[atom] for atom in atoms])
        new_dihedrals['type'].append(dihedrals['type'][ind])

for ind, atoms in enumerate(impropers['atoms']):
    if any(atom in updated_atom_id for atom in atoms):
        new_impropers['atoms'].append(
            [updated_atom_id[atom] for atom in atoms])
        new_impropers['type'].append(impropers['type'][ind])

# header info

header = sp.check_output(f"sed '/zhi/q' {datafile}",
                         shell=True).decode('utf-8').split('\n')[:-1]
for ind, line in enumerate(header):
    if 'atoms' in line:
        header[ind] = f'{len(updated_atom_id)} atoms'
    if 'bonds' in line:
        header[ind] = f'{len(new_bonds["atoms"])} bonds'
    if 'angles' in line:
        header[ind] = f'{len(new_angles["atoms"])} angles'
    if 'dihedrals' in line:
        header[ind] = f'{len(new_dihedrals["atoms"])} dihedrals'
    if 'impropers' in line:
        header[ind] = f'{len(new_impropers["atoms"])} impropers'

newfile = header
newfile += ['\nMasses\n']
for count, mass in masses.items():
    newfile += [f' {count} {mass}']
newfile += ['\nPair Coeffs\n']
for count, coeffs in pair_coeffs.items():
    newfile += [f" {count} {' '.join(coeffs)}"]
newfile += ['\nBond Coeffs\n']
for count, coeffs in bond_coeffs.items():
    newfile += [f" {count} {' '.join(coeffs)}"]
newfile += ['\nAngle Coeffs\n']
for count, coeffs in angle_coeffs.items():
    newfile += [f" {count} {' '.join(coeffs)}"]
newfile += ['\nDihedral Coeffs\n']
for count, coeffs in dihedral_coeffs.items():
    newfile += [f" {count} {' '.join(coeffs)}"]
newfile += ['\nImproper Coeffs\n']
for count, coeffs in improper_coeffs.items():
    newfile += [f" {count} {' '.join(coeffs)}"]
newfile += ['\nAtoms\n']
for _id, mol, _type, rest in zip(new_atoms['id'], new_atoms['mol'],
                                 new_atoms['type'], new_atoms['rest']):
    newfile += [f'{_id} {mol} {_type} {" ".join(rest)}']
newfile += ['\nBonds\n']
for ind, atoms in enumerate(new_bonds['atoms']):
    atoms = ' '.join([str(a) for a in atoms])
    newfile += [f" {ind + 1} {new_bonds['type'][ind]} {atoms}"]
newfile += ['\nAngles\n']
for ind, atoms in enumerate(new_angles['atoms']):
    atoms = ' '.join([str(a) for a in atoms])
    newfile += [f" {ind + 1} {new_angles['type'][ind]} {atoms}"]
newfile += ['\nDihedrals\n']
for ind, atoms in enumerate(new_dihedrals['atoms']):
    atoms = ' '.join([str(a) for a in atoms])
    newfile += [f" {ind + 1} {new_dihedrals['type'][ind]} {atoms}"]
newfile += ['\nImpropers\n']
for ind, atoms in enumerate(new_impropers['atoms']):
    atoms = ' '.join([str(a) for a in atoms])
    newfile += [f" {ind + 1} {new_impropers['type'][ind]} {atoms}"]

with open(output, 'w') as f:
    for line in newfile:
        f.write(line + '\n')
