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

if len(sys.argv) < 2 or sys.argv[1] == "-h":
    sys.exit("Syntax: keep_largest_fragment.py <lammps_datafile>")
datafile = sys.argv[1]

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

# Molecular properties
atoms = {"id": [], "mol": [], "type": [], "rest": []}

with open(datafile) as f:
    for line in f:
        if is_data(line):
            if "Atoms" in line:
                found_atoms = True
                continue
            if found_atoms and re.search("^\s*[A-Z]", line):
                break
            # Coeffs
            if found_atoms:
                if len(line.split()) > 2:
                    atom_id, mol, atom_type, *rest = line.split()
                    atoms["id"].append(atom_id)
                    atoms["mol"].append(int(mol))
                    atoms["type"].append(atom_type)
                    atoms["rest"].append(rest)
# now select fragment with the most atoms
molcounts = Counter(atoms["mol"])
mol_to_keep = molcounts.most_common(1)[0][0]
print("Add the following lines to the input file of your next run:")
print("# remove small molecules/unreacted monomers")
print(f"group keep molecule {mol_to_keep}")
print("group DEL subtract all keep")
print("delete_atoms group DEL")
print("reset_atom_ids")
print("# If LAMMPS version is earlier than 29Oct20, use reset_ids instead")
print("# write_data removed_monomers.data")
