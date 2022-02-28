#!/usr/bin/env python3
"""
Change molecule ID in a lammps data file, generated from packmol.
By default, all atoms of a file generated using lmp_gaff.py will have 
an ID of 1. Reads pack.inp to generate the new molecule ID list.
Assumes that the simulation was created from xyz files.
"""
import os
import subprocess as sp
import sys
import re


def read_file(file):
    with open(file, "r") as f:
        try:
            for line in f:
                yield line
        except UnicodeDecodeError:
            pass


if not 2 <= len(sys.argv) <= 3:
    sys.exit(f"Syntax: {os.path.basename(__file__)} lammps_file [output_file]")

lammps = sys.argv[1]
if len(sys.argv) == 2:
    output = lammps
else:
    output = sys.argv[2]

packfile = os.path.splitext(lammps)[0] + ".inp"
# map structures to numbers of each
structs = (
    sp.check_output(
        "grep '^\s*structure.*xyz' " + packfile + " | awk '{print $NF}'", shell=True
    )
    .decode("utf-8")
    .strip()
    .split("\n")
)
numbers = (
    sp.check_output(
        "grep '^\s*number' " + packfile + " | awk '{print $NF}'", shell=True
    )
    .decode("utf-8")
    .strip()
    .split("\n")
)
if len(structs) != len(numbers):
    sys.exit(
        'Error in pack.inp: Make sure all xyz files have a "number" line in the structure block.\n'
        "Cannot assign correct molecule IDs."
    )

mols = {s: {"number": int(n)} for s, n in zip(structs, numbers)}
# lengths of each molecule- first line of xyz
for mol, data in mols.items():
    data["length"] = int(
        sp.check_output(f"head -1 {mol}", shell=True).decode("utf-8").strip()
    )

# packmol retains the order in which molecules are added in pack.inp
# so loop over the dict and generate new id
mol_id = []
mol_count = 1
for s in structs:
    for _ in range(mols[s]["number"]):
        for _ in range(mols[s]["length"]):
            mol_id.append(mol_count)
        mol_count += 1
# now change data file

b4_atoms = []
atoms = []
after_atoms = []
b4_atoms_bool = True
found_atoms_bool = False
after_atoms_bool = False
for line in read_file(lammps):
    if re.search("\s*Atoms", line):
        b4_atoms_bool = False
        found_atoms_bool = True
    if re.search("^\s*Bonds", line):
        found_atoms_bool = False
        after_atoms_bool = True
    if b4_atoms_bool:
        b4_atoms.append(line)
    if found_atoms_bool:
        atoms.append(line)
    if after_atoms_bool:
        after_atoms.append(line)

newatoms = []
for line, molid in zip(atoms[2:-1], mol_id):
    line = line.split()
    line[1] = str(molid)
    newatoms.append(" ".join(line) + "\n")
atoms = atoms[:2] + newatoms + [atoms[-1]]

with open(output, "w") as f:
    f.writelines(b4_atoms)
    f.writelines(atoms)
    f.writelines(after_atoms)
