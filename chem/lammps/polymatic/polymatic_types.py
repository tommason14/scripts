#!/usr/bin/env python3

import argparse
import sys
import re

parser = argparse.ArgumentParser()
parser.add_argument("--inp", "-i", help="Lammps data file name")
parser.add_argument("--out", "-o", help="Polymatic types file name")
args = parser.parse_args()

if len(sys.argv) < 3:
    parser.print_help()
    sys.exit()

found_masses = False
found_bonds = False
found_angles = False
found_dihedrals = False
found_impropers = False

is_data = lambda line: len(line.split()) > 0

atoms = {}
bonds = {}
angles = {}
dihedrals = {}
impropers = {}

with open(args.inp) as inpfile:
    for line in inpfile:
        if is_data(line):
            if "Masses" in line:
                found_masses = True
                continue
            if found_masses and "Bond Coeffs" in line or 'Pair Coeffs' in line:
                found_masses = False
                found_bonds = True
                continue
            if found_bonds and "Angle Coeffs" in line:
                found_bonds = False
                found_angles = True
                continue
            if found_angles and "Dihedral Coeffs" in line:
                found_angles = False
                found_dihedrals = True
                continue
            if found_dihedrals and re.search('^\s*[A-Z]', line):
                found_dihedrals = False
                continue
            if 'Improper Coeffs' in line:
                found_impropers = True
                continue
            if found_impropers and re.search('^\s*[A-Z]', line):
                found_impropers = False
                continue
            if found_masses:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    atoms[count] = name
            if found_bonds:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    bonds[count] = name
            if found_angles:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    angles[count] = name
            if found_dihedrals:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    dihedrals[count] = name
            if found_impropers:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    impropers[count] = name

alldata = {
    "atom types": atoms,
    "bond types": bonds,
    "angle types": angles,
    "dihedral types": dihedrals,
    "improper types": impropers,
}

# if no impropers...
alldata = {k: v for k,v in alldata.items() if len(v) > 0}

with open(args.out, "w") as out:
    for string, data in alldata.items():
        out.write(f"{string}\n")
        for count, item in data.items():
            out.write(f"{count:<5} {item.replace('-', ',')}\n")
        out.write("#\n")
