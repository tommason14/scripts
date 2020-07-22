#!/usr/bin/env python3

import re
import sys

if len(sys.argv) < 2:
    sys.exit('Syntax: xyz_from_lammps_datafile.py lammpsfile [output.xyz]')

lammps = sys.argv[1]
if len(sys.argv) == 2:
    output = lammps.rsplit('.')[0] + '.xyz'
else:
    output = sys.argv[2]

masses = {
    'H': 1.008,
    'C': 12.01,
    'C': 12.011,
    'N': 14.01,
    'N': 14.007,
    'O': 16,
    'O': 15.999,
    'S': 32.06,
    'Cl': 35.45,
    'P': 30.974,
}

coords = []
atom_types = {}
found_m = False
found_a = False
with open(lammps) as f:
    for line in f:
        if 'Masses' in line:
            found_m = True
            continue
        if found_m and re.search('^\s*[A-Z]', line):
            found_m = False
        if found_m and not re.search('^\s*$', line):
            try:
                num, mass = line.split()
            except ValueError:
                num, mass, *_ = line.split()
            # keep num as string to lookup later
            mass = float(mass)
            for element, reference_mass in masses.items():
                if mass == reference_mass:
                    atom_types[num] = element
        if 'Atoms' in line:
            found_a = True
            continue
        if found_a and re.search('^\s*[A-Z]', line):
            found_a = False
        if found_a and not re.search('^\s*$', line):
            line = line.split()
            elem, xyz = line[2], line[4:7]
            xyz = [float(i) for i in xyz] # deal with the lammps 1e-01 format
            elem = atom_types[elem]
            coords.append([elem] + xyz)

with open(output, 'w') as f:
    f.write(f'{len(coords)}\n\n')
    for coord in coords:
        f.write('{:4} {:>15.10f} {:>15.10f} {:>15.10f}\n'.format(*coord))
