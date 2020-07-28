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
            num, mass, _hash, label  = line.split()
            # keep num as string to lookup later
            atom_types[num] = label
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
