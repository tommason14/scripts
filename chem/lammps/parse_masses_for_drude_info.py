#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    sys.exit("Syntax: print_masses_for_drude_info.py datafile")

"""
Using labels from Masses section:
1      12.0100    # C DC
2      12.0100    # C DC
3       1.0080    # H
4       1.0080    # H
5       1.0080    # H
6       0.4       # C DP
"""
atoms = []
cores = []
drudes = []
found = False

with open(sys.argv[1]) as f:
    for line in f:
        if "Masses" in line:
            found = True
            continue
        if len(atoms) > 0 and re.search("^\s*[A-Z]", line):
            break
        if found and '#' in line:
            line = line.split()
            if line[-1] == 'DP':
                drudes.append(line[0])
            elif line[-1] == 'DC':
                atoms.append(line[0])
                cores.append(line[0])
            else: # hydrogens
                atoms.append(line[0])

atoms = ' '.join(atoms)
cores = ' '.join(cores)
drudes = ' '.join(drudes)
print(f'atoms = {atoms}')
print(f'cores = {cores}')
print(f'drudes = {drudes}')
