#!/usr/bin/env python3

"""
File: add_atom_types_from_datafile.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: takes atom types defined in the masses section of one LAMMPS datafile and adds them to
another.
"""

import argparse
import sys
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--from-file', help='Datafile with desired atom types', required=True)
parser.add_argument('-t', '--to', help='Datafile with no atom types that you wish to modify', required=True)

args = parser.parse_args()


types = []
found = False
with open(args.from_file) as f:
    for line in f:
        if 'Masses' in line:
            found = True
            continue
        if re.search('^\s*[A-Z]', line) and found:
            break
        if found:
            line = line.split()
            if len(line) > 2:
                types.append(line[-1])

with open(args.to) as f:
    new = f.readlines()

found_masses = False
for ind, line in enumerate(new):
    if 'Masses' in line:
        found_masses = True
        continue
    if re.search('^\s*[A-Z]', line) and found_masses:
        break
    if found_masses:
        line = line.split()
        if len(line) == 2:
            new[ind] = ' '.join(line + ['#', types.pop(0), '\n'])

with open(args.to, 'w') as f:
    for line in new:
        f.write(line)
