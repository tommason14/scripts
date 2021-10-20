#!/usr/bin/env python3

"""
File: add_atomic_symbols_to_lammps_datafile.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Just takes the Masses section of a datafile
and adds a # atomname to each atom type.

i.e. 
Masses

1 1.008
2 16

becomes
1 1.008 # H
2 16    # O

so that ovito knows which colour to provide
"""

import sys
import re

masses = {
    '1.0075': 'H', 
    '1.00797': 'H', 
    '1.008':  'H', 
    '1.0080':  'H', 
    '12.01':  'C', 
    '12.0100': 'C',
    '12.011': 'C', 
    '12.0112': 'C', 
    '14.01':  'N', 
    '14.0100':  'N', 
    '14.007': 'N', 
    '16':     'O', 
    '16.0000':     'O', 
    '15.9990': 'O', 
    '15.999': 'O', 
    '22.99': 'Na',
    '22.9900': 'Na',
    '30.974': 'P',
    '32.06':  'S', 
    '32.0600':  'S', 
    '32.065': 'S', 
    '32.066': 'S', 
    '35.45':  'Cl',
    '35.4500':  'Cl',
    '55.847': 'Fe'
}

filename = sys.argv[1]

new = []

found = False
with open(filename, 'r') as f:
    for line in f:
        if 'Masses' in line:
            found = True
            new.append(line)
            continue
        if re.search('^\s*[A-Z]', line):
            found = False
        if found and not re.search('^\s*$', line):
            line = line.split()
            changed = '{} {} # {}\n'.format(*line, masses[line[-1]])
            new.append(changed)
        else:
            new.append(line)

with open(filename, 'w') as f:
    f.writelines(new)
