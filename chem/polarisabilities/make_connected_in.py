#!/usr/bin/env python3
from autochem import Molecule
import sys

if len(sys.argv) not in [2,3] or '-h' in sys.argv:
    print(f'Syntax: {sys.argv[0]} xyzfile [output]')
    sys.exit(1)

try:
    output = sys.argv[2]
except IndexError:
    output = 'connected.in'

mol = Molecule(sys.argv[1])
bonds = []
for atom in mol.coords:
    for a in atom.connected_atoms:
        if f'bond {a.index} {atom.index}' not in bonds:
            bonds.append(f'bond {atom.index} {a.index}')
with open(output, 'w') as f:
    for b in bonds:
        f.write(f'{b}\n')
