#!/usr/bin/env python3
from autochem import Molecule
import sys

if len(sys.argv) not in [2,3] or '-h' in sys.argv:
    print(f'Syntax: {sys.argv[0]} xyzfile')
    sys.exit(1)

try:
    output = sys.argv[2]
except IndexError:
    output = 'connected.in'

mol = Molecule(sys.argv[1])
bonds = []
for atom in mol.coords:
    for a in atom.connected_atoms:
        if f'{a.symbol}{a.index} - {atom.symbol}{atom.index}' not in bonds:
            bonds.append(f'{atom.symbol}{atom.index} - {a.symbol}{a.index}')
for b in bonds:
    print(b)
