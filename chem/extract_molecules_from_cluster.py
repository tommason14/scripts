#!/usr/bin/env python3

from autochem import Molecule, write_xyz
from glob import glob
import os
import sys

base = os.path.basename(__file__)

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print('Extract all ions of a particular name from all xyz files in the parent directory.')
    print(f'i.e. `{base} dhp` pulls out all dhp anions from xyz in the directory above.')
    print(f'Syntax: {base} molecule')
    sys.exit(1)

for f in glob('../*xyz'):
    newname = f.replace('../', '')
    mol = Molecule(f)
    extracted = [atom for atom in mol.coords if sys.argv[1] in atom.fragment]
    write_xyz(extracted, newname) 
