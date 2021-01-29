#!/usr/bin/env python3

from autochem import Molecule, write_xyz
from glob import glob
import os
import sys
from tqdm import tqdm

base = os.path.basename(__file__)

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    print('Extract all ions of a particular name from all xyz files in the current directory.')
    print(f'i.e. `{base} dhp` pulls out all dhp anions from all xyz files and place them in a subdirectory')
    print(f'Syntax: {base} molecule')
    sys.exit(1)

if not os.path.isdir('extracted'):
    os.mkdir('extracted')

for f in tqdm(glob('*xyz')):
    mol = Molecule(f)
    newname = os.path.join(os.getcwd(), 'extracted', f)
    extracted = [atom for atom in mol.coords if sys.argv[1] in atom.fragment]
    write_xyz(extracted, newname) 
