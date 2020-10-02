#!/usr/bin/env python3

import argparse
import os
import sys
"""
Converting xyz created by packmol to pdb format and changing residues
to the name of each molecule.
"""

parser = argparse.ArgumentParser(
    description=('Convert xyz created by packmol into pdb. '
                 'If no pdb given, uses xyz and replaces extension with pdb'))
parser.add_argument('-x', '--xyz', help='xyz file to convert')
parser.add_argument('-i',
                    '--input',
                    help='packmol input file i.e. pack.inp',
                    default='pack.inp')
parser.add_argument('-o', '--output', help='name of pdb file to create')
args = parser.parse_args()

if len(sys.argv) == 1 or '-h' in sys.argv:
    parser.print_help()
    sys.exit()

if not args.output:
    args.output = args.xyz.replace('xyz', 'pdb')

os.system(f"obabel -ixyz {args.xyz} -opdb > tmp")

mols = {}
with open(args.input) as f:
    for line in f:
        if line.startswith("structure"):
            if '_pack' in line:
                # fftool- from ch_pack.xyz, take ch
                name = line.split()[1].split("_")[0]
            else:
                # otherwise lose extension and include the first two
                # characters
                name = line.split()[1].rsplit('.')[0][:3]
        if "number" in line:
            num = int(line.split()[1])
        if "end structure" in line:
            mols[name] = num

residues = []
for name, num in mols.items():
    for _ in range(num):
        residues.append(name)
residues = {ind: name for ind, name in enumerate(residues, 1)}

with open("tmp") as f:
    contents = f.readlines()

newfile = []
for line in contents:
    if "HETATM" in line:
        tmp = line.split()
        tmp[3] = residues[int(tmp[4])]
        # pdb formatting retained
        tmp = "{:<8}{:>3}  {:<4}{:<8}{:<7}{:<8}{:<8}{:<8}{:<6}{:<15}{}\n".format(
            *tmp)
        newfile.append(tmp)
    else:
        newfile.append(line)

with open("simbox.pdb", "w") as f:
    for line in newfile:
        f.write(line)

os.system("rm tmp")
