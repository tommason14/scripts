#!/usr/bin/env python3
import sys
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--datafile', default='pdata.lmp', required=True, help='Lammps datafile')
parser.add_argument('-x', '--exclude', nargs='+', help='Atom types to exclude i.e. HO for hydroxyl groups. If not given, all bonds involving a hydrogen will be fixed.')

args = parser.parse_args()

""" 
Checking Bond Coeffs section for # HC-CA or # CA-HC etc...
"""
shake = []
found = False
with open(args.datafile) as f:
    for line in f:
        if "Bond Coeffs" in line:
            found = True
            continue
        if "Angle Coeffs" in line:
            break
        if found and len(line.split()) > 1:
            line = line.split()
            if line[-1].startswith('H') or '-H' in line[-1]:
                if args.exclude and any(a in line[-1] for a in args.exclude):
                    pass
                else:
                    shake.append(line[0])

print(" ".join(shake))
