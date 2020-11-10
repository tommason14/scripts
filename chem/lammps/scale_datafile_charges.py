#!/usr/bin/env python3
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='LAMMPS datafile- required')
parser.add_argument('-c',
                    '--charge',
                    help='Scaling factor. Default = 0.8',
                    type=float,
                    default=0.8)
parser.add_argument('-o',
                    '--output',
                    help='Filename of output. Default is scaled_charges.data',
                    default='scaled_charges.data')
args = parser.parse_args()

if args.input is None:
    parser.print_help()
    sys.exit()

with open(args.input) as f:
    contents = f.readlines()

new = []
found = False
for line in contents:
    if 'Atoms' in line:
        found = True
        new.append(line)
        continue
    if 'Bonds' in line:
        found = False
    if found and len(line.split()) > 2:
        line = line.split()
        line[3] = str(float(line[3]) * args.charge)
        line = ' '.join(line) + '\n'
    new.append(line)

with open(args.output, 'w') as f:
    for line in new:
        f.write(f'{line}')
