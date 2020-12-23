#!/usr/bin/env python3

import os
import sys
import re
"""
print tip4p information to the screen
- oxygen/hydrogen atoms
- oh bond/angle number
- sodium (NAM)
- chloride (CL)
"""


def check_args():
    if not len(sys.argv) == 2 or sys.argv[1] == '-h':
        print('Syntax: print_tip4p_info.py datafile')
        sys.exit()


def attributes(datafile):
    """
    Returns the atom types, H-O bond ID and H-O-H angle ID for water.
    Assumes that the IDs are HW and OW.
    Also returns sodium (NAM) and chloride (CL) atom numbers.
    """
    otype = ''
    htype = ''
    bond = ''
    angle = ''
    na = ''
    cl = ''
    in_coeffs_section = False
    with open(datafile) as f:
        for line in f:
            line = line.upper().strip()
            if 'ZLO' in line:
                in_coeffs_section = True
            if 'ATOMS' in line and in_coeffs_section:
                break
            if line.endswith('# NAM'):
                na = line.split()[0]
            if line.endswith('# CL'):
                cl = line.split()[0]
            if line.endswith('# OW'):
                otype = line.split()[0]
            if line.endswith('# HW'):
                htype = line.split()[0]
            if line.endswith('# HW-OW') or line.endswith('# OW-HW'):
                bond = line.split()[0]
            if line.endswith('# HW-OW-HW'):
                angle = line.split()[0]
    return otype, htype, bond, angle, na, cl


def main():
    check_args()
    otype, htype, ohbond, ohangle, na, cl = attributes(sys.argv[1])
    print(f'otype = {otype}')
    print(f'htype = {htype}')
    print(f'ohbond = {ohbond}')
    print(f'ohangle = {ohangle}')
    print(f'na = {na}')
    print(f'cl = {cl}')


if __name__ == "__main__":
    main()
