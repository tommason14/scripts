#!/usr/bin/env python3

from autochem import (Atom, Molecule,
read_file, get_files, write_csv_from_nested)
import os
import sys
import re


def file_is_gamess(file):
    """ Check first line of file for 'rungms' string """
    with open(file, 'r') as f:
        return 'rungms' in f.readline()

atom_regex = '^\s[A-Za-z]{1,2}\s*[0-9]*.[0-9]*(\s*-?[0-9]*.[0-9]*){3}$'
charge_regex = '^\s[A-Za-z]{1,2}(\s*-?[0-9]*.[0-9]*){2}$'

results = []

files = get_files('.', ['log', 'out'])

for logfile in files:
    if file_is_gamess(logfile):
        path, filename = os.path.split(logfile)
        print(path)
        inpfile = logfile[:-3] + 'inp'

        res = []
        assigned = []

        for line in read_file(inpfile):
            if re.search(atom_regex, line):
                sym, atnum, x, y, z = line.split()
                x, y, z = map(float, (x, y, z))
                res.append([path, Atom(sym, coords = (x, y, z))]) # new key for each coord
        found = False
        counter = 0
        for line in read_file(logfile):
            if 'NET CHARGES:' in line:
                found = True
            if 'RMS DEVIATION' in line:
                break
            if found: 
                if re.search(charge_regex, line):
                    res[counter].append(float(line.split()[1]))
                    counter += 1
        coords = [atom[1] for atom in res]
        mol = Molecule(atoms = coords)
        mol.separate()
        # bad practice, should use for i, j in zip(mol.coords, res):
        for atom, r in zip(mol.coords, res):
            path, _, geodesic_charge = r
            results.append([path, atom.index, atom.symbol, geodesic_charge, atom.x, atom.y, atom.z, f"{mol.fragments[atom.mol]['name']}_{atom.mol}"])

write_csv_from_nested(results, col_names = ('Path', 'Index', 'Element', 'Geodesic', 'Rx', 'Ry', 'Rz', 'Fragment'), filename='charges.csv')
