#!/usr/bin/env python3

"""
For each atom of a GAMESS log file, find the corresponding atom type in the
given monomer data file and average the charge across each type.
Then find the right line in the force field file.

IMPORTANT: USE THE SAME XYZ FILE TO RUN THE GEODESIC CALC AND
MAKE THE MONOMER DATA FILE- OR THIS WILL PROBABLY ALLOCATE 
CHARGES TO THE WRONG ATOMS
"""
import re
import argparse
import sys

def eof(file, percFile):
    # OPEN IN BYTES
    with open(file, "rb") as f:
        f.seek(0, 2)  # Seek @ EOF
        fsize = f.tell()  # Get size
        Dsize = int(percFile * fsize)
        f.seek(max(fsize - Dsize, 0), 0)  # Set pos @ last n chars lines
        lines = f.readlines()  # Read to end

    # RETURN DECODED LINES
    for i in range(len(lines)):
        try:
            lines[i] = lines[i].decode("utf-8")
        except:
            lines[i] = "CORRUPTLINE"
            print("eof function passed a corrupt line in file ", File)
        # FOR LETTER IN SYMBOL
    return lines

def read_file(file):
    with open(file, "r") as f:
        try:
            for line in f:
                yield line
        except UnicodeDecodeError:
            pass

parser = argparse.ArgumentParser()
parser.add_argument('--charges', '-c', help='Log file of GAMESS geodesic charge calculation')
parser.add_argument('--data',    '-d', help='Data file containing info for one molecule')
parser.add_argument('--forcefield', '-f', help='Force field file to be updated')
parser.add_argument('--output', '-o', help='Name of the modified force field (Optional)')
args = parser.parse_args()

if len(sys.argv) < 3:
    parser.print_help()
    sys.exit()

if not args.output:
    args.output = args.forcefield.replace('.ff', '_modified.ff')

def find_charges(log):
    """
    Searches log file for charges, returns nested list of dicts
    """
    charges = []
    for line in eof(log, 0.2):
        match = re.match('^\s*([A-z]+)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)', line)
        if match is not None:
            atom = match.group(1)
            charge = match.group(2)
            charges.append({'atom': atom, 'charge': charge})
    return charges

def find_atom_types(datafile):
    """
    Finds

    1    1   3  0.1919    0.071155   -0.039661    0.031933   # HC
    2    1   1  0.3835    0.066881    0.117109    1.111691   # C2
    3    1   1  0.2559    1.138454   -0.204681    1.843109   # C2
    4    1   3  0.1727    2.017804   -0.624621    1.346867   # HC
    5    1   4  0.1870   -0.839329    0.542180    1.548774   # HA

    and returns the last column
    """
    types = []
    found = False
    for line in read_file(datafile):
        if 'Atoms' in line:
            found = True
        if 'Bonds' in line:
            break
        if found:
            line = line.split()
            if len(line) == 9:
                types.append(line[-1])
    return types

def combine_charges_and_types(atoms, types):
    """
    Assume atom 1 of charges list is atom 1 of the types list,
    and performs an element-wise combination, where the type 
    becomes a['type']
    """
    # Assumes same ordering of atoms!!!
    for a, t in zip(atoms, types):
        a['type'] = t
    return atoms

def average_charge_per_type(atoms):
    """
    Now groups atoms by 'type' and then averages the 
    charge
    """
    types = {}
    for atom in atoms:
        if atom['type'] not in types:
            types[atom['type']] = [atom['charge']]
        else:
            types[atom['type']].append(atom['charge'])
    types = {
        k : list(map(float, v)) for k,v in types.items()
    }
    return {k: sum(v) / len(v) for k,v in types.items()}

def add_charges_to_ff(ff, ave_charge, new_ff):
    """
    Takes the mean charges and adds them to the force field
    file in the correct place.
    """
    newfile = []
    found_sec = False
    for line in read_file(ff):
        if 'ATOMS' in line:
            found_sec = True
        if 'BONDS' in line:
            found_sec = False
        if found_sec:
            line = line.split()
            if len(line) > 7:
                if line[0] in ave_charge.keys():
                    line[3] = ave_charge[line[0]]
                    line[7] = ' '.join(line[7:])
                else:
                    line[3] = float(line[3])
                    line[7] = ' '.join(line[7:])
                line = "{:<5} {:<5} {:>5} {:>8.4f} {} {:>8} {:^18} {}\n".format(*line)
            else:
                line = ' '.join(line + ['\n'])
        newfile.append(line)

    with open(new_ff, 'w') as f:
        for line in newfile:
            f.write(line)

def main():
    atoms = find_charges(args.charges)
    types = find_atom_types(args.data)
    atoms = combine_charges_and_types(atoms, types)
    ave_charge = average_charge_per_type(atoms)
    add_charges_to_ff(args.forcefield, ave_charge, new_ff = args.output)

main()
