#!/usr/bin/env python3
"""
For each atom of a GAMESS log file, find the corresponding atom type in the
given monomer data file and average the charge across each type.
Then find the right line in the force field file.
"""
import argparse
import re
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


def find_charges(log):
    """
    Searches log file for charges, returns nested list of dicts
    """
    charges = []
    found = False
    for line in eof(log, 0.2):
        if 'NET CHARGES' in line:
            found = True
        if 'RMS DEVIATION' in line:
            found = False
        match = re.match(
            '^\s*([A-z]+)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)', line)
        if found and match is not None:
            atom = match.group(1)
            charge = match.group(2)
            charges.append({'atom': atom, 'charge': charge})
    return charges


def find_atom_types(xyz):
    """
    Returns the atom labels for each atom of an xyz
    """
    types = []
    for num, line in enumerate(read_file(xyz)):
        if num > 1:
            types.append(line.split()[0])
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
    types = {k: list(map(float, v)) for k, v in types.items()}
    return {k: sum(v) / len(v) for k, v in types.items()}


def average_charges_for_all_molecules(charges):
    """
    'charges' contains charges for each atom type in each given molecule.
    This function averages the charges across all molecules.
    """
    overall = {}
    for charge_dict in charges:
        for atomtype, charge in charge_dict.items():
            if atomtype not in overall:
                overall[atomtype] = [charge]
            else:
                overall[atomtype].append(charge)
    return {k: sum(v) / len(v) for k, v in overall.items()}


def add_charges_to_ff(ff, ave_charge, new_ff):
    """
    Takes the average charges and adds them to the force field
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
                line = "{:<5} {:<5} {:>5} {:>8.4f} {} {:>8} {:^18} {}\n".format(
                    *line)
            else:
                line = ' '.join(line + ['\n'])
        newfile.append(line)

    with open(new_ff, 'w') as f:
        for line in newfile:
            f.write(line)


def main():
    parser = argparse.ArgumentParser(description=(
        'Merges partial charges from GAMESS geodesic charge '
        'calculations into a force field file\n'
        'Pass in xyz files along with the corresponding GAMESS log '
        'files. Example usage: average_partial_charges_across_molecules.py '
        '-f gaff.ff amps.xyz amps_spec.log gdma.xyz gdma_spec.log -o gaff.amps.gdma.ff'
    ))

    parser.add_argument('-f',
                        '--forcefield',
                        help='Force field parameter file')
    parser.add_argument('-o',
                        '--output',
                        help='Filename of modified force field')
    parser.add_argument(
        'files',
        nargs='+',
        help=
        ('xyz1 log1 [xyz2 log2 ...], where xyz_i is a labelled xyz of one molecule '
         'and log_i is the GEODESIC calculation of that same structure'))
    args = parser.parse_args()

    if any(x is None for x in (args.files, args.forcefield,
                               args.output)) or len(args.files) % 2 != 0:
        parser.print_help()
        sys.exit()

    xyzs = args.files[::2]
    logs = args.files[1::2]

    ave_per_molecule = []
    for xyz, log in zip(xyzs, logs):
        charges = find_charges(log)
        types = find_atom_types(xyz)
        combined = combine_charges_and_types(charges, types)
        ave_charge = average_charge_per_type(combined)
        ave_per_molecule.append(ave_charge)

    ave_charge = average_charges_for_all_molecules(ave_per_molecule)
    add_charges_to_ff(args.forcefield, ave_charge, new_ff=args.output)


main()
