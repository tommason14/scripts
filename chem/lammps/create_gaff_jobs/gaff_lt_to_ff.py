#!/usr/bin/env python3

from autochem import read_file
from pprint import pprint
from datetime import date
import re

atoms = []
bonds = []
angles = []
dihs = []
imps = []

def find_masses(f, atoms):
    """
    Parses 
    @atom:pc 30.97   # Sp2 P in non-pure aromatic systems
    for the atomic symbol, mass and comment (after #)
    Coincidentally, the atom labels are the same as atom groups, as 
    seen in the pair coeffs section, so the symbol is duplicated to
    provide labels and types
    """
    found = False
    for line in f:
        if 'Data Masses' in line:
            found = True
        if found and '@' in line:
            sym, mass, *comment = line.split()
            comment = ' '.join(comment)
            sym = sym.split(':')[1].upper()
            atoms.append([sym, sym, mass])
        if 'end of masses' in line:
            break

def pair_coeffs(f, atoms):
    """
    Parses
    pair_coeff @atom:h1 @atom:h1 lj/charmm/coul/long 0.0157 2.47135304412   #
    Veenstra et al JCC,8,(1992),963
    for the epsilon (0.0157) and sigma (2.47135304412) terms.
    """
    for line in f:
        if 'pair_coeff' in line and '@' in line:
            line = line.split()
            sym, eps, sigma, comment = line[1], line[4], line[5], line[6:]
            sym = sym.split(':')[1].upper()
            comment = ' '.join(comment)
            for a in atoms:
                if a[0] == sym:
                    a += [eps, sigma, comment]
        if 'end of pair_coeffs' in line:
            break

def bond_coeffs(f, bonds):
    """
    Parses 
    bond_coeff @bond:ow-hw harmonic 553.0 0.9572   # TIP3P_Water 1
    for the bonding atoms, force constant, equil length and comment
    """
    for line in f:
        if 'bond_coeff' in line and '@' in line:
            line = line.split()
            bond, k, length, comment = line[1], line[3], line[4], line[5:]
            atoms = bond.split(':')[1].upper()
            comment = ' '.join(comment)
            bonds.append([atoms, k, length, comment])
        if 'end of bond_coeffs' in line:
            break

def angle_coeffs(f, angles):
    """
    Parses 
    angle_coeff @angle:br-c1-br harmonic 57.76 180.00   # Guess 0
    for the atoms, force constant, equil angle and comment
    """
    for line in f:
        if 'angle_coeff' in line and '@' in line:
            line = line.split()
            atoms, k, angle, comment = line[1], line[3], line[4], line[5:]
            atoms = atoms.split(':')[1].upper()
            comment = ' '.join(comment)
            angles.append([atoms, k, angle, comment])
        if 'end of angle_coeffs' in line:
            break

def dihedral_coeffs(f, dihs):
    """
    Parses 
    dihedral_coeff @dihedral:X-c-cg-X fourier 1 0.0 2 180.0    # same as X-c-c1-X
    for the atoms involved, the arguments to the dihedral/fourier function and
    comment.
    There are a variable number of arguments to collect here.
    """
    for line in f:
        if 'dihedral_coeff' in line and '@' in line:
            line = line.split()
            atoms = line[1].split(':')[1].upper()
            for idx, val in enumerate(line):
                if '#' in val:
                    splitter = idx
                    break
            args = line[3:splitter]
            comment = ' '.join(line[splitter:])
            dihs.append([atoms] + args + [comment])
        if 'end of dihedral_coeffs' in line:
            break

def improper_coeffs(f, imps):
    """
    Parses
    improper_coeff @improper:X-o-c-o cvff 1.1 -1 2   # JCC,7,(1986),230
    for the atoms, args and comment
    """
    for line in f:
        if 'improper_coeff' in line and '@' in line:
            line = line.split()
            atoms = line[1].split(':')[1].upper()
            args = line[3:6]
            comment = ' '.join(line[6:])
            imps.append([atoms] + args + [comment])
        if 'end of improper_coeffs' in line:
            break

def write_file(atoms, bonds, angles, dihs, imps):
    """
    Write a readable force field parameter file
    """
    new = []
    new.append(f'# GAFF.FF, converted on {date.today().strftime("%d %B %Y")}')
    new.append('')
    new.append('ATOMS')
    new.append('')
    new.append('# charge needs filling in ...')
    for a in atoms:
        new.append("{:<5} {:<5} {:>5}  0.00000  lj {:>8} {:^18} {}".format(*a))
    new.append('')
    new.append('BONDS')
    new.append('')
    for b in bonds:
        bond_type = 'cons' if '-H' in b[0] or re.search('^H', b[0]) else 'harm'
        b.insert(1, bond_type)
        b = b[0].split('-') + b[1:]
        new.append("{:<4} {:<4} {:>7} {:>7} {:>7} {}".format(*b))
    new.append('')
    new.append('ANGLES')
    new.append('')
    for a in angles:
        a.insert(1, 'harm')
        a = a[0].split('-') + a[1:]
        new.append("{:<4} {:<4} {:<4} {:>8} {:>7} {:>7} {}".format(*a))
    new.append('')
    new.append('DIHEDRALS')
    new.append('')
    for d in dihs:
        args = '  '.join(d[1:-1])
        d = d[0].split('-') + ['fourier', args, d[-1]]
        new.append("{:<4} {:<4} {:<4} {:<4} {} {}   {}".format(*d))
    new.append('')
    new.append('IMPROPERS')
    new.append('')
    for i in imps:
        args = "  ".join(i[1:-1])
        i = i[0].split('-') + ['cvff', args, i[-1]] 
        new.append("{:<4} {:<4} {:<4} {:<4} {:>6} {}   {}".format(*i))
    with open('gaff.ff', 'w') as n:
        for line in new:
            n.write(line + '\n')

def main():
    file = 'gaff.lt'
    file = read_file(file)
    find_masses(file, atoms)
    pair_coeffs(file, atoms)  
    bond_coeffs(file, bonds)
    angle_coeffs(file, angles)
    dihedral_coeffs(file, dihs)
    improper_coeffs(file, imps)
    write_file(atoms, bonds, angles, dihs, imps)


if __name__ == "__main__":
    main()
