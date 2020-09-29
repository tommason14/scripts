#!/usr/bin/env python3

import argparse
import re
"""
Modify a lammps input file.
Assumes that the file has 
- a line containing the lj/cut/tip4p/... pair style,
- a line that sets a 'tip4p' group,
- a line that sets a shake fix for that tip4p group
- a line setting elements for a dump command
"""

parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-d',
                      '--datafile',
                      help='LAMMPS datafile to take info from',
                      required=True)
required.add_argument('-i',
                      '--input',
                      help='Input file to place info into',
                      required=True)
optional.add_argument(
    '-o',
    '--output',
    help='Output filename. If not given, file is printed to stdout')
args = parser.parse_args()


def traj_elements():
    """
    Extracting elements from
    1      12.0100    # C
    2      12.0100    # C3
    3       1.0080    # HC
    4       1.0080    # H1
    5       1.0080    # HO
    6       1.0080    # HN
    """
    # need to allow Na and N, while trimming off numbers
    # also correct for Polymatic linkers (LC -> carbon)
    # do this by checking masses and not the force field types
    masses = []
    found = False
    with open(args.datafile) as f:
        for line in f:
            if 'Masses' in line:
                found = True
                continue
            if len(masses) > 0 and re.search('^\s*[A-Z]', line):
                break
            if found and len(line.split()) > 1:  # element or #element
                masses.append(float(line.split()[1]))
    ptable = {}
    #atomic weights from: http://www.ciaaw.org/atomic-weights.htm
    ptable['Xx'] = 0.00000
    ptable['H'] = 1.00798
    ptable['He'] = 4.00260
    ptable['Li'] = 6.96750
    ptable['Be'] = 9.01218
    ptable['B'] = 10.81350
    ptable['C'] = 12.01060
    ptable['N'] = 14.00685
    ptable['O'] = 15.99940
    ptable['F'] = 18.99840
    ptable['Ne'] = 20.17970
    ptable['Na'] = 22.98977
    ptable['Mg'] = 24.30550
    ptable['Al'] = 26.98154
    ptable['Si'] = 28.08500
    ptable['P'] = 30.97376
    ptable['S'] = 32.06750
    ptable['Cl'] = 35.45150
    ptable['Ar'] = 39.94800
    ptable['K'] = 39.09830
    ptable['Ca'] = 40.07800
    ptable['Sc'] = 44.95591
    ptable['Ti'] = 47.86700
    ptable['V'] = 50.94150
    ptable['Cr'] = 51.99610
    ptable['Mn'] = 54.93804
    ptable['Fe'] = 55.84500
    ptable['Co'] = 58.93319
    ptable['Ni'] = 58.69340
    ptable['Cu'] = 63.54600
    ptable['Zn'] = 65.38000
    ptable['Ga'] = 69.72300
    ptable['Ge'] = 72.63000
    ptable['As'] = 74.92159
    ptable['Se'] = 78.97100
    ptable['Br'] = 79.90400
    ptable['Kr'] = 83.79800
    ptable['Rb'] = 85.46780
    ptable['Sr'] = 87.62000
    ptable['Y'] = 88.90584
    ptable['Zr'] = 91.22400
    ptable['Nb'] = 92.90637
    ptable['Mo'] = 95.95000
    ptable['Tc'] = 98.00000
    ptable['Ru'] = 101.07000
    ptable['Rh'] = 102.90550
    ptable['Pd'] = 106.42000
    ptable['Ag'] = 107.86820
    ptable['Cd'] = 112.41400
    ptable['In'] = 114.81800
    ptable['Sn'] = 118.71000
    ptable['Sb'] = 121.76000
    ptable['Te'] = 127.60000
    ptable['I'] = 126.90447
    ptable['Xe'] = 131.29300
    ptable['Cs'] = 132.90545
    ptable['Ba'] = 137.32700
    ptable['La'] = 138.90547
    ptable['Ce'] = 140.11600
    ptable['Pr'] = 140.90766
    ptable['Nd'] = 144.24200
    ptable['Pm'] = 145.00000
    ptable['Sm'] = 150.36000
    ptable['Eu'] = 151.96400
    ptable['Gd'] = 157.25000
    ptable['Tb'] = 158.92535
    ptable['Dy'] = 162.50000
    ptable['Ho'] = 164.93033
    ptable['Er'] = 167.25900
    ptable['Tm'] = 168.93422
    ptable['Yb'] = 173.04500
    ptable['Lu'] = 174.96680
    ptable['Hf'] = 178.49000
    ptable['Ta'] = 180.94788
    ptable['W'] = 183.84000
    ptable['Re'] = 186.20700
    ptable['Os'] = 190.23000
    ptable['Ir'] = 192.21700
    ptable['Pt'] = 195.08400
    ptable['Au'] = 196.96657
    ptable['Hg'] = 200.59200
    ptable['Tl'] = 204.38350
    ptable['Pb'] = 207.20000
    ptable['Bi'] = 208.98040
    ptable['Po'] = 209.00000
    ptable['At'] = 210.00000
    ptable['Rn'] = 222.00000
    ptable['Fr'] = 223.00000
    ptable['Ra'] = 226.00000
    ptable['Ac'] = 227.00000
    ptable['Th'] = 232.03770
    ptable['Pa'] = 231.03588
    ptable['U'] = 238.02891
    ptable['Np'] = 237.00000
    ptable['Pu'] = 244.00000
    ptable['Am'] = 243.00000
    ptable['Cm'] = 247.00000
    ptable['Bk'] = 247.00000
    ptable['Cf'] = 251.00000
    ptable['Es'] = 252.00000
    ptable['Fm'] = 257.00000
    ptable['Md'] = 258.00000
    ptable['No'] = 259.00000
    ptable['Lr'] = 266.00000
    ptable['Rf'] = 267.00000
    ptable['Db'] = 268.00000
    ptable['Sg'] = 269.00000
    ptable['Bh'] = 270.00000
    ptable['Hs'] = 277.00000
    ptable['Mt'] = 278.00000
    ptable['Ds'] = 281.00000
    ptable['Rg'] = 282.00000
    ptable['Cn'] = 285.00000
    ptable['Nh'] = 286.00000
    ptable['Fl'] = 289.00000
    ptable['Mc'] = 290.00000
    ptable['Lv'] = 293.00000
    ptable['Ts'] = 294.00000
    ptable['Og'] = 294.00000

    elements = []
    for mass in masses:
        for el, reference in ptable.items():
            # int rounds down, so int(15.999) = 15 nort 16, fix by adding 0.5
            if int(mass + 0.5) == int(reference + 0.5):
                elements.append(el)
    return ' '.join(elements)


def water_attributes():
    """
    Returns the atom types, H-O bond ID and H-O-H angle ID for water.
    Assumes that the IDs are HW and OW
    """
    otype = ''
    htype = ''
    bond = ''
    angle = ''
    with open(args.datafile) as f:
        for line in f:
            line = line.strip()
            if 'Atoms' in line:
                break
            if line.endswith('# OW'):
                otype = line.split()[0]
            if line.endswith('# HW'):
                htype = line.split()[0]
            if line.endswith('# HW-OW') or line.endswith('# OW-HW'):
                bond = line.split()[0]
            if line.endswith('# HW-OW-HW'):
                angle = line.split()[0]
    return otype, htype, bond, angle


def modify_input(vals):
    with open(args.input) as f:
        inp = f.readlines()
    newinp = []
    for line in inp:
        if 'lj/cut/tip4p' in line:
            tmp = line.split()
            # find element with the phrase, and only modify the 4 elements after that
            ID = 0
            for ind, el in enumerate(tmp):
                if 'lj/cut/tip4p' in el:
                    ID = ind
                    break
            newvals = [
                vals['otype'], vals['htype'], vals['wbond'], vals['wangle']
            ]
            tmp[ID + 1:ID + 5] = newvals
            newinp.append(' '.join(tmp))
        elif 'group tip4p' in line:
            newinp.append(f'group tip4p type {vals["otype"]} {vals["htype"]}')
        elif all(x in line for x in ('fix', 'tip4p', 'shake')):
            # assuming fix tip4p shake ... b bondtype a angletype
            tmp = line.split()
            tmp[-3] = vals['wbond']
            tmp[-1] = vals['wangle']
            newinp.append(' '.join(tmp))
        elif 'dump_modify' in line and 'element' in line:
            tmp = line.split()
            # replace after element
            ID = 0
            for ind, el in enumerate(tmp):
                if 'element' in el:
                    ID = ind
            tmp = tmp[:ID + 1] + [vals['elements']]
            newinp.append(' '.join(tmp))
        else:
            newinp.append(line.strip())
    if args.output:
        with open(args.output, 'w') as f:
            for line in newinp:
                f.write(f'{line}\n')
    else:
        for line in newinp:
            print(line)


def main():
    elements = traj_elements()
    otype, htype, wbond, wangle = water_attributes()
    modify_input(locals())


if __name__ == "__main__":
    main()
