#!/usr/bin/env python3
import argparse
import os
import sys


def writepackmol(mols, files, box):
    lines = ['tolerance 2.0', 'filetype xyz', 'output pack.xyz', '']

    for n, xyz in zip(mols, files):
        lines += [f'structure {xyz}']
        lines += [f'  number {n}']
        lines += [f'  inside box 0 0 0 {box} {box} {box}']
        lines += ['end structure']
        lines += ['']
    with open('pack.inp', 'w') as f:
        for line in lines:
            f.write(f'{line}\n')


def main():
    parser = argparse.ArgumentParser(
        description=('Replicating the setup procedure in fftool. '
                     'Pass in xyz files along with a box size in Å'))

    parser.add_argument('-b', '--box', help='box length in Å')
    parser.add_argument(
        'files',
        nargs='+',
        help=('n1 xyz1 [n2 xyz2 ...], where n_i is the number of molecules '
              'and xyz_i is a xyz structure for one molecule'))
    args = parser.parse_args()
    if any(x is None
           for x in (args.files, args.box)) or len(args.files) % 2 != 0:
        parser.print_help()
        sys.exit()

    nmols = args.files[::2]
    xyzs = args.files[1::2]

    writepackmol(nmols, xyzs, args.box)
    os.system('packmol < pack.inp')


if __name__ == "__main__":
    main()
