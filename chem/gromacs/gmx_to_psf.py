#!/usr/bin/env python3
import argparse
import parmed as pmd
import sys

parser = argparse.ArgumentParser(
    description="Convert gromacs coordinates and topology into a charmm PSF format"
)
parser.add_argument("-g", "--gro", help="gro file", required=True)
parser.add_argument("-t", "--top", help="topology file", required=True)
parser.add_argument("-o", "--out", help="psf format", required=True)
args = parser.parse_args()
pmd.load_file(args.top, xyz=args.gro).save(args.out)
