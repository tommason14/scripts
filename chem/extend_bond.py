#!/usr/bin/env python3
import sys
import math
import argparse
from autochem import PeriodicTable as PT, Atom, read_xyz, write_xyz

parser = argparse.ArgumentParser()
parser.add_argument(
    "-s",
    "--still",
    help="Atom(s) remaining stationary. Give the positions in the \
list of coordinates as a comma-separated list. e.g. 1,2,3,4",
    type=lambda s: [int(item) for item in s.split(",")],
)
parser.add_argument(
    "-m",
    "--move",
    help="Atom(s) to move. Give the positions in the \
list of coordinates as a comma-separated list. e.g. 1,2,3,4",
    type=lambda s: [int(item) for item in s.split(",")],
)
parser.add_argument(
    "-b",
    "--by",
    help="Distance(s) to extend the bond by, in Å. Give the distances \
as a comma-separated list. e.g. 0.1,0.2,0.3",
    type=lambda s: [float(item) for item in s.split(",")],
)
parser.add_argument("-u", "--using", help="XYZ file to modify", action="store")

args = parser.parse_args()


if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

if len(args.still) != len(args.move):
    print("Error: Number of stationary and moving atoms must be the same.")
    sys.exit()


def extend_bond(atom1, atom2, atoms, dist):
    """
    Finding vector between atoms 1 and 2, then scale that vector by a ratio of desired length over
    the original, and moving atom2 by adding the scaled vector back onto atom1's coordinates and
    reassigning them to atom2. Returns a new list of atoms
    """
    a_one = atoms[atom1 - 1]
    a_two = atoms[atom2 - 1]
    vec = a_one.vector_to(a_two)
    original_distance = (sum(i ** 2 for i in vec)) ** 0.5
    scaling_factor = (original_distance + dist) / original_distance
    scaled = tuple((coord * scaling_factor) for coord in vec)
    added = tuple((i + j) for i, j in zip(a_one.coords, scaled))
    atoms[atom2 - 1] = Atom(a_two.symbol, coords=added)
    print(f"{a_one.symbol}-{a_two.symbol} ({atom1}-{atom2}) extended by {dist} Å")
    return atoms


def group_items(list1, list2):
    """Zips two lists together into a list of tuples"""
    return [(i, j) for i, j in zip(list1, list2)]


def make_changes(list_of_atoms, atoms_changing, dist):
    """Return altered coordinates"""
    for group in atoms_changing:  # [(1,2), (3,4)]
        static, mobile = group
        coords = extend_bond(static, mobile, list_of_atoms, dist)
    return coords


def main():
    file, stationary, moving, extensions = args.using, args.still, args.move, args.by
    atoms = group_items(stationary, moving)
    for dist in extensions:
        print("-" * 30)
        print(f"Moving by {dist} Å")
        coordinates = make_changes(read_xyz(file), atoms, dist)
        write_xyz(coordinates, f"modded_{file[:-4]}_{dist}.xyz")
        print(f"Written to modded_{file[:-4]}_{dist}.xyz")


main()
