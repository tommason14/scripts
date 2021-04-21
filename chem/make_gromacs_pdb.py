#!/usr/bin/env python3
import argparse
import os
import sys


def writepackmol(mols, files, box):
    lines = [
        "tolerance 2.0",
        "filetype pdb",
        "output pack.pdb",
        "add_amber_ter",
        "add_box_sides 1.0",
        "",
    ]

    for n, pdb in zip(mols, files):
        lines += [f"structure {pdb}"]
        lines += [f"  number {n}"]
        lines += [f"  inside box 0 0 0 {box} {box} {box}"]
        lines += ["end structure"]
        lines += [""]
    with open("pack.inp", "w") as f:
        for line in lines:
            f.write(f"{line}\n")


def add_box_dimensions(fname, box):
    with open(fname) as f:
        inp = f.readlines()
    index = 0
    for ind, line in enumerate(inp):
        if "ATOM" in line:
            index = ind
            break
    box = float(box)
    inp.insert(
        index, f"CRYST1   {box:.3f}   {box:.3f}   {box:.3f}  90.00  90.00  90.00 P 1\n"
    )
    with open(fname, "w") as f:
        for line in inp:
            f.write(line)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Replicating the setup procedure in fftool. "
            "Pass in pdb files with residues labelled, "
            "along with a box size in Å"
        )
    )

    parser.add_argument("-b", "--box", help="box length in Å")
    parser.add_argument(
        "files",
        nargs="+",
        help=(
            "n1 pdb1 [n2 pdb2 ...], where n_i is the number of molecules "
            "and pdb_i is a pdb structure for one molecule"
        ),
    )
    args = parser.parse_args()
    if any(x is None for x in (args.files, args.box)) or len(args.files) % 2 != 0:
        parser.print_help()
        sys.exit()

    nmols = args.files[::2]
    pdbs = args.files[1::2]

    writepackmol(nmols, pdbs, args.box)
    os.system("packmol < pack.inp")
    # add_box_dimensions("pack.pdb", args.box)


if __name__ == "__main__":
    main()
