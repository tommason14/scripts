#!/usr/bin/env python3

import os
import re
import sys
from autochem import Molecule
from glob import glob


def read_file(fname):
    with open(fname) as f:
        return f.readlines()


def replace_name(inp, name):
    newinp = []
    for line in inp:
        if "name" in line:
            newinp.append(line.replace("name", name))
        else:
            newinp.append(line)
    return newinp


def add_mol_info(inp, mol):
    """
    Add molecule info to the psi4 template
    """

    atoms = [frag["atoms"] for frag in mol.fragments.values()]
    charges = [frag["charge"] for frag in mol.fragments.values()]
    mult = [frag["multiplicity"] for frag in mol.fragments.values()]

    if any(len(x) == 1 for x in (atoms, charges, mult)):
        raise AttributeError(
            (
                "Two molecules are needed, but autochem can only find one. Your molecules may be too close together."
                "If not, a molecule definition might need to be added to ~/.config/autochem/molecules.txt"
            )
        )

    count = 0
    start = 0
    end = 0
    # split on '0 1' or '-1 1'
    chunks = []
    for ind, line in enumerate(inp, 1):
        if re.search("-?[0-9]+ [0-9]+", line):
            count += 1
            end = ind
            chunks.append(inp[start : end - 1])
            start = ind
        if count == 2:
            chunks.append(inp[start:])

    newinp = (
        chunks[0]
        + [f"{charges[0]} {mult[0]}\n"]
        + [a.as_xyz() for a in atoms[0]]
        + chunks[1]
        + [f"{charges[1]} {mult[1]}\n"]
        + [a.as_xyz() for a in atoms[1]]
        + ["\n"]
        + chunks[2]
    )
    return newinp


def main():
    if "template.inp" not in os.listdir(".") or "template.job" not in os.listdir("."):
        print("Error: make sure template.inp and template.job are in this folder")
        sys.exit(1)
    for xyz in glob("*xyz"):
        name = xyz.replace(".xyz", "")
        newdir = os.path.join(os.getcwd(), name + "-sapt")
        if not os.path.isdir(newdir):
            os.mkdir(newdir)
        mol = Molecule(xyz)
        inp = read_file("template.inp")
        inp = replace_name(inp, name)
        inp = add_mol_info(inp, mol)
        job = read_file("template.job")
        job = replace_name(job, name)
        inpname = os.path.join(newdir, name + ".inp")
        jobname = os.path.join(newdir, name + ".job")
        with open(inpname, "w") as f:
            for line in inp:
                f.write(line)
        with open(jobname, "w") as f:
            for line in job:
                f.write(line)


if __name__ == "__main__":
    main()
