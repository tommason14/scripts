#!/usr/bin/env python3
from autochem import Settings, GaussJob, Molecule
import os
from glob import glob

cwd = os.getcwd()

sett = Settings()
sett.input.method = "M062X"
sett.input.basis = "Gen"
sett.input.polar = "Dipole"

# monarch settings
sett.meta.ncpus = 16
sett.meta.mem = "64gb"
sett.meta.nodemem = "64gb"
sett.meta.time = "2:00:00"
sett.meta.partition = "comp,short"


def get_basis(
    atoms,
    fname=os.path.expanduser("~/.local/scripts/chem/polarisabilities/sadlej.gbs"),
):
    """
    Returns the required basis functions as a list
    """
    basis_set = []
    tmp = []
    found = False
    with open(fname) as f:
        for line in f:
            if "****" in line:
                if len(tmp) > 0:
                    basis_set.append(tmp)
                    tmp = []
                continue
            tmp.append(line)
    # basis_set is a list of basis functions for each basis set,
    # turn it into a dict of {atom: basis}
    basis_set = {bset[0].split()[0]: bset for bset in basis_set}
    basis = []
    for atom, bset in basis_set.items():
        if atom in atoms:
            basis += bset + ["****\n"]
    return basis


get_basis(["H", "C", "O"])
fields = {
    "plusX": "X+8",
    "minusX": "X-8",
    "plusY": "Y+8",
    "minusY": "Y-8",
    "plusZ": "Z+8",
    "minusZ": "Z-8",
}

for xyz in glob("*xyz"):
    xyzdir = xyz.replace(".xyz", "")
    mol = Molecule(xyz)
    atoms = list(set([a.symbol for a in mol.coords]))
    try:
        basis = get_basis(sorted(atoms)) + ["\n"]
    except ValueError:
        print("Some or all atoms are not included in the SadleJ basis set")
    if os.path.isdir(xyzdir):
        print(f"Skipping {xyzdir}")
        continue
    os.mkdir(xyzdir)
    for name, field in fields.items():
        os.mkdir(f"{xyzdir}/{name}")
        os.chdir(f"{xyzdir}/{name}")
        _sett = sett.copy()
        _sett.input.field = field
        GaussJob(f"../../{xyz}", settings=_sett)  # makes spec.job
        # add basis set
        with open("spec.job") as f:
            job = f.readlines()
        job = job[:-1] + basis + job[-1:]
        with open("spec.job", "w") as new:
            for line in job:
                new.write(line)
        os.chdir(cwd)
