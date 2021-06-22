#!/usr/bin/env python3
from mstools.forcefield import ForceField, PaduaLJScaler
from mstools.topology import Topology
from mstools.simsys import System
import argparse


def make_system(molecules, numbers, ff, box=None, virtual_sites=False, drudes=False):
    from shutil import which
    from mstools.wrapper import Packmol

    packmol_path = which("packmol")
    packmol = Packmol(packmol_path)

    # convert to nm - nargs gives a list
    if len(box) == 1:
        box = [b / 10 for b in box] * 3
    else:
        box = [b / 10 for b in box]
    top = Topology(molecules)
    top.generate_angle_dihedral_improper()
    if virtual_sites:
        top.generate_virtual_sites(ff)
    if drudes:
        top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)
    top.cell.set_box(box)
    top.scale_with_packmol(numbers, packmol)
    return System(top, ff)


def read_args():
    parser = argparse.ArgumentParser(
        description="Create coordinate and topology files for OpenMM using OPLS-based forcefields"
    )
    parser.add_argument(
        "-f",
        "--files",
        help="xyz/pdb/zmat structure files to combine",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--numbers",
        help="Number of each structure file",
        nargs="+",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-ff", "--forcefield", help="Forcefield files to use", nargs="+", required=True
    )
    parser.add_argument(
        "-b",
        "--box",
        help=(
            "Box size in angstroms. Can be either one or three values; "
            "if one is given, this will be used as the length of each dimension"
        ),
        nargs="+",
        type=float,
        required=True,
    )

    parser.add_argument(
        "-v", "--virtuals", help="Include virtual sites", action="store_true"
    )
    parser.add_argument(
        "-d", "--drudes", help="Include drude particles", action="store_true"
    )

    return parser.parse_args()


def read_mols(files):
    mols = []
    for f in files:
        mol = Topology.open(f).molecules[0]
        if f.endswith("xyz"):
            mol.guess_connectivity_from_ff(ff)
        mols.append(mol)


def main():
    args = read_args()
    mols = read_mols(args.files)
    ff = ForceField.open(*args.forcefield)
    system = make_system(
        mols,
        args.numbers,
        ff,
        box=args.box,
        virtual_sites=args.virtuals,
        drudes=args.drudes,
    )

    system.export_gromacs(gro_out="conf.gro", top_out=None, mdp_out=None)
    system.export_charmm(psf_out="topol.psf", prm_out="ff.prm", pdb_out=None)


if __name__ == "__main__":
    main()
