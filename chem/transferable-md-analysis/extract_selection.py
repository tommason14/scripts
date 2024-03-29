#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis import transformations
import argparse
from utils import completion


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sel", help="MDAnalysis selection", required=True)
    parser.add_argument(
        "files",
        nargs="+",
        help=(
            "Coordinate/topology file followed trajectory file(s). "
            "If no trajectory files are passed in, a selection of the topology is returned in gro "
            "format. If trajectory files are passed in, an xtc file is saved with the selection over "
            "the trajectory."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output name without extension. Default = selection, "
            "meaning that selection.gro and selection.xtc will be saved."
        ),
        default="selection",
    )
    parser.add_argument(
        "-wrap",
        help="Use this parameter to output a wrapped trajectory.",
        action="store_true",
    )
    return parser.parse_args()


def main():
    args = arguments()
    gro = args.output + ".gro"
    xtc = args.output + ".xtc"
    u = mda.Universe(*args.files)
    selection = u.select_atoms(args.sel)
    if args.wrap:
        ag = selection.atoms
        transform = mda.transformations.wrap(ag, compound="residues")
        u.trajectory.add_transformations(transform)
    if len(args.files) == 1:
        selection.write(gro)
    else:
        selection.write(gro)
        with mda.Writer(xtc, selection.n_atoms) as f:
            for frame in completion(u.trajectory):
                f.write(selection)


if __name__ == "__main__":
    main()
