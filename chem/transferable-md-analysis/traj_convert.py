#!/usr/bin/env python3

import MDAnalysis as mda
from MDAnalysis import transformations
import argparse
from utils import completion


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert trajectories to different formats via MDAnalysis"
    )
    parser.add_argument(
        "-c", "--coords", help="Coordinate file (.gro/.psf/.pdb)", required=True
    )
    parser.add_argument(
        "-t", "--traj", help="Trajectory file(s) (.dcd/.xtc)", nargs="+", required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename to write to, default is a pdb with the same name as the trajectory",
    )
    parser.add_argument(
        "-nodrudes", help="Remove drude particles from trajectory", action="store_true"
    )
    parser.add_argument(
        "-novirtuals", help="Remove virtual sites from trajectory", action="store_true"
    )
    parser.add_argument("--sel", help="Selection of atoms, for example 'resname bf4'")
    parser.add_argument(
        "-s",
        "--start",
        help="Frame to start at, default is to include all frames",
        type=int,
    )
    parser.add_argument(
        "-e",
        "--end",
        help="Frame to end at, default is to include all frames",
        type=int,
    )
    parser.add_argument(
        "--wrap",
        help="Wrap trajectory so that each residue remains intact",
        action="store_true",
    )
    return parser.parse_args()


def convert_traj(args):
    u = mda.Universe(args.coords, *args.traj)
    if args.output is None:
        args.output = args.traj[-1].rsplit(".")[0] + ".pdb"

    sel = mda.AtomGroup(u.atoms)  # so we can use sel.n_atoms later
    if args.nodrudes:
        sel = sel.select_atoms("not name D*")
    if args.novirtuals:
        sel = sel.select_atoms("not name V*")
    if args.sel:
        sel = sel.select_atoms(args.sel)

    if args.start and args.end:
        frames = u.trajectory[args.start : args.end]
    elif args.start:
        frames = u.trajectory[args.start :]
    elif args.end:
        frames = u.trajectory[: args.end]
    else:
        frames = u.trajectory

    if args.wrap:
        print("Wrapping trajectory")
        ag = sel.atoms
        transform = mda.transformations.wrap(ag, compound="residues")
        # can't add transformations to a slice...
        u.trajectory.add_transformations(transform)

    with mda.Writer(args.output, sel.n_atoms) as f:
        for _ in completion(frames):
            f.write(sel)


def main():
    args = parse_args()
    convert_traj(args)


if __name__ == "__main__":
    main()
