#!/usr/bin/env python3

import MDAnalysis as mda
import argparse


def completion(iterable):
    """
    Write percentage completion of a loop to the screen
    """
    import sys

    total = len(iterable)
    for num, val in enumerate(iterable):
        yield val
        sys.stdout.write(f"\r{num/total*100:.2f} % complete")
    sys.stdout.write("\n")


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
    return parser.parse_args()


def convert_traj(args):
    u = mda.Universe(args.coords, *args.traj)
    if args.output is None:
        args.output = args.traj.rsplit(".")[0] + ".pdb"

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

    with mda.Writer(args.output, sel.n_atoms) as f:
        for _ in completion(frames):
            f.write(sel)


def main():
    args = parse_args()
    convert_traj(args)


if __name__ == "__main__":
    main()
