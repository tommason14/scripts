#!/usr/bin/env python3

import MDAnalysis as mda
import argparse
import numpy as np
import pandas as pd
import sys
from utils import completion


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute dipole moments via MDAnalysis"
    )
    parser.add_argument(
        "-p",
        "--topology",
        help="Topology file that must contain charges (.top/.psf)",
        required=True,
    )
    parser.add_argument(
        "-t", "--traj", help="Trajectory file (.dcd/.xtc)", required=True
    )
    parser.add_argument(
        "-o",
        "--csv",
        help="csv to save all vectors to, default=dipole_vectors.csv",
        default="dipole_vectors.csv",
    )
    parser.add_argument("--sel", help="MDAnalysis atomgroup.select_atoms() expression")
    parser.add_argument(
        "-nodrudes", help="Remove drude particles from trajectory", action="store_true"
    )
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
        "-m",
        "--maxval",
        help=(
            "Maximum dipole moment used when creating a histogram "
            "of moments for each molecule"
        ),
        default=10,
        type=float,
    )
    parser.add_argument(
        "-n",
        "--nbins",
        help="Number of intervals to create for the dipole moment histogram",
        default=50,
        type=int,
    )
    parser.add_argument(
        "-b",
        "--histname",
        help=(
            "Filename of the csv containing the histogram of total dipole moments, "
            "default = dipole_moments.csv"
        ),
        default="dipole_moments.csv",
    )
    return parser.parse_args()


def _calc_dipoles(selection):
    """
    Returns a dictionary with list items containing the dipole moments (in the x, y, z dimensions) of each molecule 
    in the trajectory. The dipole moment is calculated relative to the centre of mass of each
    molecule, given in e.Å.
    TODO: Compute dipole moments via numpy instead of iterating over each atom
    """
    dipoles = {res: [] for res in np.unique(selection.atoms.resnames)}
    for res in selection.residues:
        dipole_moment = np.array([0, 0, 0])
        residue_com = res.atoms.center_of_mass()
        for atom in res.atoms:
            dipole_moment = dipole_moment + atom.charge * (atom.position - residue_com)
        dipoles[res.resname].append(dipole_moment)
    return dipoles


def _calc_dipoles_all_frames(selection, frames):
    """
    Loops over frames and computes dipoles for every molecule, returning them as 
    lists of 1x3 arrays for each molecule, stored inside a dictionary of {res: [[x1 y1 z1], [x2, y2 z2]]}
    """
    dipole_moments = {res: [] for res in np.unique(selection.atoms.resnames)}
    numframes = len(frames)
    for ts in completion(frames):
        dipoles = _calc_dipoles(selection)
        for res, dipole_moment_list in dipoles.items():
            dipole_moments[res] += dipole_moment_list
    return dipole_moments


def compute_moments(args):
    """
    Computes dipole moments for each molecule (residue) in every frame of the trajectory, 
    and then returns the values in a pandas DataFrame with headers:
    resname, Dx, Dy, Dz, Dtot.
    For example, for a system of 500 molecules, analysing 10 frames will produce a DataFrame
    of 5000 rows.
    All values are given in Debyes.
    """
    u = mda.Universe(args.topology, args.traj)

    selection = u.select_atoms("all")
    if args.nodrudes:
        selection = selection.select_atoms("not name D*")
    if args.sel:
        selection = selection.select_atoms(args.sel)

    if args.start and args.end:
        frames = u.trajectory[args.start : args.end]
    elif args.start:
        frames = u.trajectory[args.start :]
    elif args.end:
        frames = u.trajectory[: args.end]
    else:
        frames = u.trajectory

    dipole_moments = _calc_dipoles_all_frames(selection, frames)
    df = {"resname": [], "Dx": [], "Dy": [], "Dz": []}
    for res, moments in dipole_moments.items():
        for moment in moments:
            df["resname"].append(res)
            df["Dx"].append(moment[0])
            df["Dy"].append(moment[1])
            df["Dz"].append(moment[2])

    df = pd.DataFrame(df)
    eA_to_debye = 0.20819434  # 0.208 e.angstrom -> 1 Debye
    df["Dx"] = df["Dx"] / eA_to_debye
    df["Dy"] = df["Dy"] / eA_to_debye
    df["Dz"] = df["Dz"] / eA_to_debye
    df["Dtot"] = (df["Dx"] ** 2 + df["Dy"] ** 2 + df["Dz"] ** 2) ** 0.5
    return df


def group_moments_into_intervals(df, maxval=10, n=50):
    """
    Group dipole moments of each molecule into n intervals for easier plotting, with a much smaller file
    produced.
    """
    minval = 0
    bins = np.linspace(minval, maxval, n)
    counts = (
        df.groupby(["resname", pd.cut(df["Dtot"], bins)])
        .size()
        .unstack()
        .T.reset_index()
    )
    counts["Dtot"] = counts["Dtot"].apply(lambda x: round(x.mid, 2))
    # move to first column
    col = counts.pop("Dtot")
    counts.insert(0, "Dtot", col)
    return counts


def main():
    args = parse_args()
    dipoles = compute_moments(args)
    dipoles.to_csv(args.csv, index=False)
    binned = group_moments_into_intervals(dipoles, maxval=args.maxval, n=args.nbins)
    binned.to_csv(args.histname, index=False)


if __name__ == "__main__":
    main()
