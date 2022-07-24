#!/usr/bin/env python3
import argparse
from matplotlib import use
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns

use("Agg")

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--coords", help="Coordinate file (.gro/.pdb)", required=True)
parser.add_argument(
    "-t", "--trajectory", help="Trajectory file (.xtc/.dcd/.pdb)", required=True,
)
parser.add_argument(
    "-s",
    "--selection",
    help="MDTraj selection of a subset of atoms, default=backbone",
    default="backbone",
)
parser.add_argument(
    "-o", "--output", help="Prefix of saved files, default=dssp", default="dssp"
)
parser.add_argument(
    "-p",
    "--plot",
    help="Create plot of types of secondary structure against time",
    action="store_true",
)
parser.add_argument(
    "-r",
    "--rolling",
    help="Number of frames to apply a rolling average over, default=1, meaning no averaging",
    default=1,
    type=int,
)


def analyse(coords, trajectory, selection):
    """
    Run mdtraj.compute_dssp on the trajectory
    """
    traj = md.load(trajectory, top=coords)
    traj_selection = traj.atom_slice(traj.topology.select(selection))
    dssp = md.compute_dssp(traj_selection)
    # codes for each residue in each frame
    # summarise code per frame
    codes = np.unique(dssp.ravel())
    data = {i: [] for i in codes}
    for frame in dssp:
        codes_, counts = np.unique(frame, return_counts=True)
        # if a code doesn't appear in the frame, add a 0
        counted = []
        for code, count in zip(codes_, counts):
            data[code].append(count)
            counted.append(code)
        if len(counted) < len(codes):
            for code in codes:
                if code not in counted:
                    data[code].append(0)

    # remove the residues that aren't proteins
    return (
        pd.DataFrame(data)
        .drop("NA", axis=1)
        .rename(columns={"H": "Helix", "E": "Strand", "C": "Coil"})
        .assign(time_ns=traj.time / 1000)  # MDTraj default is ps
    )


if __name__ == "__main__":
    args = parser.parse_args()
    df = analyse(args.coords, args.trajectory, args.selection)
    df.to_csv(args.output + ".csv", index=False)

    if args.plot:
        sns.set(style="ticks", font="DejaVu Sans", font_scale=1.1)
        p = (
            df.set_index("time_ns")
            .rolling(window=args.rolling)
            .mean()
            .apply(
                lambda x: 100 * x / x.sum(), axis=1
            )  # percentages computed for each frame
            .plot(color=["#545775", "#FFA686", "#53B3CB"])
        )
        p.set_xlabel("Time (ns)")
        p.set_ylabel("Proportion (%)")
        plt.tight_layout()
        plt.savefig(args.output + ".png", dpi=300)
