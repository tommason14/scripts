#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import get_font, completion
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--coords", help="Coordinate file (.gro/.pdb)", required=True)
parser.add_argument(
    "-t",
    "--trajectory",
    help="Trajectory file(s) (.xtc/.dcd/.pdb)",
    nargs="+",
    required=True,
)
parser.add_argument("-s", "--selection", help="MDAnalysis selection")
parser.add_argument(
    "-o", "--output", help="Prefix of saved files, default=rg", default="rg"
)
parser.add_argument(
    "-p", "--plot", help="Create plot of Rg against time", action="store_true"
)
parser.add_argument(
    "-r",
    "--per-residue",
    help="Compute radius of gyration for each residue in the selection",
    action="store_true",
)
args = parser.parse_args()
u = mda.Universe(args.coords, *args.trajectory)
selection = u.select_atoms(args.selection)

if args.per_residue:
    resids = np.unique(selection.resids)
    rg = []
    for ts in completion(u.trajectory):
        for res in resids:
            residue = selection.select_atoms(f"resid {res}")
            rg.append((ts.time, res, residue.radius_of_gyration()))
    df = pd.DataFrame(rg, columns=["Time (ps)", "Resid", "Rg"])
else:
    rg = []
    for ts in completion(u.trajectory):
        rg.append((ts.time, selection.radius_of_gyration()))
    df = pd.DataFrame(rg, columns=["Time (ps)", "Rg"])

df["Time (ns)"] = df["Time (ps)"] / 1000
df.drop("Time (ps)", axis=1, inplace=True)
df.to_csv(args.output + ".csv", index=False)

if args.plot:
    sns.set(
        style="ticks",
        font=get_font(),
        font_scale=1.2,
        rc={"mathtext.default": "regular"},
    )
    if args.per_residue:
        p = sns.lineplot(x="Time (ns)", y="Rg", hue="Resid", ci=None, data=df)
    else:
        p = sns.lineplot(x="Time (ns)", y="Rg", ci=None, data=df)
    p.set_ylabel(r"R$_{gyr}$ ($\AA$)")
    plt.savefig(args.output + ".pdf", dpi=300)
