#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import get_font
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
    "-o", "--output", help="Prefix of saved files, default=rmsf", default="rmsf"
)
parser.add_argument(
    "-p", "--plot", help="Create plot of RMSF per residue", action="store_true"
)

args = parser.parse_args()
u = mda.Universe(args.coords, *args.trajectory)

# create average structure of the first frame
average = align.AverageStructure(u, u, select=args.selection, ref_frame=0).run()
ref = average.results.universe
aligner = align.AlignTraj(u, ref, select=args.selection, in_memory=True).run()
selection = u.select_atoms(args.selection)
rmsf = rms.RMSF(selection, verbose=True).run()
# RMSF is computed per atom, but the convention is to average over each residue
# selection.resnums = 1 per atom, so averaged = RMSF for each residue
df = (
    pd.DataFrame({"resid": selection.resnums, "rmsf": rmsf.results.rmsf})
    .groupby("resid")["rmsf"]
    .mean()
    .reset_index()
)
df.to_csv(args.output + ".csv", index=False)
if args.plot:
    sns.set(
        style="ticks",
        font=get_font(),
        font_scale=1.2,
        rc={"mathtext.default": "regular"},
    )
    p = sns.lineplot(x="resid", y="rmsf", ci=None, data=df)
    p.set_xlabel("Residue number")
    p.set_ylabel(r"RMSF ($\AA$)")
    plt.savefig(args.output + ".pdf", dpi=300)
