#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import warnings


warnings.filterwarnings(
    "ignore", message="warning: two consecutive residues with same number.*"
)
# warnings of residues with same number - but that's not the case in the gro file, and also not
# the
# case when looking at the resid ... so just ignore the warning
# this is the code:
# if residue is none or residue.resseq != thisresnum or residue.name != thisresname:
#     if residue is not none and residue.name != thisresname:
#         warnings.warn("warning: two consecutive residues with same number (%s, %s)" %
#         (thisresname, residue. name))
# issue with resseq i guess?

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--coords", help="Coordinate file (.gro/.pdb)", required=True)
parser.add_argument(
    "-t",
    "--trajectory",
    help="Trajectory file(s) (.xtc/.dcd/.pdb)",
    nargs="+",
    required=True,
)
parser.add_argument("-s", "--selection", help="Optional MDTraj selection")
parser.add_argument(
    "-o", "--output", help="Prefix of saved files, default=sasa", default="sasa"
)

parser.add_argument(
    "-p", "--plot", help="Create plot of SASA against time", action="store_true"
)

parser.add_argument(
    "-r",
    "--per-residue",
    help="Compute SASA for each residue in the selection",
    action="store_true",
)
args = parser.parse_args()
traj = md.load(args.trajectory, top=args.coords)
selection = traj.top.select(args.selection)
traj = traj.atom_slice(selection)

print("Starting SASA calculation")
if args.per_residue:
    sasa = md.shrake_rupley(traj, mode="residue")
    df = pd.DataFrame(sasa, columns=list(traj.top.residues))
else:
    sasa = md.shrake_rupley(traj)
    # look at total area
    sasa = sasa.sum(axis=1)
    df = pd.DataFrame(sasa, columns=["SASA"])
print("Finished SASA calculation")

df.insert(0, "Time (ns)", traj.time / 1000)
df.to_csv(args.output + ".csv", index=False)

if args.plot:
    sns.set(
        style="ticks",
        font="DejaVu Sans",
        font_scale=1.2,
        rc={"mathtext.default": "regular"},
    )
    if args.per_residue:
        tidy = df.melt(id_vars="Time (ns)", var_name="Residue", value_name="SASA")
        # normally only want a specific residue or two over time, but could have a lot...
        # so just plot averages
        tidy["Resid"] = tidy["Residue"].astype(str).str[3:].astype(int)
        p = sns.pointplot(
            x="Resid", y="SASA", ci="sd", data=tidy, join=False, errwidth=2, capsize=0.4
        )
        if len(tidy["Resid"].unique()) > 10:
            # presumably looking at all residues of a protein; show every 10
            p.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
            p.xaxis.set_major_locator(ticker.MultipleLocator(base=10))
            # also make plot wide so that all error bars are visible
            p.figure.set_figwidth(10)
    else:
        p = sns.lineplot(x="Time (ns)", y="SASA", ci=None, data=df)
    p.set_ylabel(r"SASA (nm$^2$)")
    plt.savefig(args.output + ".pdf", dpi=300)
