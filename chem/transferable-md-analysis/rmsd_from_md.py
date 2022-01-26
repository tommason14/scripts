#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import MDAnalysis.transformations as trans
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
    "-w", "--mass-weighted", help="Mass-weight RMSD per atom", action="store_true"
)
parser.add_argument(
    "--centre",
    help="Pre-process the trajectory so that the selected atoms are centred in the each frame before computing RMSDs",
    action="store_true",
)
parser.add_argument(
    "-o", "--output", help="Prefix of saved files, default=rmsd", default="rmsd"
)
parser.add_argument(
    "-p", "--plot", help="Create plot of RMSD against time", action="store_true"
)

args = parser.parse_args()
weights = "mass" if args.mass_weighted else None
u = mda.Universe(args.coords, *args.trajectory)
if args.centre:
    selection = u.select_atoms(args.selection)
    notsel = mda.AtomGroup([a for a in u.atoms if a not in selection])
    transforms = [
        trans.unwrap(selection),
        trans.center_in_box(selection, wrap=True),
        trans.wrap(notsel),
    ]
    u.trajectory.add_transformations(*transforms)
rmsd = rms.RMSD(u, u, args.selection, weights=weights).run()

df = (
    pd.DataFrame(rmsd.results.rmsd, columns=["Frame", "Time (ps)", "RMSD"])
    .assign(**{"Time (ns)": lambda x: x["Time (ps)"] / 1000})  # expects keyword args
    .drop("Time (ps)", axis=1)
)
df.to_csv(args.output + ".csv", index=False)
if args.plot:
    sns.set(
        style="ticks",
        font='DejaVu Sans',
        font_scale=1.2,
        rc={"mathtext.default": "regular"},
    )
    p = sns.lineplot(x="Time (ns)", y="RMSD", ci=None, data=df)
    p.set_ylabel(r"RMSD ($\AA$)")
    plt.savefig(args.output + ".pdf", dpi=300)
