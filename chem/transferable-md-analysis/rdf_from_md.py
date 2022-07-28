#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import use
import matplotlib.pyplot as plt

use("Agg")


class COM:
    def __init__(self, ag):
        self._groups = list(ag.groupby("resids").values())

    def __len__(self):
        return len(self._groups)

    @property
    def positions(self):
        return np.array([g.center_of_mass() for g in self._groups], dtype=np.float32)

    @property
    def universe(self):
        return self._groups[0].universe


parser = argparse.ArgumentParser()
parser.add_argument(
    "-c",
    "--coords",
    help="Coordinate/topology file (.gro/.pdb/.psf/.top etc)",
    required=True,
)
parser.add_argument(
    "-t",
    "--trajectory",
    help="Trajectory file(s) (.dcd/.pdb/.xtc etc)",
    nargs="+",
    required=True,
)
parser.add_argument("-ref", help="Reference MDAnalysis selection", required=True)
parser.add_argument("-sel", help="MDAnalysis selection to search for", required=True)
parser.add_argument(
    "-o", "--output", help="Output prefix. Default = rdf", default="rdf",
)
parser.add_argument(
    "-p", "--plot", help="Plot data, saved with the output prefix", action="store_true",
)
parser.add_argument(
    "-s",
    "--start",
    help="Smallest distance to consider in Angstroms, default=0",
    default=0,
    type=float,
)
parser.add_argument(
    "-e",
    "--end",
    help="Largest distance to consider in Angstroms, default=15",
    default=15,
    type=float,
)
parser.add_argument(
    "-com",
    "--centre-of-mass",
    help="Compute RDFs using the centre of mass of each selection",
    action="store_true",
)
parser.add_argument(
    "-i",
    "--intermolecular",
    help="Disregard interactions of atoms in the same residue",
    action="store_true",
)
parser.add_argument(
    "-b",
    "--bins",
    help="Number of intervals to bin distances into, default=200",
    default=200,
    type=int,
)

args = parser.parse_args()

u = mda.Universe(args.coords, *args.trajectory)
ref = u.select_atoms(args.ref)
sel = u.select_atoms(args.sel)

print(f"Reference selection contains {ref.n_atoms} atoms")
print(f"Mobile selection contains {sel.n_atoms} atoms")

if args.centre_of_mass:
    ref = COM(ref)
    sel = COM(sel)

if args.intermolecular:
    ref_atoms = len(ref.residues[0].atoms)
    sel_atoms = len(sel.residues[0].atoms)
    rdf = InterRDF(
        ref,
        sel,
        nbins=args.bins,
        range=(args.start, args.end),
        exclusion_block=(ref_atoms, sel_atoms),
    )
else:
    rdf = InterRDF(ref, sel, nbins=args.bins, range=(args.start, args.end))
rdf.run(verbose=True)
df = pd.DataFrame({"bins": rdf.bins, "rdf": rdf.rdf})
df.to_csv(args.output + ".csv", index=False)

if args.plot:
    sns.set(
        style="white",
        font="DejaVu Sans",
        font_scale=1.2,
        rc={"mathtext.default": "regular"},
    )
    p = sns.lineplot(x="bins", y="rdf", ci=None, data=df)
    p.set_xlabel(r"Distance (${\AA}$)")
    p.set_ylabel("g(r)")
    plt.tight_layout()
    plt.savefig(args.output + ".png", dpi=300)
