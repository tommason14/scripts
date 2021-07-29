#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c", "--coords", help="Coordinate/topology file (.gro/.pdb/.psf/.top etc)"
)
parser.add_argument(
    "-t", "--trajectory", help="Trajectory file(s) (.dcd/.pdb/.xtc etc)", nargs="+"
)
parser.add_argument("-ref", help="Reference MDAnalysis selection")
parser.add_argument("-sel", help="MDAnalysis selection to search for")
parser.add_argument(
    "-o",
    "--csv",
    help="Filename of csv containing RDF data to write to. Default = rdf.csv",
    default="rdf.csv",
)
parser.add_argument(
    "-p",
    "--plot",
    help=(
        "Create a plot that is saved as rdf.pdf. If a csv "
        "filename is given, say example.csv, the plot will be named as example.pdf"
    ),
    action="store_true",
)

args = parser.parse_args()

u = mda.Universe(args.coords, *args.trajectory)
ref = u.select_atoms(args.ref)
sel = u.select_atoms(args.sel)

rdf = InterRDF(ref, sel, nbins=200)
rdf.run(verbose=True)
df = pd.DataFrame({"bins": rdf.bins, "rdf": rdf.rdf})
df.to_csv(args.csv, index=False)

if args.plot:
    sns.set(style="white", font="Nimbus Sans", rc={"mathtext.default": "regular"})
    p = sns.lineplot(x="bins", y="rdf", ci=None, data=df)
    p.set_xlabel(r"Distance (${\AA}$)")
    p.set_ylabel("g(r)")
    plt.tight_layout()
    plt.savefig(args.csv.replace("csv", "pdf"), dpi=300)
