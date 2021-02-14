#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csv", help="CSV written by Travis", required=True)
parser.add_argument("-o", "--output", help="Filename to save image to")
parser.add_argument(
    "-d", "--distance", help="Maximum distance to plot in Angstroms", type=int
)
parser.add_argument(
    "-l",
    "--levels",
    help="Number of contours to plot. Default = 30",
    default=30,
    type=int,
)
args = parser.parse_args()

df = pd.read_csv(args.csv, sep=";\s+", engine="python")
x = df["Distance / pm"]
y = df["# Tau / ps"]
z = df["Occurrence"]

# convert distance to angstroms
x = x / 100
# convert time to ns
y = y / 1000

xmin = x.min()
if args.distance:
    xmax = args.distance
else:
    xmax = x.max()

colour = "bwr"
sns.set(style="whitegrid", font="Helvetica")
# display mathematical symbols with a regular font (the default is to write them in italics)
plt.rcParams["mathtext.default"] = "regular"
ax = plt.tricontourf(x, y, z, levels=args.levels, cmap=colour, linestyles="solid")
bar = plt.tricontourf(x, y, z, levels=args.levels, cmap=colour)
cbar = plt.colorbar(bar)
cbar.set_label("Occurrence", rotation=270, labelpad=15)
plt.xlim((xmin, xmax))
plt.xlabel(r"Distance ($\AA$)")
plt.ylabel("Time (ns)")
if args.output:
    plt.savefig(args.output, dpi=300, width=5, height=4, bbox_inches="tight")
else:
    plt.show()
