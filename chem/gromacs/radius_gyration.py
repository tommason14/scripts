#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from tqdm import tqdm

if len(sys.argv) < 3 or any("-h" in x for x in sys.argv[1:]):
    print("Syntax: radius_of_gyration.py file.tpr file.xtc [file2.xtc...]")
    sys.exit(1)

u = mda.Universe(sys.argv[1], *sys.argv[2:])

pol = u.select_atoms("resname pol")

# Rg for multiple polymers
resids = np.unique(pol.resids)
rg = []
for ts in tqdm(u.trajectory):
    for res in resids:
        polymer = pol.select_atoms(f"resid {res}")
        rg.append((res, ts.time, polymer.radius_of_gyration()))
df = pd.DataFrame(rg, columns=["Residue", "Time (ps)", "Rg"])
df["Time (ns)"] = df["Time (ps)"] / 1000
df.drop("Time (ps)", axis=1, inplace=True)
df.to_csv("rg.csv", index=False)

sns.set(style="white", font="Helvetica")
p = sns.relplot(x="Time (ns)", y="Rg", hue="Residue", kind="line", data=df)
p.set_axis_labels("Time (ns)", r"R$_{gyr}$ (${\AA}$)")
plt.savefig("rg.png", dpi=300, bbox_inches="tight")
