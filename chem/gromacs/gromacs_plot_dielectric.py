#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

if not os.path.isfile('epsilon.xvg'):
    sys.exit(
        'Error: assumes that a file named epsilon.xvg exists in the current directory.'
    )

df = pd.read_csv(
    "epsilon.xvg",
    sep="\s+",
    engine="python",
    skiprows=25,
    names=["time", "eps", "kirkwood_G", "kirkwood_g"],
)
df["time"] /= 1000  # ps to ns

sns.set(style="whitegrid", font_scale=1.5)  # font size for keynote
sns.lineplot(x="time", y="eps", data=df, ci=None)
plt.xlabel("Time (ns)")
plt.ylabel("Dielectric Constant")
plt.savefig("dielectric.png", dpi=300, bbox_inches="tight")
