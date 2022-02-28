#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import re

if len(sys.argv) != 2 or "-h" in sys.argv:
    print("Plot potential energy and density over time")
    print(f"Syntax: {os.path.basename(__file__)} omm.log")
    sys.exit(1)

fname = sys.argv[1]
data = []
header = None
found = False
with open(fname) as f:
    for line in f.readlines():
        if '#"Step' in line:
            header = line.strip("#").replace('"', "").split(",")
            found = True
            continue
        if found and re.search("^[0-9]", line):  # also check if error
            data.append(line.split(","))
df = (
    pd.DataFrame(data, columns=header)
    .apply(lambda col: pd.to_numeric(col))
    .melt(
        id_vars="Time (ps)", value_vars=["Potential Energy (kJ/mole)", "Density (g/mL)"]
    )
)
sns.set(style="ticks", font="DejaVu Sans")

p = sns.relplot(
    x="Time (ps)",
    y="value",
    col="variable",
    kind="line",
    ci=None,
    data=df,
    facet_kws={"sharex": False, "sharey": False},
)
p.set_titles("{col_name}")
p.set_axis_labels("", "")
plt.show()
