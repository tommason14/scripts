from polfuncs import PolSim
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

# run from the "root" folder:
# .
# ├── anneal_1
# ├── anneal_2
# ├── anneal_3
# ├── anneal_4
# ├── anneal_5
# ├── em
# ├── npt_1
# ├── npt_2
# ├── npt_3
# ├── npt_4
# ├── npt_5
# └── pack
# scripts expects npt_X folders one level deep

dirs = [x for x in os.listdir() if os.path.isdir(x) and "npt" in x]

data = []

# should keep z_coord in case z-dimension changes

for d in tqdm(dirs):
    os.chdir(d)
    cycle = int(d.split("_")[1])
    data.append(
        PolSim("npt.tpr", "npt.xtc")
        .compute_partial_densities(molecules=["pol"])
        .assign(cycle=cycle)
    )
    os.chdir("..")

data = pd.concat(data, ignore_index=True)

print("Plotting...")

sns.set(style="ticks", font="DejaVu Sans")
sns.lineplot(data=data, x="z_coord", y="pol", hue="cycle", ci=None).set(
    xlabel="Z-coordinate (Å)", ylabel="Density (g cm⁻³)"
)
plt.savefig("polymer_density.png", dpi=300)
