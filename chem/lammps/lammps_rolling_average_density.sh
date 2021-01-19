#!/usr/bin/env bash

[[ $# -ne 2 || $1 == '-h' ]] &&
echo "Syntax: $(basename $0) num_steps filename" &&
exit 1

grep_lammps_data.sh $2 | csvcut -c Step,Density | 
python3 -c "
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
f = lambda x, pos: f'{x/10**6:,.1f}x10$^6$'
df = pd.read_csv(sys.stdin)
df['rolling'] = df['Density'].rolling($1).mean()
rolling = df['rolling'].values[-1]
fig, ax = plt.subplots()
ax.plot(df['Step'], df['Density'], 'k.')
ax.plot(df['Step'], df['rolling'], 'r')
ax.xaxis.set_major_formatter(FuncFormatter(f))
ax.set_xlabel('Timesteps')
ax.set_ylabel('Density')
ax.set_title(f'$1 step average density = {rolling:.3f}')
plt.show()"
