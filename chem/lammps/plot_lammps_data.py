#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

df = pd.read_csv(sys.stdin)
df = df.melt(id_vars='Step', var_name='data')
g = sns.relplot(x='Step',
                y='value',
                col='data',
                col_wrap=4,
                data=df,
                kind='line',
                ci=None,
                facet_kws={'sharey': False})
g.set_titles('{col_name}')
plt.show()
