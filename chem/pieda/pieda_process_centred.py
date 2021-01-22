#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from funcs import boltz
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Takes raw_data.csv produced by k_ij.py '
    '(https://github.com/tommason14/scripts/blob/master/chem/pieda/gadi/k_ij.py) '
    'and plots k_ij against the difference in total energy of each cation or anion '
    'centred configuration. This script also finds mean and Boltzmann-weighted average '
    'k_ij parameters for each centred ion')
parser.add_argument('-c', '--csv', help='csv file produced by k_ij.py', default='raw_data.csv',
        required=True)
parser.add_argument('-p', '--png', help='filename of png created', default='k_ij.png')
parser.add_argument('-o', '--output', help='filename of csv containing mean/BW k_ij values',
        default='average_results.csv')
args = parser.parse_args()

H_to_kJ = 2625.5

df = pd.read_csv(args.csv)
df['centred'] = df['Path'].str.split('/').str[1].str.split(
    '-').str[:2].str.join(' ').str.capitalize()
df['diffs'] = df.groupby('centred')['MP2/SRS'].transform(
    lambda x: (x - x.min()) * H_to_kJ)
df['weights'] = df.groupby('centred')['diffs'].transform(lambda x: boltz(x))
sns.set(style='white')
plt.rcParams['mathtext.default'] = 'regular'
sns.lmplot(x='diffs',
           y='k_ij',
           hue='centred',
           palette='Dark2',
           data=df,
           legend=False,
           ci=None)
plt.xlabel(r'$\Delta$E$_{Tot}$ (kJ mol$^{-1}$)')
plt.ylabel(r'k$_{ij}$')
plt.xlim(df['diffs'].min() - 10,
         df['diffs'].max() + 10)  # extend x axis so that points aren't cut off
plt.legend()
plt.savefig(args.png, dpi=300, bbox_inches='tight')

bw = df.groupby('centred').apply(lambda g: (g['weights'] * g['k_ij']).sum()
                                 ).reset_index(name='boltz_weighted')
means = df.groupby('centred')['k_ij'].mean().reset_index(name='mean')
merged_means = bw.merge(means, on='centred')
merged_means = merged_means.rename(columns={'centred': 'grouping'})
no_groupby = df.drop(['diffs', 'weights'], axis=1).copy()
no_groupby['diffs'] = (no_groupby['MP2/SRS'] -
                       no_groupby['MP2/SRS'].min()) * H_to_kJ
no_groupby['weights'] = boltz(no_groupby['diffs'])

bw_ = (no_groupby['weights'] * no_groupby['k_ij']).sum()
mean_ = no_groupby['k_ij'].mean()
merged_means.append(
    {
        'grouping': 'all configs',
        'boltz_weighted': bw_,
        'mean': mean_
    },
    ignore_index=True).to_csv(args.output, index=False)
