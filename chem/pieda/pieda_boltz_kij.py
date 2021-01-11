#!/usr/bin/env python3

from numpy import exp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

cwd = os.path.dirname(os.path.abspath(__file__))

if len(sys.argv) == 2 and sys.argv[1] == '-h':
    print(f"""Produces a Boltzmann-averaged k_ij value, weighted by the total energy of each system.
Run {cwd}/fmo3_energies.sh beforehand to generate a csv with FMO3
energies (fmo3.csv) that is used by this script.
A plot of total energy against k_ij is saved as kij_vs_fmo3.png,
but to view the plot in a popup window, type:
$ pieda_boltz_kij.py --show""")
    sys.exit(1)


def k_ij_values(path):
    folders = []
    values = []
    with open(path) as f:
        for line in f:
            if 'k_ij' in line:
                folder, value = line.strip().split(': k_ij = ')
                folders.append(folder)
                values.append(float(value))
    return pd.DataFrame({'Path': folders, 'k_ij': values})


def boltz(diffs):
    """
    Pass in energy differences in kJ/mol
    """
    R = 8.3145
    T = 298.15
    kJ_to_J = 1000
    numerator = exp((-1 * kJ_to_J * diffs) / (R * T))
    return numerator / numerator.sum()

def plot_graph(df, show_fig = False):
    sns.set(style='white')
    plt.rcParams['mathtext.default'] = 'regular'
    sns.scatterplot(x='diffs', y='k_ij', data=df)
    plt.xlabel(r'$\Delta$E$_{Tot}$ (kJ mol$^{-1}$)')
    plt.ylabel(r'k$_{ij}$')
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.savefig('kij_vs_fmo3.png', dpi=300, bbox_inches='tight')
    print('Plot of kij against energy saved as kij_vs_fmo3.png')
    if show_fig:
        plt.show()

def main(show_fig=False):
    H_to_kJ = 2625.5
    kij = k_ij_values('results.txt')
    energies = pd.read_csv('fmo3.csv')
    energies.Path = energies.Path.str.replace('\./', '')
    df = energies.merge(kij, on='Path')
    df['diffs'] = (df['MP2/SRS'] - df['MP2/SRS'].min()) * H_to_kJ
    plot_graph(df, show_fig)
    df['weights'] = boltz(df['diffs'])
    print(f'Boltzmann averaged value = {(df.weights * df.k_ij).sum():.3f}')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        show_fig = sys.argv[1] == '--show'
    else:
        show_fig = False
    main(show_fig)
