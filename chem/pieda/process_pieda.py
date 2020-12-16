#!/usr/bin/env python3

from numpy import exp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

cwd = os.path.dirname(os.path.abspath(__file__))

if len(sys.argv) == 2 and sys.argv[1] == '-h':
    print(f"""Plot total energies against k_ij by typing:
$ process_pieda.py --plot
Otherwise a Boltzmann-averaged k_ij value is produced, weighted by the total energy of each system.
Run {cwd}/fmo3_energies.sh beforehand to generate a csv with FMO3
energies (fmo3.csv) that is used by this script.""")
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


def plot_graph(df):
    sns.set(style='white', font='Helvetica')
    sns.scatterplot(x='diffs', y='k_ij', data=df)
    plt.xlabel(r'$\Delta$E$_{Tot}$ (kJ mol$^{-1}$)')
    plt.ylabel(r'k$_{ij}$')
    plt.savefig('energies_v_kij.png',
                dpi=300,
                width=5,
                height=4,
                bbox_inches='tight')
    plt.show()
    print('saved image as energies_v_kij.png')


def main(plot=False):
    H_to_kJ = 2625.5
    kij = k_ij_values('results.txt')
    energies = pd.read_csv('fmo3.csv')
    energies.Path = energies.Path.str.replace('\./', '')
    df = energies.merge(kij, on='Path')
    df['diffs'] = (df['MP2/SRS'] - df['MP2/SRS'].min()) * H_to_kJ
    if plot:
        plot_graph(df)
    else:
        df['weights'] = boltz(df['diffs'])
        print(f'Boltzmann averaged value = {(df.weights * df.k_ij).sum():.3f}')


if __name__ == "__main__":
    try:
        plot = sys.argv[1] == '--plot'
    except IndexError:
        plot = False
    main(plot)
