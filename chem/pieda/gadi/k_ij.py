#!/usr/bin/env python3
"""
File: k_ij.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Uses FMO3 calculations to extract ratios of
dispersion / dispersion + induction, used to scale lj interactions
of OPLS force fields in order to create a polarisable force field

Expects a file structure of:
.
├── cluster-0
│   ├── cluster-0.xyz
│   ├── fmo0
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── free-state-polarisation
│   │   ├── spec.F40.000
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── full-pieda
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── spec.inp
│   ├── spec.job
│   └── spec.log
├── cluster-1
│   ├── cluster-1.xyz
│   ├── fmo0
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── free-state-polarisation
│   │   ├── spec.F40.000
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── full-pieda
│   │   ├── spec.inp
│   │   ├── spec.job
│   │   └── spec.log
│   ├── spec.inp
│   ├── spec.job
│   └── spec.log

which you'll get by running autochem -d on a set of xyz files to make FMO3 calcs, then
- pieda_make_fmo0.sh
- pieda_make_pl0.sh
- pieda_make_full_pieda.sh
from the base dir

Then run pieda_fmo3_energies.sh (https://github.com/tommason14/scripts/blob/master/chem/pieda/pieda_fmo3_energies.sh)
to write all fmo3_energies to fmo3.csv, and then you can run this script.
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import exp
import re
import subprocess as sp
from tqdm import tqdm
import sys


def eof(filename, percentage):
    """
    return percentage of file from the bottom, with the percentage
    given as a decimal:
    eof(file, 0.2) returns last 20% of file
    """
    with open(filename, "rb") as f:
        f.seek(0, 2)  # Seek @ EOF
        fsize = f.tell()  # Get size
        Dsize = int(percentage * fsize)
        f.seek(max(fsize - Dsize, 0), 0)  # Set pos @ last n chars lines
        lines = f.readlines()  # Read to end

    # RETURN DECODED LINES
    for i in range(len(lines)):
        try:
            lines[i] = lines[i].decode("utf-8")
        except:
            lines[i] = "CORRUPTLINE"
            print("eof function passed a corrupt line in file ", File)
        # FOR LETTER IN SYMBOL
    return lines


def fmo0_energies(filename):
    """
    To calculate the extent of polarisation at the monomer
    level, need monomer energies without the influence of the
    electric field (Coulomb bath)
    """
    cols = ["Frags", "Ecorr"]
    vals = []
    found = False
    for line in eof(filename, 0.2):
        if "Ecorr           Euncorr       DX       DY       DZ" in line:
            found = True
            continue
        if re.match("^\s+$", line):
            found = False
        if found:
            vals.append(line.split()[:2])
    df = pd.DataFrame(vals, columns=cols)
    df["Frags"] = df["Frags"].str.split("(").str[0]
    return df.apply(pd.to_numeric)


def polarised_monomer_energies(filename):
    """
    Return correlated FMO1 energies
    """
    cols = ["Frags", "Ecorr", "Euncorr"]
    vals = []
    found = False
    for line in eof(filename, 0.2):
        if 'E"corr          E"uncorr       DX        DY        DZ' in line:
            found = True
            continue
        if re.match("^\s+$", line):
            found = False
        if found:
            vals.append(line.split()[:3])
    df = pd.DataFrame(vals, columns=cols)
    df["Frags"] = df["Frags"].str.split("(").str[0]
    return df.apply(pd.to_numeric)


def two_body_table(filename):
    """
    Search the last 20% of a file for the table of two-body interactions.
    Returns the table as a pandas dataframe.
    """
    cols = []
    vals = []
    found = False
    for line in eof(filename, 0.2):
        if ('I    J DL  Z    R   Q(I->J)        E"corr          E"uncorr    E"IJ-E"I-E"J,corr/uncorr dDIJ*VIJ,unc   Gsol  tot,corr'
                in line):
            found = True
            for v in line.split():
                # E"IJ-E"I-E"J,corr/uncorr need to be split up
                if v.startswith('E"IJ'):
                    cols += v.split(",")
                else:
                    cols.append(v)
            continue
        if found and re.search("^\s*$", line):
            break
        if found and not re.search("^\s*-+$", line):
            vals.append(line.split())
    df = pd.DataFrame(vals, columns=cols)
    conv = df.columns.drop("DL")
    df[conv] = df[conv].apply(pd.to_numeric)
    return df


def three_body_table(filename):
    """
    Search the last 20% of a file for the table of three-body interactions.
    Returns the table as a pandas dataframe.
    """
    cols = []
    vals = []
    found = False
    for line in eof(filename, 0.2):
        if ('I   J   K DL   RMIN   RMAX       E"corr      deltaE"IJK,corr/uncorr  dDIJK*VIJK     Gsol     tot'
                in line):
            found = True
            for v in line.split():
                # deltaE"IJK,corr/uncorr need to be split up
                if v.startswith("delta"):
                    cols += v.split(",")
                else:
                    cols.append(v)
            continue
        if found and re.search("^\s*$", line):
            break
        if found and not re.search("^\s*-+$", line):
            vals.append(line.split())
    df = pd.DataFrame(vals, columns=cols)
    conv = df.columns.drop("DL")
    df[conv] = df[conv].apply(pd.to_numeric)
    return df


def full_pieda_table(piedalog):
    """
    Searches the last 10% of the file for a table of two-body decomposed interactions.
    Returns a pandas dataframe
    """
    found = False
    cols = None
    vals = []
    for line in eof(piedalog, 0.1):
        if ("I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ    total     Ees      Eex    Ect+mix   Erc+di    Gsol"
                in line):
            found = True
            cols = line.split()
            continue
        if found and re.search("^\s*$", line):
            break
        if found and not re.search("^\s*-+$", line):
            vals.append(line.split())
    df = pd.DataFrame(vals, columns=cols)
    conv = df.columns.drop("DL")
    df[conv] = df[conv].apply(pd.to_numeric)
    return df


def pol_from_pieda(pieda_log):
    with open(pieda_log) as f:
        for line in f:
            if 'Polarisation (total)' in line:
                return float(line.split()[-1])


def calc_kij(fmo0=None, fmo3=None, pieda=None):
    CAL_TO_J = 4.184
    HART_TO_KJ = 2625.5
    three_body = three_body_table(fmo3)
    pieda_df = full_pieda_table(pieda)

    pieda_pol = pol_from_pieda(pieda)
    three_body_pol = three_body["dDIJK*VIJK"].sum() * HART_TO_KJ
    three_body_hf = three_body["corr/uncorr"].sum() * HART_TO_KJ
    three_body_srs_corr = (three_body['deltaE"IJK'].sum() -
                           three_body["corr/uncorr"].sum()) * HART_TO_KJ
    two_body_ctmix = pieda_df["Ect+mix"].sum() * CAL_TO_J
    two_body_rc_di = pieda_df["Erc+di"].sum() * CAL_TO_J

    indPI = pieda_pol + two_body_ctmix
    indFMO = three_body_hf + three_body_pol
    disp = two_body_rc_di + three_body_srs_corr
    ind = indPI + indFMO
    return (disp / (disp + ind))

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

def check_files():
    if 'fmo3.csv' not in os.listdir('.'):
        print('Run pieda_fmo3_energies.sh '
        '(https://github.com/tommason14/scripts/blob/master/chem/pieda/pieda_fmo3_energies.sh)\n'
        'from the top-level directory to generate fmo3.csv, then try again')
        sys.exit(1)

def main(show_fig=False):
    check_files()
    dirs = (sp.check_output(
        "find . -maxdepth 1 -type d | tail -n +2 | sort",
        shell=True,
    ).decode("utf-8").split("\n")[:-1])
    cwd = os.getcwd()
    folders = []
    res = []
    for d in tqdm(dirs):
        folders.append(d)
        os.chdir(d)
        kij = calc_kij(fmo3="spec.log",
                       fmo0="fmo0/spec.log",
                       pieda="full-pieda/spec.log")
        res.append(kij)
        os.chdir(cwd)
    kij = pd.DataFrame({'Path': folders, 'k_ij': res})

    H_to_kJ = 2625.5
    energies = pd.read_csv('fmo3.csv')
    df = energies.merge(kij, on='Path')
    df['diffs'] = (df['MP2/SRS'] - df['MP2/SRS'].min()) * H_to_kJ
    plot_graph(df, show_fig)
    df['weights'] = boltz(df['diffs'])

    with open('results.txt', 'w') as f:
        for d, r in zip(folders, res):
            f.write(f'{d.replace("./", "")}: {r:.3f}\n')
        f.write(f'Mean value = {df.k_ij.mean():.3f}\n')
        f.write(f'Boltzmann averaged value = {(df.weights * df.k_ij).sum():.3f}')
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        show_fig = sys.argv[1] == '--show'
    else:
        show_fig = False
    main(show_fig)
