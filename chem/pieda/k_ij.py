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

which you'll get by running chem_assist -d on a set of xyz files to make FMO3 calcs, then
- pieda_make_fmo0.sh
- pieda_make_pl0.sh
- pieda_make_full_pieda.sh
from the base dir
"""

import os
import pandas as pd
import re
import subprocess as sp
from tqdm import tqdm


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
        if (
            'I    J DL  Z    R   Q(I->J)        E"corr          E"uncorr    E"IJ-E"I-E"J,corr/uncorr dDIJ*VIJ,unc   Gsol  tot,corr'
            in line
        ):
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
        if (
            'I   J   K DL   RMIN   RMAX       E"corr      deltaE"IJK,corr/uncorr  dDIJK*VIJK     Gsol     tot'
            in line
        ):
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
        if (
            "I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ    total     Ees      Eex    Ect+mix   Erc+di    Gsol"
            in line
        ):
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


def calc_kij(fmo0=None, fmo3=None, pieda=None):
    CAL_TO_J = 4.184
    HART_TO_KJ = 2625.5
    zero_body = fmo0_energies(fmo0)
    one_body = polarised_monomer_energies(fmo3)
    two_body = two_body_table(fmo3)
    three_body = three_body_table(fmo3)
    pieda = full_pieda_table(pieda)

    one_body_pol = (one_body["Ecorr"].sum() - zero_body["Ecorr"].sum()) * HART_TO_KJ
    two_body_pol = two_body["dDIJ*VIJ,unc"].sum() * HART_TO_KJ
    three_body_pol = three_body["dDIJK*VIJK"].sum() * HART_TO_KJ
    three_body_hf = (
        three_body["corr/uncorr"].sum() - two_body["corr/uncorr"].sum()
    ) * HART_TO_KJ
    three_body_srs_corr = (
        three_body['deltaE"IJK'].sum()
        - three_body["corr/uncorr"].sum()
        - two_body['E"IJ-E"I-E"J'].sum()
    ) * HART_TO_KJ
    two_body_ctmix = pieda["Ect+mix"].sum() * CAL_TO_J
    two_body_rc_di = pieda["Erc+di"].sum() * CAL_TO_J

    disp = one_body_pol + two_body_pol + two_body_ctmix + three_body_pol + three_body_hf
    ind = two_body_rc_di + three_body_srs_corr
    return disp / (disp + ind)


def main():
    dirs = (
        sp.check_output(
            "find . -maxdepth 1 -type d -name 'cluster*' | sort -t'-' -k2 -n",
            shell=True,
        )
        .decode("utf-8")
        .split("\n")[:-1]
    )
    cwd = os.getcwd()
    # collect results for better printing
    res = []
    for d in tqdm(dirs):
        os.chdir(d)
        kij = calc_kij(
            fmo3="spec.log", fmo0="fmo0/spec.log", pieda="full-pieda/spec.log"
        )
        res.append(f"{d.rsplit('/')[1]}: k_ij = {kij:.3f}")
        os.chdir(cwd)
    ave = sum([float(r.split()[-1]) for r in res]) / len(res)
    print('Writing to results.txt')
    with open('results.txt', 'w') as f:
        for r in res:
            f.write(r + '\n')
        f.write(f'Average: {ave:.3f}')


if __name__ == "__main__":
    main()
