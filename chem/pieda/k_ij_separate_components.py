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
"""

import os
import pandas as pd
from numpy import exp
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


def fmo_total_energy(log):
    """
    Reads the FMO energy block and extracts the final total energy reported with
    scaled MP2 parameters
    """
    energy = 0.0
    for line in eof(log, 0.1):
        if "Total SCS" in line:
            energy = float(line.split()[-1])
    return energy


def total_pol_from_pieda(pieda_log):
    with open(pieda_log) as f:
        for line in f:
            if "Polarisation (total)" in line:
                return float(line.split()[-1])


def k_ij_value(pieda_df, three_body, subset=None, total_pieda_pol=None):
    """
    pieda_df = Pandas dataframe from the full pieda calculation
    three_body = Pandas dataframe from the FMO3 calculation
    subset = list of indices of cations/anions. If given, only the terms
    where all components contain those indices will be used. I.e. for two body terms
    of anions ([1,2,3,4]), both the i and j fragments will be fragment 1, 2, 3 or 4.
    For all interactions (cation-anion k_ij), pass in subset='all'
    """
    CAL_TO_J = 4.184
    HART_TO_KJ = 2625.5
    if subset == "all":
        three_rows = three_body.copy()  # all rows
        pieda_rows = pieda_df.copy()
    else:
        three_rows = three_body[
            (three_body["I"].isin(subset))
            & (three_body["J"].isin(subset))
            & (three_body["K"].isin(subset))
        ]
        pieda_rows = pieda_df[
            (pieda_df["I"].isin(subset)) & (pieda_df["J"].isin(subset))
        ]

    three_body_pol = three_rows["dDIJK*VIJK"].sum() * HART_TO_KJ
    three_body_hf = three_rows["corr/uncorr"].sum() * HART_TO_KJ
    three_body_srs_corr = (
        three_rows['deltaE"IJK'].sum() - three_rows["corr/uncorr"].sum()
    ) * HART_TO_KJ
    two_body_ctmix = pieda_rows["Ect+mix"].sum() * CAL_TO_J
    two_body_rc_di = pieda_rows["Erc+di"].sum() * CAL_TO_J
    if total_pieda_pol is not None:
        indPI = total_pieda_pol + two_body_ctmix  # cation-anion
    else:
        indPI = two_body_ctmix  # rest

    indFMO = three_body_hf + three_body_pol
    disp = two_body_rc_di + three_body_srs_corr
    ind = indPI + indFMO
    return disp / (disp + ind)


def get_ions(inp):
    """
    Returns a list of indices of cations and anions in the calculation.
    Note that no check for neutral fragments is performed - the function
    assumes that each fragment is either a singly-charged cation
    or singly-charged anion.
    """
    charges = []
    with open(inp) as f:
        for line in f:
            if "ICHARG(1)" in line:
                line = line.strip().split("=")[1]
                charges = [int(i) for i in line.split(",")]
                break
    if len(charges) == 0:
        raise AttributeError(f"Failed to find charges of each ion in {inp}")
    anions = []
    cations = []
    neutrals = []
    for ind, charge in enumerate(charges, 1):  # gamess says first frag is 1, not 0
        if charge == -1:
            anions.append(ind)
        elif charge == 1:
            cations.append(ind)
        elif charge == 0:
            neutrals.append(ind)
        else:
            raise AttributeError(
                "Multivalent ions deliberately not accounted for - sometimes multiple cations/anions "
                "of the same charge may be included in the same calculation and the current code will "
                "not account for this."
            )
    return cations, anions, neutrals


def calc_kij(fmo0=None, fmo3=None, pieda=None, inp=None):
    three_body = three_body_table(fmo3)
    pieda_df = full_pieda_table(pieda)
    cations, anions, neutrals = get_ions(inp)
    total_pieda_pol = total_pol_from_pieda(pieda)

    # deal with different parameters by returning values with a dict
    res = {
        "Cation-Anion k_ij": None,
        "Cation-Cation k_ij": None,
        "Anion-Anion k_ij": None,
        "Cation-Neutral k_ij": None,
        "Anion-Neutral k_ij": None,
        "Neutral-Neutral k_ij": None,
    }

    res["Cation-Cation k_ij"] = k_ij_value(pieda_df, three_body, subset=cations)
    res["Anion-Anion k_ij"] = k_ij_value(pieda_df, three_body, subset=anions)
    if len(neutrals) == 0:
        res["Cation-Anion k_ij"] = k_ij_value(
            pieda_df, three_body, subset="all", total_pieda_pol=total_pieda_pol
        )
    else:
        # Can't include total polarisation as the Coulomb bath will also include neutral terms
        res["Cation-Anion k_ij"] = k_ij_value(
            pieda_df, three_body, subset=cations + anions
        )
        res["Cation-Neutral k_ij"] = k_ij_value(
            pieda_df, three_body, subset=cations + neutrals
        )
        res["Anion-Neutral k_ij"] = k_ij_value(
            pieda_df, three_body, subset=anions + neutrals
        )
        res["Neutral-Neutral k_ij"] = k_ij_value(pieda_df, three_body, subset=neutrals)
    return res


def obtain_results(dirs):
    """
    Run calc_kij on each subdirectory passed in, returning a Pandas
    DataFrame with the results. 
    Also includes the total energy in kJ/mol (from the FMO3 calculation), 
    in order to Boltzmann-weight the k_ij values.
    """
    results = {
        "Configuration": [],
        "Cation-Anion k_ij": [],
        "Cation-Cation k_ij": [],
        "Anion-Anion k_ij": [],
        "Cation-Neutral k_ij": [],
        "Anion-Neutral k_ij": [],
        "Neutral-Neutral k_ij": [],
        "Total Energy": [],
    }
    cwd = os.getcwd()
    for d in tqdm(dirs):
        os.chdir(d)
        files = {
            "fmo3": "spec.log",
            "fmo0": "fmo0/spec.log",
            "pieda": "full-pieda/spec.log",
            "inp": "spec.inp",
        }
        if all(os.path.isfile(x) for x in files.values()):
            res = calc_kij(**files)
            results["Total Energy"].append(fmo_total_energy("spec.log"))
            results["Configuration"].append(d.replace("./", ""))
            for key, val in res.items():
                results[key].append(val)
        os.chdir(cwd)
    # drop columns if all NA - will happen if no neutrals present
    return pd.DataFrame(results).dropna(axis=1, how="all")


def _boltz(energies):
    """
    Pass in energies in Hartrees (just the raw total energies)
    """
    HART_TO_KJ = 2625.5
    R = 8.3145
    T = 298.15
    kJ_to_J = 1000
    diffs = (energies - energies.min()) * HART_TO_KJ
    numerator = exp((-1 * kJ_to_J * diffs) / (R * T))
    return numerator / numerator.sum()


def compute_boltzmann_weighted_values(df):
    """
    Pass in a DataFrame created with 'obtain_results()'
    """
    df["Weights"] = _boltz(df["Total Energy"])
    df_with_weights = df.copy()
    df = df.drop("Total Energy", axis="columns")
    df = df.melt(
        id_vars=["Configuration", "Weights"], var_name="Type", value_name="kij"
    )
    df = (
        df.groupby("Type")[["Weights", "kij"]]
        .apply(lambda x: x.prod(axis=1).sum())
        .reset_index()
    )
    return df_with_weights, df


def save_all_kij(df):
    """
    Pass in a DataFrame created with 'obtain_results()'
    """
    print("All results saved to all_k_ij.csv")
    df.to_csv("all_k_ij.csv", index=False)


def save_mean_kij(df):
    """
    Pass in a DataFrame created with 'obtain_results()'
    """
    print("Averages saved to average_k_ij.csv")
    cols = [
        "Cation-Anion k_ij",
        "Cation-Cation k_ij",
        "Anion-Anion k_ij",
        "Cation-Neutral k_ij",
        "Anion-Neutral k_ij",
        "Neutral-Neutral k_ij",
    ]

    cols_to_agg = [col for col in cols if col in df.columns]
    ave = df.agg({col: "mean" for col in cols_to_agg})
    ave.to_csv("average_k_ij.csv", header=False)


def save_bw_kij(df):
    """
    Pass in a DataFrame created with 'compute_boltzmann_weighted_values()'
    """
    print("Boltzmann-weighted values saved to bw_k_ij.csv")
    df.to_csv("bw_k_ij.csv", index=False, header=False)


def main():
    dirs = (
        sp.check_output("find . -type d -name '*cluster*' | sort", shell=True,)
        .decode("utf-8")
        .split("\n")[:-1]
    )
    results = obtain_results(dirs)
    results_with_weights, bw = compute_boltzmann_weighted_values(results)
    save_all_kij(results_with_weights)
    save_mean_kij(results)
    save_bw_kij(bw)


if __name__ == "__main__":
    main()
