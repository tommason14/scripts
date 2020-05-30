#!/usr/bin/env python3

"""
File: interaction_energies_psi4.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Uses the output of chem_assist -r (which writes to energies.csv).
Note splits path on first directory of `Path`. Make sure these are different
for each different system!
"""


import pandas as pd
from dfply import *

# make sure the parent folder of each path is used for each group

@make_symbolic
def if_else(bools, val_if_true, val_if_false):
    return np.where(bools, val_if_true, val_if_false)

df = pd.read_csv('energies.csv')

print((df >>
    mutate(Config = X.Path.str.split('/').str[0]) >>
    mutate(Type = if_else(X.Path.str.contains('frag'), 'frag', 'complex')) >>
    mutate(HF = X['HF/DFT']) >>
    mutate(SRS = X['HF/DFT'] + 1.64 * X.MP2_opp) >>
    mutate(Disp = X.SRS - X.HF) >>
    gather('energy', 'values', [X.HF, X.Disp]) >>
    mutate(energy_type = X.energy + '-' + X.Type) >>
    spread('energy_type', 'values') >>
    group_by(X.Config) >>
    # summing the complexes removes NaN, and need to sum frags anyway for int energy
    summarise(disp_complex = X['Disp-complex'].sum(),
              disp_frags = X['Disp-frag'].sum(),
              hf_complex = X['HF-complex'].sum(),
              hf_frags = X['HF-frag'].sum()) >>
    mutate(hf_int_kj = (X.hf_complex - X.hf_frags) * 2625.5,
           disp_int_kj = (X.disp_complex - X.disp_frags) * 2625.5)
).to_string(index=False, justify='left'))
