#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-b", "--box", help="Box length in angstroms", type=float, required=True
)
parser.add_argument(
    "-d",
    "--diffusion",
    help="Diffusion coefficient from the simulation in cm2/s",
    type=float,
    required=True,
)
parser.add_argument(
    "-v",
    "--viscosity",
    help="Shear viscosity of the solvent, in mPas. Default = 0.896 (assuming water, DOI:10.1063/1.3330544)",
    default=0.896,
    type=float,
)
parser.add_argument(
    "-t",
    "--temperature",
    help="Temperature of the simulation in K. Default = 300",
    default=300,
    type=float,
)

args = parser.parse_args()
"""
Following procedure laid out in DOI:10.1021/jp0477147, describing the correction that should be
applied to diffusion coefficients found when periodic boundary conditions are applied
"""

from scipy.constants import pi, k

CONS = 2.837297
MPAS_TO_PAS = 1e-3
ANGSTROMS_TO_M = 1e-10
SQUARE_METRES_TO_SQUARE_CM = 1e4
correction = (
    (k * args.temperature * CONS)
    / (6 * pi * args.viscosity * MPAS_TO_PAS * args.box * ANGSTROMS_TO_M)
) * SQUARE_METRES_TO_SQUARE_CM

corrected = args.diffusion + correction

print(f"Diffusion coefficient from simulation    = {args.diffusion} cm²/s")
print(f"Diffusion coefficient corrected for PBCs = {corrected} cm²/s")
