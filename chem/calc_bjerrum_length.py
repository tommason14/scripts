#!/usr/bin/env python3
from scipy.constants import pi, e, k, epsilon_0
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-s",
    "--sep",
    help="Separation between charged groups in Angstroms",
    required=True,
    type=float,
)
parser.add_argument(
    "-d",
    "--dielectric",
    help="Dielectric constant of polymer",
    required=True,
    type=float,
)
parser.add_argument(
    "-t",
    "--temp",
    help="Desired temperature in Kelvin, default = 300",
    default=300,
    type=float,
)

args = parser.parse_args()


def bjerrum_length(eps, b, T=300):
    b = b / 1e10  # angstroms to m
    return e ** 2 / (4 * pi * epsilon_0 * eps * k * T * b)


bjerrum = bjerrum_length(args.dielectric, args.sep, args.temp)
print(f"Bjerrum length = {bjerrum:.2f} Å")
