#!/usr/bin/env python3

import MDAnalysis as mda
import argparse
from utils import completion

parser = argparse.ArgumentParser(
    description="Concatenate multiple trajectory files into one file using MDAnalysis"
)
parser.add_argument("files", help="Coordinate and trajectory files", nargs="+")
parser.add_argument(
    "-o", "--output", help="Filename of whole trajectory to be created", required=True
)
args = parser.parse_args()

u = mda.Universe(*args.files)
with mda.Writer(args.output, len(u.atoms)) as f:
    for frame in completion(u.trajectory):
        f.write(u)
