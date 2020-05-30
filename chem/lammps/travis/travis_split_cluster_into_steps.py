#!/usr/bin/env python3

"""
Takes in the output from the travis program after a 'cut'
operation has been performed, which generates one xyz file
containing all cluster geometries. This script takes those
geometries and separates them into individual files in a 
'cluster' subfolder.
"""

import os
import re
import sys

if len(sys.argv) != 2:
    sys.exit(f"Syntax: {os.path.basename(__file__)} <cluster.xyz>")

files = {}

with open(sys.argv[1]) as f:
    steps = [re.search("Step ([0-9]+)", line).group(1) for line in f if "#" in line]

clusters = []
tmp = []
with open(sys.argv[1]) as f:
    for line in f:
        if re.search("^\s*[0-9]+$", line) and len(tmp) != 0:
            clusters.append(tmp)
            tmp = []
        tmp.append(line)

if not os.path.isdir("clusters"):
    os.mkdir("clusters")

for step, cluster in zip(steps, clusters):
    with open(f"clusters/step-{step}.xyz", "w") as new:
        for line in cluster:
            new.write(line)
