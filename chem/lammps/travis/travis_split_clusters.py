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

clusters = []
tmp = []
with open(sys.argv[1]) as f:
    for line in f:
        if re.search('^\s*[0-9]+$', line) and len(tmp) != 0:
            clusters.append(tmp)
            tmp = []
        tmp.append(line)

if not os.path.isdir("clusters"):
    os.mkdir("clusters")

for idx, cluster in enumerate(clusters):
    with open(f"clusters/cluster-{idx}.xyz", "w") as new:
        for line in cluster:
            new.write(line)
