#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    sys.exit("Syntax: lammps_print_drude_fix.py datafile")

"""
Extracting labels from Masses section:
1      12.0100    # C DC
2      12.0100    # C DC
3       1.0080    # H
4       1.0080    # H
5       1.0080    # H
6       0.4       # C DP
"""
labels = []
found = False
with open(sys.argv[1]) as f:
    for line in f:
        if "Masses" in line:
            found = True
            continue
        if len(labels) > 0 and re.search("^\s*[A-Z]", line):
            break
        if found and len(line.split()) > 1:
            labels.append(line.split()[-1])

def convert(label):
    if label == 'DC':
        return 'C'
    elif label == 'DP':
        return 'D'
    else: # hydrogens
        return 'N'

labels = map(convert, labels)
print(" ".join(labels))
