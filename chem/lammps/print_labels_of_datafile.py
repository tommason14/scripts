#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    sys.exit("Syntax: print_labels_of_datafile.py datafile")

"""
Extracting labels from Masses section:
1      12.0100    # C
2      12.0100    # C3
3       1.0080    # HC
4       1.0080    # H1
5       1.0080    # HO
6       1.0080    # HN
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
print(" ".join(labels))
