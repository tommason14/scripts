#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    sys.exit("Syntax: print_shake.py datafile")

""" 
Checking Bond Coeffs section for # HC-CA or # CA-HC etc...
"""
shake = []
found = False
with open(sys.argv[1]) as f:
    for line in f:
        if "Bond Coeffs" in line:
            found = True
            continue
        if "Angle Coeffs" in line:
            break
        if found and len(line.split()) > 1:
            line = line.split()
            if line[-1].startswith('H') or '-H' in line[-1]:
                shake.append(line[0])

print(" ".join(shake))
