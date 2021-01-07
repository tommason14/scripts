#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    sys.exit("Syntax: lammps_print_elements.py datafile")

atomic_wt = ['H',  'Li', 'B',  'C',                  
             'N',  'O',  'F',  'Ne',                
             'Na', 'Mg', 'Al', 'Si',            
             'P',  'S',  'Cl', 'Ar',               
             'K',  'Ca', 'Ti', 'Fe',
             'Zn', 'Se', 'Br', 'Kr',              
             'Mo', 'Ru', 'Sn', 'Te',
             'I',  'Xe', 'D',  'M']                                       

def atomic_symbol(name):
    if name[:2] in atomic_wt:                            
        return name[:2]                                  
    elif name[0] in atomic_wt:                           
        return name[0]                                   
    else:
        return name                                      

types = []
found = False
with open(sys.argv[1]) as f:
    for line in f:
        if "Masses" in line:
            found = True
            continue
        if len(types) > 0 and re.search("^\s*[A-Z]", line):
            break
        if found and len(line.split()) > 1:
            if line.strip().endswith('DC'):
                types.append(line.split()[-2])
            elif line.strip().endswith('DP'):
                types.append('D')
            else:
                types.append(line.split()[-1])

types = [atomic_symbol(name) for name in types]
print(" ".join(types))
